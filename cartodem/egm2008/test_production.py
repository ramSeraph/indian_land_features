"""Small end-to-end checks of mosaics, categorical coverage and packaged tiles."""
import io
import json
import shutil
import sqlite3
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

import numpy as np
import rasterio
from PIL import Image
from rasterio.transform import from_origin
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.warp import transform
from rasterio.windows import Window

from build import atomic_json,create_mosaic,setup,download
from production_tiles import Sampler,tiles,package,validate_rgba


class ProductionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=setup(Path(self.temp.name)/'build')
        self.a=np.full((8,8),10,dtype=np.float32)
        self.b=np.full((8,8),20,dtype=np.float32)
        self.a[2:5,2:5]=-255
        paths=[]
        for i,data in enumerate((self.a,self.b)):
            p=self.root/'sheets/coastal'/f'{i}.tif'
            with rasterio.open(p,'w',driver='GTiff',count=1,dtype='float32',
                               width=8,height=8,nodata=-255,crs='EPSG:4326',
                               transform=from_origin(75+i*8/3600,13,1/3600,1/3600)) as dst:
                dst.write(data,1)
            paths.append(p)
        self.height=self.root/'mosaics/pilot.vrt'
        self.mask=self.root/'mosaics/pilot-validity.vrt'
        create_mosaic(paths,self.height,self.mask)

    def test_small_download_is_flushed_before_size_check(self):
        payload=b'small complete download'
        response=io.BytesIO(payload)
        response.headers={'Content-Length':str(len(payload)),'ETag':'test'}
        path=self.root/'inputs/small.bin'
        with patch('urllib.request.urlopen',return_value=response) as fetch:
            result=download('https://example.invalid/input',path,expected_size=len(payload))
            self.assertEqual(fetch.call_count,1)
        self.assertEqual(path.read_bytes(),payload)
        self.assertEqual(result['bytes'],len(payload))
        with patch('urllib.request.urlopen') as fetch:
            download('https://example.invalid/input',path,expected_sha=result['sha256'])
            fetch.assert_not_called()

    def test_portable_mosaic_keeps_validity_zero_as_data(self):
        expected=np.concatenate((self.a,self.b),axis=1)
        with rasterio.open(self.height) as h,rasterio.open(self.mask) as m:
            np.testing.assert_array_equal(h.read(1),expected)
            np.testing.assert_array_equal(m.read(1),np.where(expected==-255,0,255))
            self.assertIsNone(m.nodata)
        moved=self.root.with_name('moved')
        shutil.move(self.root,moved)
        with rasterio.open(moved/'mosaics/pilot.vrt') as h:
            np.testing.assert_array_equal(h.read(1),expected)

    def test_resampling_preserves_pixel_centre_holes_and_excludes_sentinel(self):
        expected=np.concatenate((self.a,self.b),axis=1)
        valid_count=invalid_count=0
        with Sampler(self.height,self.mask,12) as sampler:
            g=sampler.grid
            for x in range(g['xmin'],g['xmax']+1):
                for y in range(g['ymin'],g['ymax']+1):
                    h,valid=sampler.read(x,y)
                    valid_count+=int(valid.sum())
                    invalid_count+=int((~valid).sum())
                    self.assertTrue(np.all(h[valid]>=10-1e-5))
                    self.assertTrue(np.all(h[valid]<=20+1e-5))
                    rows,cols=np.nonzero(valid)
                    mx,my=rasterio.transform.xy(g['transform'],
                            rows+(y-g['ymin'])*512,cols+(x-g['xmin'])*512)
                    lon,lat=transform('EPSG:3857','EPSG:4326',mx,my)
                    rr=np.floor((13-np.array(lat))*3600).astype(int)
                    cc=np.floor((np.array(lon)-75)*3600).astype(int)
                    self.assertTrue(np.all((rr>=0)&(rr<8)&(cc>=0)&(cc<16)))
                    self.assertTrue(np.all(expected[rr,cc]!=-255))
        self.assertGreater(valid_count,0)
        self.assertGreater(invalid_count,0)

    def test_tiles_resume_and_pmtiles_keep_exact_rgba(self):
        first=tiles(self.root,'pilot',11,12)
        again=tiles(self.root,'pilot',11,12)
        self.assertEqual(first[0]['tiles'],again[0]['tiles'])
        output=package(self.root,'pilot',max_bytes=1000000)
        self.assertEqual(sum(x['verified_tiles'] for x in output),sum(z['tiles'] for z in first))
        with sqlite3.connect(self.root/'tiles/pilot/z12.sqlite') as db,Sampler(self.height,self.mask,12) as sampler:
            for x,y,data,valid in db.execute('SELECT x,y,data,valid FROM tiles'):
                rgba=np.asarray(Image.open(io.BytesIO(data)).convert('RGBA'))
                _,coverage=sampler.read(x,y)
                self.assertEqual(validate_rgba(rgba,12,coverage),valid)
        recipe=json.loads((self.root/'tiles/pilot/recipe.json').read_text())
        self.assertNotIn(str(self.root),json.dumps(recipe))

    def test_chunked_tiles_match_whole_grid_warp_across_sheet_and_tile_edges(self):
        # A native 30 m DEM straddles an XYZ boundary and a source-sheet join.
        rows,cols=np.indices((96,192))
        height=(100+rows*.7+cols*.3+80*np.sin(cols/7)).astype('float32')
        height[25:45,85:110]=-255
        paths=[]
        for i in range(2):
            path=self.root/'sheets/coastal'/f'edge{i}.tif'
            with rasterio.open(path,'w',driver='GTiff',count=1,dtype='float32',
                               width=96,height=96,nodata=-255,crs='EPSG:4326',
                               transform=from_origin(75.234375+(i-1)*96/3600,13,1/3600,1/3600)) as dst:
                dst.write(height[:,i*96:(i+1)*96],1)
            paths.append(path)
        create_mosaic(paths,self.height,self.mask)
        for zoom in (9,12):
            with self.subTest(zoom=zoom),Sampler(self.height,self.mask,zoom) as sampler:
                g=sampler.grid;count=0
                common=dict(src_crs='EPSG:4326',crs='EPSG:3857',transform=g['transform'],
                            width=g['width'],height=g['height'],warp_mem_limit=32,tolerance=1e-6)
                with rasterio.open(self.height) as source,rasterio.open(self.mask) as mask, \
                     WarpedVRT(source,**common,dtype='float32',src_nodata=-255,nodata=np.nan,
                               resampling=sampler.method) as expected_heights, \
                     WarpedVRT(mask,**common,dtype='uint8',src_nodata=None,nodata=0,
                               resampling=Resampling.nearest) as expected_mask:
                    for x in range(g['xmin'],g['xmax']+1):
                        for y in range(g['ymin'],g['ymax']+1):
                            actual,valid=sampler.read(x,y)
                            window=Window((x-g['xmin'])*512,(y-g['ymin'])*512,512,512)
                            expected=expected_heights.read(1,window=window)
                            want=(expected_mask.read(1,window=window)==255)&np.isfinite(expected)
                            np.testing.assert_array_equal(valid,want)
                            np.testing.assert_allclose(actual[valid],expected[valid],rtol=0,atol=.0001)
                            count+=int(valid.sum())
                self.assertGreater(g['xmax'],g['xmin'])
                self.assertGreater(count,0)

    def test_packaged_river_metadata_identifies_the_actual_mask_source(self):
        from pmtiles.reader import Reader
        atomic_json(self.root/'rivers/pilot/cogs.json',
                    dict(river_source='glo30',quality_status='Candidate test data'))
        tiles(self.root,'pilot',12,12)
        output=package(self.root,'pilot',max_bytes=1000000)
        data=(self.root/'output/pilot'/output[0]['file']).read_bytes()
        metadata=Reader(lambda offset,length:data[offset:offset+length]).metadata()
        self.assertEqual(metadata['river_mask'],
                         'Copernicus GLO-30 WBM class 3 (River); class 2 (Lake) excluded')
        self.assertFalse(metadata['nodata_preserved'])
        self.assertEqual(metadata['nodata_height_m'],0)
        self.assertNotIn('nodata_alpha',metadata)
        self.assertNotIn('WRIS',metadata['attribution'])
        recipe=json.loads((self.root/'tiles/pilot/recipe.json').read_text())
        self.assertEqual(recipe['river_source'],'glo30')

    def test_parallel_tiles_match_serial_and_resume_an_incomplete_column(self):
        rows,cols=np.indices((96,96))
        values=(rows*.7+cols*1.3).astype('float32');values[30:50,30:50]=-255
        path=self.root/'sheets/coastal/parallel.tif'
        with rasterio.open(path,'w',driver='GTiff',count=1,dtype='float32',
                           width=96,height=96,nodata=-255,crs='EPSG:4326',
                           transform=from_origin(75.045,13,1/3600,1/3600)) as dst:dst.write(values,1)
        create_mosaic([path],self.height,self.mask)
        # One row per task exercises partial-column commits on a tiny 2x2 grid.
        with patch('production_tiles.ROWS_PER_TASK',1):
            tiles(self.root,'pilot',12,12,output_prefix='serial',workers=1)
            result=tiles(self.root,'pilot',12,12,output_prefix='parallel',workers=2)
            self.assertEqual(result[0]['workers'],2)
            def contents(name):
                with sqlite3.connect(self.root/f'tiles/{name}/z12.sqlite') as db:
                    return db.execute('SELECT tileid,data,valid,error FROM tiles ORDER BY tileid').fetchall()
            expected=contents('serial')
            self.assertEqual(contents('parallel'),expected)
            database=self.root/'tiles/parallel/z12.sqlite'
            with sqlite3.connect(database) as db:
                x=db.execute('SELECT min(x) FROM tiles').fetchone()[0]
                # Simulate a partial column whose last batch marker was never written.
                db.execute('DELETE FROM batches WHERE x=?',(x,))
                db.execute('DELETE FROM tiles WHERE x=? AND y=(SELECT min(y) FROM tiles WHERE x=?)',(x,x))
                db.execute('UPDATE tiles SET data=? WHERE x=?',(b'incomplete',x))
            tiles(self.root,'pilot',12,12,output_prefix='parallel',workers=1)
            self.assertEqual(contents('parallel'),expected)
        serial=package(self.root,'serial');parallel=package(self.root,'parallel')
        self.assertEqual([p['sha256'] for p in serial],[p['sha256'] for p in parallel])




if __name__=='__main__':
    unittest.main()
