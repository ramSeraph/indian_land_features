"""Checks for rebuilding display tiles without changing COGs or invoking legacy packaging."""
import io
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import rasterio
from rasterio.shutil import copy as rio_copy
from rasterio.transform import from_origin
from PIL import Image

from build import atomic_json,read_json,sha256
from production_tiles import Sampler,validate_rgba
from rebuild_tiles import rebuild


class RebuildTests(unittest.TestCase):
    def test_existing_cogs_are_unchanged_and_display_tiles_resume(self):
        with tempfile.TemporaryDirectory() as temp:
            parent=Path(temp);cogs=parent/'old/cogs';cogs.mkdir(parents=True)
            entries=[]
            for i in range(2):
                raw=parent/f'input{i}.tif';path=cogs/f'cdn{i}.tif'
                data=np.full((16,16),i*16,dtype='float32');data[5:10,5:10]=-255
                with rasterio.open(raw,'w',driver='GTiff',count=1,dtype='float32',
                                   width=16,height=16,nodata=-255,crs='EPSG:4326',
                                   transform=from_origin(75+i*16/3600,13,1/3600,1/3600)) as dst:dst.write(data,1)
                rio_copy(raw,path,driver='COG',compress='DEFLATE',blocksize=512)
                entries.append(dict(file=path.name,sha256=sha256(path),bytes=path.stat().st_size))
            manifest=parent/'old/cogs.json'
            atomic_json(manifest,dict(sheets=entries,river_source='glo30',quality_status='Test'))
            root=parent/'new'
            result=rebuild(root,cogs,manifest,expected_sheets=2,min_zoom=11,max_zoom=12)
            self.assertEqual(result['sheets'],2)
            database=root/'tiles/display/z12.sqlite';before=sha256(database)
            with Sampler(root/'mosaics/display.vrt',root/'mosaics/display-validity.vrt',12) as sampler,sqlite3.connect(database) as db:
                for x,y,data,count in db.execute('SELECT x,y,data,valid FROM tiles'):
                    _,valid=sampler.read(x,y)
                    rgba=np.asarray(Image.open(io.BytesIO(data)).convert('RGBA'))
                    self.assertEqual(validate_rgba(rgba,12,valid),count)
            again=rebuild(root,cogs,manifest,expected_sheets=2,min_zoom=11,max_zoom=12)
            self.assertEqual(again,result);self.assertEqual(sha256(database),before)
            for entry in entries:
                self.assertEqual(sha256(cogs/entry['file']),entry['sha256'])
                with rasterio.open(cogs/entry['file']) as src:self.assertEqual(src.nodata,-255)
            meta=read_json(root/'output/display/archives.json')['metadata']
            self.assertEqual(meta['alpha'],255);self.assertFalse(meta['nodata_preserved'])
            with patch('rebuild_tiles.package',side_effect=AssertionError('Legacy packager must not run')):
                tile_result=rebuild(root,cogs,manifest,expected_sheets=2,min_zoom=11,max_zoom=12,tiles_only=True)
            self.assertEqual(tile_result['tiles'],result['tiles'])
            self.assertEqual(tile_result['archives'],0)
            self.assertEqual(sha256(database),before)
            self.assertEqual(read_json(root/'output/display/archives.json')['partitions'],[])
            with self.assertRaisesRegex(ValueError,'Expected 543'):
                rebuild(parent/'wrong',cogs,manifest)
            with (cogs/entries[0]['file']).open('ab') as f:f.write(b'changed')
            with self.assertRaisesRegex(ValueError,'differs'):
                rebuild(root,cogs,manifest,expected_sheets=2,min_zoom=11,max_zoom=12)



if __name__=='__main__':unittest.main()
