"""Checks for rebuilding without changing COGs and resuming verified uploads."""
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
from upload_release import assets,upload


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
            self.assertEqual(len(assets(root,'display','cogs')),4)
            self.assertEqual(len(assets(root,'display','pmtiles')),result['archives']+5)
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

    def test_upload_skips_matching_assets_and_checks_new_digests(self):
        old=dict(name='old.tif',path=Path('old.tif'),bytes=5,sha256='a'*64)
        new=dict(name='new.tif',path=Path('new.tif'),bytes=6,sha256='b'*64)
        old_remote=dict(name=old['name'],size=5,digest='sha256:'+old['sha256'])
        new_remote=dict(name=new['name'],size=6,digest='sha256:'+new['sha256'])
        with patch('upload_release.gh_json',side_effect=[{'id':7},[[old_remote]],[[old_remote,new_remote]]]),patch('upload_release.subprocess.run') as run:
            upload('gh','owner/repo','tag',[old,new])
            run.assert_called_once_with(['gh','release','upload','tag','new.tif','--repo','owner/repo'],check=True)
        with patch('upload_release.gh_json',side_effect=[{'id':7},[[dict(old_remote,digest='sha256:wrong')]]]),patch('upload_release.subprocess.run') as run:
            with self.assertRaisesRegex(ValueError,'differs'):
                upload('gh','owner/repo','tag',[new,old])
            run.assert_not_called()


if __name__=='__main__':unittest.main()
