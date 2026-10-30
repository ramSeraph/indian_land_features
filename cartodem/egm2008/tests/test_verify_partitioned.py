"""Verify archived bytes and keep the application-check index private."""
import hashlib
import io
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest.mock import patch

from PIL import Image
from pmtiles.tile import Compression, TileType, zxy_to_tileid
from pmtiles.writer import Writer

from validation.verify_partitioned import main


class PartitionVerificationTests(unittest.TestCase):
    def test_exact_tiles_create_private_index_and_changed_bytes_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for sub in ('tiles/display', 'reports', 'output/display', 'pmtiles'):
                (root/sub).mkdir(parents=True)
            image = io.BytesIO()
            Image.new('RGB', (2, 2), (1, 134, 160)).save(image, format='WEBP', lossless=True)
            data = image.getvalue()
            digest = hashlib.sha256(data).hexdigest()
            tile = zxy_to_tileid(3, 4, 3)
            database = root/'tiles/display/z3.sqlite'
            with sqlite3.connect(database) as db:
                db.execute('CREATE TABLE tiles(tileid INTEGER PRIMARY KEY,sha256 TEXT)')
                db.execute('INSERT INTO tiles VALUES (?,?)', (tile, digest))
            (root/'reports/display-tiles.json').write_text(json.dumps({'zooms': [{'zoom': 3}]}))
            (root/'output/display/archives.json').write_text(json.dumps({'metadata': {'format': 'webp'}}))
            archive = root/'pmtiles/example.pmtiles'
            with archive.open('wb') as stream:
                writer = Writer(stream)
                writer.write_tile(tile, data)
                writer.finalize(dict(tile_type=TileType.WEBP, tile_compression=Compression.NONE,
                                     min_zoom=3, max_zoom=3, min_lon_e7=0, min_lat_e7=0,
                                     max_lon_e7=450000000, max_lat_e7=409798981), {})
            args = ['validation.verify_partitioned', '--workdir', str(root), '--pmtiles-dir', str(root/'pmtiles')]
            with patch('sys.argv', args):
                main()
            info = json.loads((root/'pmtiles/.validation/partition-info.json').read_text())
            self.assertEqual(info['example.pmtiles']['tiles'], {'3,4,3': len(data)})
            self.assertEqual(list((root/'pmtiles').glob('*.partition_info.json')), [])
            self.assertEqual(json.loads((root/'pmtiles/archives.json').read_text())['verified_tiles'], 1)
            with sqlite3.connect(database) as db:
                db.execute('UPDATE tiles SET sha256=?', ('0'*64,))
            with patch('sys.argv', args), self.assertRaisesRegex(ValueError, 'Changed, duplicate or unexpected'):
                main()


if __name__ == '__main__':
    unittest.main()
