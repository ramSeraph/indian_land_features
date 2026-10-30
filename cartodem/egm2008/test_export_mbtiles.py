import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest

from export_mbtiles import export_zoom


class ExportTests(unittest.TestCase):
    def test_xyz_rows_become_tms_and_image_bytes_are_unchanged(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); source=root/'z3.sqlite'; target=root/'z3.mbtiles'
            with sqlite3.connect(source) as db:
                db.execute('CREATE TABLE tiles(tileid INTEGER PRIMARY KEY,x INTEGER,y INTEGER,data BLOB,sha256 TEXT)')
                for i,(x,y,data) in enumerate([(2,1,b'image A'),(2,6,b'image B')]):
                    db.execute('INSERT INTO tiles VALUES(?,?,?,?,?)',(i,x,y,data,hashlib.sha256(data).hexdigest()))
            first=export_zoom(source,target,3,dict(format='webp',name='test'),2)
            with sqlite3.connect(target) as db:
                self.assertEqual(db.execute('SELECT * FROM tiles ORDER BY tile_row').fetchall(),
                                 [(3,2,1,b'image B'),(3,2,6,b'image A')])
                self.assertEqual(dict(db.execute('SELECT * FROM metadata'))['minzoom'],'3')
            self.assertEqual(export_zoom(source,target,3,dict(format='webp',name='test'),2),first)
            with sqlite3.connect(target) as db:
                db.execute("UPDATE tiles SET tile_data=? WHERE tile_row=1",(b'bad',))
            with self.assertRaisesRegex(ValueError,'coordinates or image bytes'):
                export_zoom(source,target,3,dict(format='webp',name='test'),2)


if __name__=='__main__': unittest.main()
