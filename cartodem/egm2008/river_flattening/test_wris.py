import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import Polygon, box

from river_flattening.wris import build_index, river_mask, sheet_mask


class WrisTests(unittest.TestCase):
    def test_river_holes_and_sand_islands_override_overlapping_water(self):
        t = from_origin(0, 10, 1, 1)
        river = Polygon([(0,0),(10,0),(10,10),(0,10)],
                        holes=[[(2,2),(4,2),(4,4),(2,4)]])
        features = [(1,'01',river),(2,'03',box(2,2,4,4)),(3,'02',box(6,6,8,8)),
                    (4,None,box(2,2,4,4)),(5,' ',box(0,0,1,1))]
        mask, report = river_mask(features, (10,10), t)
        self.assertEqual(int(mask.sum()),92)
        self.assertFalse(mask[6:8,2:4].any())
        self.assertFalse(mask[2:4,6:8].any())
        self.assertEqual(report['excluded_overlap_pixels'],8)

    def test_pixel_centres_stitch_across_sheet_boundary(self):
        features = [(1,'01',box(1.2,1.2,18.2,8.2)), (2,'03',box(8.2,3.2,12.2,5.2))]
        whole, _ = river_mask(features, (10,20), from_origin(0,10,1,1))
        left, _ = river_mask(features, (10,10), from_origin(0,10,1,1))
        right, _ = river_mask(features, (10,10), from_origin(10,10,1,1))
        np.testing.assert_array_equal(whole,np.concatenate((left,right),axis=1))
        self.assertFalse(whole[0].any())

    def test_unknown_codes_and_invalid_geometry_fail(self):
        with self.assertRaisesRegex(ValueError,'ripcode'):
            river_mask([(1,'99',box(0,0,1,1))],(10,10),from_origin(0,10,1,1))
        bow = Polygon([(1,1),(9,9),(9,1),(1,9),(1,1)])
        with self.assertRaisesRegex(ValueError,'Invalid WRIS'):
            river_mask([(1,'01',bow)],(10,10),from_origin(0,10,1,1))
        with self.assertRaisesRegex(ValueError,'more than one native pixel'):
            river_mask([(1,'01',bow)],(10,10),from_origin(0,10,1,1),normalize_invalid=True)

    def test_self_touch_normalization_preserves_native_membership(self):
        # Two squares sharing a vertex, encoded as an invalid single ring.
        polygon=Polygon([(1,1),(3,1),(3,3),(5,3),(5,5),(3,5),(3,3),(1,3),(1,1)])
        self.assertFalse(polygon.is_valid)
        transform=from_origin(0,10,1,1)
        fixed,report=river_mask([(7,'01',polygon)],(10,10),transform,normalize_invalid=True)
        expected,_=river_mask([(1,'01',box(1,1,3,3)),(2,'01',box(3,3,5,5))],(10,10),transform)
        np.testing.assert_array_equal(fixed,expected)
        self.assertEqual(report['normalized_invalid_polygons'][0]['objectid'],7)
        self.assertEqual(report['normalized_invalid_polygons'][0]['native_membership_changes'],0)

    def test_parquet_index_and_native_classification_precedence(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=root/'river.parquet';index=root/'river.sqlite'
            table=pa.table({'objectid':[1,2,3],'ripcode':['01','03',None],
                            'geometry':[box(0,0,10,10).wkb,box(3,3,5,5).wkb,box(20,0,21,1).wkb]})
            metadata={'version':'1.1.0','primary_column':'geometry',
                      'columns':{'geometry':{'encoding':'WKB','geometry_types':['Polygon']}}}
            table=table.replace_schema_metadata({b'geo':json.dumps(metadata).encode()})
            pq.write_table(table,source)
            report=build_index(source,index)
            self.assertEqual(report['features'],3)
            self.assertEqual(build_index(source,index),report)
            classes=np.zeros((10,10),dtype=np.uint8);classes[0]=1;classes[1]=2
            classes[8]=3
            glo=root/'glo.tif';out=root/'classes.tif'
            with rasterio.open(glo,'w',driver='GTiff',height=10,width=10,count=1,
                               dtype='uint8',nodata=0,crs='EPSG:4326',transform=from_origin(0,10,1,1)) as dst:
                dst.write(classes,1)
            sheet_mask(index,glo,out)
            with rasterio.open(out) as src:
                result=src.read(1)
                self.assertEqual(src.nodata,255)
            np.testing.assert_array_equal(result[0],classes[0])
            # A GLO lake classification cannot override WRIS's chosen river.
            self.assertTrue((result[1]==3).all())
            self.assertTrue((result[5:7,3:5]==0).all())
            self.assertEqual(int((result==3).sum()),86)


if __name__=='__main__': unittest.main()
