import tempfile
import unittest
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin

from cog import export_cog
from build import create_mosaic


class COGTests(unittest.TestCase):
    def test_cog_keeps_native_values_and_pixel_centre_holes(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            values=np.full((1024,1024),12.5,dtype=np.float32)
            values[1::2,1::2]=-255
            values[50,50]=0
            values[60,60]=-7.25
            source=root/'source.tif';target=root/'final.tif'
            with rasterio.open(source,'w',driver='GTiff',width=1024,height=1024,count=1,
                               dtype='float32',crs='EPSG:4326',transform=from_origin(75,13,1/3600,1/3600),
                               nodata=-255) as dst:
                dst.write(values,1)
            report=export_cog(source,target)
            self.assertTrue(report['native_pixels_verified'])
            self.assertEqual(report['overviews'],[2])
            create_mosaic([target],root/'heights.vrt',root/'validity.vrt')
            with rasterio.open(root/'validity.vrt') as valid:
                mask=valid.read(1,out_shape=(512,512))
            # Nearest-native samples land on odd pixel centres, all nodata.
            # The averaged numeric COG overview contains valid elevations there.
            np.testing.assert_array_equal(mask,np.zeros((512,512),dtype=np.uint8))


if __name__=='__main__':
    unittest.main()
