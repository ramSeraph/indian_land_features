import unittest

import numpy as np

from river_flattening.geometry import extract_window,compress_network
from .helpers import allocator


class GeometryTests(unittest.TestCase):
    def test_closed_loop_without_stations_remains_unsupported(self):
        # The four-node loop that stopped the national fit in cdnh43m.tif.
        ids=np.array([4490232075,4490368874,4490368876,4490505675],dtype=np.int64)
        network=compress_network(ids,136800,30.,30.)
        np.testing.assert_array_equal(network['ids'],ids)
        self.assertEqual(network['station_pixels'].size,0)
        self.assertEqual(network['station_pixels'].dtype,np.dtype('int64'))
        self.assertEqual(network['reaches'],[])
        self.assertEqual(network['unsupported_skeleton_pixels'],4)
        np.testing.assert_array_equal(network['left'],[-1]*4)
        np.testing.assert_array_equal(network['right'],[-1]*4)

    def test_halo_windows_share_skeleton_coordinates_and_owners(self):
        classes=np.zeros((50,80),dtype=np.uint8)
        classes[21:28,3:77]=3
        heights=np.full(classes.shape,10.,dtype=np.float32)
        allocate=allocator(self)
        whole=extract_window(heights,classes,dx=30,dy=30,grid_width=80,allocator=allocate)
        west=extract_window(heights[:,:56],classes[:,:56],dx=30,dy=30,grid_width=80,core=(0,50,0,40),allocator=allocate)
        east=extract_window(heights[:,24:],classes[:,24:],dx=30,dy=30,grid_width=80,col_offset=24,core=(0,50,16,56),allocator=allocate)
        joined=np.sort(np.r_[west['nodes'],east['nodes']])
        np.testing.assert_array_equal(joined,whole['nodes'])
        pixels=np.r_[west['pixels'],east['pixels']]
        owners=np.r_[west['owners'],east['owners']]
        order=np.argsort(pixels)
        np.testing.assert_array_equal(pixels[order],whole['pixels'])
        np.testing.assert_array_equal(owners[order],whole['owners'])

    def test_ownership_does_not_jump_across_land_between_parallel_channels(self):
        classes=np.zeros((50,70),dtype=np.uint8)
        classes[15:21,3:67]=3
        classes[26:32,3:67]=3
        heights=np.full(classes.shape,20.,dtype=np.float32)
        result=extract_window(heights,classes,dx=30,dy=30,grid_width=70,allocator=allocator(self))
        north=result['pixels']//70<23
        self.assertTrue(np.all(result['owners'][north]//70<23))
        self.assertTrue(np.all(result['owners'][~north]//70>23))


if __name__=='__main__':
    unittest.main()
