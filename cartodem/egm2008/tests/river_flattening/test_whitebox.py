import unittest

import numpy as np

from river_flattening.whitebox import TO_DIRECTION, uncertain_upstream
from river_flattening.routing import OFFSETS


class WhiteboxEvidenceTests(unittest.TestCase):
    def test_compass_convention(self):
        # NE, E, SE, S, SW, W, NW, N in Whitebox order.
        actual = OFFSETS[TO_DIRECTION[[1,2,4,8,16,32,64,128]]]
        np.testing.assert_array_equal(actual, [(-1,1),(0,1),(1,1),(1,0),(1,-1),(0,-1),(-1,-1),(-1,0)])

    def test_unknown_upstream_area_propagates_downstream_from_crop_and_holes(self):
        pointer = np.zeros((9,9),np.uint8)
        valid = np.ones(pointer.shape,bool)
        pointer[4,:5] = 2  # West crop edge flows east to an interior terminal.
        valid[2,6] = False
        pointer[3:6,6] = 8  # Missing interior data has a downstream influence.
        marked = uncertain_upstream(pointer,valid)
        self.assertTrue(marked[4,:6].all())
        self.assertTrue(marked[3:7,6].all())
        self.assertFalse(marked[2,2])
        self.assertFalse(marked[2,6])

    def test_contamination_trace_terminates_on_a_cycle(self):
        pointer = np.zeros((7,7),np.uint8)
        valid = np.ones(pointer.shape,bool)
        valid[2,2] = False
        pointer[3,3] = 2
        pointer[3,4] = 32
        marked = uncertain_upstream(pointer,valid)
        self.assertTrue(marked[3,3:5].all())
        self.assertFalse(marked[3,5])

    def test_invalid_codes_are_rejected(self):
        pointer=np.zeros((3,3),np.uint8);pointer[1,1]=3
        with self.assertRaises(ValueError):
            uncertain_upstream(pointer,np.ones(pointer.shape,bool))


if __name__=='__main__':unittest.main()
