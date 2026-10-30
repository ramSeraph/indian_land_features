import tempfile
import unittest
from pathlib import Path

import numpy as np
import rasterio

from river_flattening.allocation import read_owners, write_inputs


class AllocationTests(unittest.TestCase):
    def test_block_label_mapping_preserves_holes_and_seed_coordinates(self):
        # Both source labels and ownership cross TIFF block boundaries.
        river = np.zeros((270,280), bool)
        river[2:268,4:277] = True
        river[100:110,120:130] = False
        seeds = np.zeros_like(river)
        seeds[5,5] = seeds[260,270] = True
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            info = write_inputs(folder,river,seeds,23.,31.)
            self.assertEqual(info['river_pixels'],int(river.sum()))
            with rasterio.open(folder/'sources.tif') as src:
                labels = src.read(1);profile = src.profile.copy()
                self.assertEqual(labels[5,5],1);self.assertEqual(labels[260,270],2)
                self.assertTrue(np.all(labels[~river]==-1))
                self.assertEqual(src.res,(23.,31.))
            labels[river]=np.where(np.indices(river.shape)[0][river]<135,1,2)

            def output():
                with rasterio.open(folder/'owners.tif','w',**profile) as dst:
                    dst.write(labels,1)

            output()
            expected=np.flatnonzero(seeds)[labels[river]-1]
            np.testing.assert_array_equal(read_owners(folder,river,seeds),expected)
            labels[105,125]=1
            output()
            with self.assertRaisesRegex(ValueError,'crossed the river mask'):
                read_owners(folder,river,seeds)
            labels[105,125]=-1;labels[5,5]=2
            output()
            with self.assertRaisesRegex(ValueError,'seed was assigned'):
                read_owners(folder,river,seeds)
            labels[5,5]=1;labels[50,50]=0
            output()
            with self.assertRaisesRegex(ValueError,'no valid centreline owner'):
                read_owners(folder,river,seeds)

    def test_rejects_seeds_outside_river_and_invalid_metric(self):
        with tempfile.TemporaryDirectory() as directory:
            river=np.eye(8,dtype=bool);seeds=river.copy();seeds[0,1]=True
            with self.assertRaisesRegex(ValueError,'Invalid centreline seeds'):
                write_inputs(directory,river,seeds,30.,30.)
            with self.assertRaisesRegex(ValueError,'Invalid allocation inputs'):
                write_inputs(directory,river,river.copy(),np.nan,30.)


if __name__=='__main__':unittest.main()
