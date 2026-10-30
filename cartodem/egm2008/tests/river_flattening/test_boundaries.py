import unittest
import numpy as np
from river_flattening.boundaries import adjacent_contacts
from validation.river_audit import audit_native_joins


class BoundaryTests(unittest.TestCase):

    def test_shared_sheet_border_checked_without_wrapping_raster_rows(self):
        data=adjacent_contacts(np.array([8,9,10,18,19]),np.arange(5),np.ones(5),10,
                               sheets=np.array([0,0,1,1,1]))
        self.assertEqual({tuple(p) for p in data['owners']},{(0,3),(1,4)})


    def test_raster_audit_detects_a_bad_written_join_across_sheet_border(self):
        args=(np.array([8,9,10]),np.array([5,4,3],dtype='float32'),
              np.array([7,4,3],dtype='float32'),np.array([1,0,0],bool),np.array([0,-1,-1]),10)
        result=audit_native_joins(*args,sheets=np.array([0,1,1]))
        self.assertEqual(result['boundary_pairs'],1)
        self.assertEqual(result['maximum_treated_boundary_correction_m'],2)
        self.assertEqual(result['maximum_signed_step_change_m'],2)


if __name__=='__main__': unittest.main()
