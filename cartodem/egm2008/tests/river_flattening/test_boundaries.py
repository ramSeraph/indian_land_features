import unittest
import numpy as np
from river_flattening.boundaries import adjacent_contacts,guard,constrain_joins
from validation.river_audit import audit_native_joins


class BoundaryTests(unittest.TestCase):
    def test_reverts_whole_shared_profile_at_bad_native_transition(self):
        net=dict(ids=np.arange(5),left=np.arange(5),right=np.arange(5),fraction=np.zeros(5),
                 fitted=np.array([0,0,np.nan,5,5],dtype='float32'),native_supported=np.array([1,1,0,1,1],bool),
                 component=np.array([0,0,-1,1,1]),supported_edges=np.array([[0,1],[3,4]]))
        table=dict(owners=np.array([[1,2],[3,4]]),heights=np.array([[14,14],[5,5]],dtype='float32'))
        result=guard(net,[table])
        self.assertEqual(result['rejected_profile_groups'],[0])
        self.assertTrue(np.isnan(net['fitted'][:3]).all())
        np.testing.assert_array_equal(net['fitted'][3:],[5,5])
        np.testing.assert_array_equal(net['supported_edges'],[[3,4]])
        self.assertEqual(result['remaining_bad_pairs'],0)

    def test_shared_sheet_border_checked_without_wrapping_raster_rows(self):
        data=adjacent_contacts(np.array([8,9,10,18,19]),np.arange(5),np.ones(5),10,
                               sheets=np.array([0,0,1,1,1]))
        self.assertEqual({tuple(p) for p in data['owners']},{(0,3),(1,4)})

    def test_constrained_boundary_preserves_original_slope_without_copying_neighbour(self):
        net=dict(ids=np.arange(3),left=np.arange(3),right=np.arange(3),fraction=np.zeros(3),
                 fitted=np.array([3,2,np.nan],dtype='float32'),native_supported=np.array([1,1,0],bool),
                 component=np.array([0,0,-1]),supported_edges=np.array([[0,1]]))
        table=dict(owners=np.array([[1,2]]),heights=np.array([[5,4]],dtype='float32'))
        calls=[]
        def refit(group,samples):
            calls.append((group,samples.copy()))
            net['fitted'][:2]=[6,5]
        result=constrain_joins(net,[table],refit)
        self.assertEqual(result['rejected_profile_groups'],[])
        np.testing.assert_array_equal(calls[0][1],[[1,1,0,5]])
        self.assertEqual(float(net['fitted'][1])-4,1)  # Original 1 m gradient remains.
        self.assertEqual(result['maximum_boundary_mismatch_m'],0)

    def test_incompatible_constraint_retains_whole_group_and_rechecks_contacts(self):
        net=dict(ids=np.arange(3),left=np.arange(3),right=np.arange(3),fraction=np.zeros(3),
                 fitted=np.array([3,2,np.nan],dtype='float32'),native_supported=np.array([1,1,0],bool),
                 component=np.array([0,0,-1]),supported_edges=np.array([[0,1]]))
        table=dict(owners=np.array([[1,2]]),heights=np.array([[5,4]],dtype='float32'))
        def refit(group,samples):raise ValueError('Infeasible bank cap')
        result=constrain_joins(net,[table],refit)
        self.assertEqual(result['rejected_profile_groups'],[0])
        self.assertTrue(np.isnan(net['fitted']).all())
        self.assertFalse(net['native_supported'].any())
        self.assertEqual(len(result['iterations']),2)

    def test_raster_audit_detects_a_bad_written_join_across_sheet_border(self):
        args=(np.array([8,9,10]),np.array([5,4,3],dtype='float32'),
              np.array([7,4,3],dtype='float32'),np.array([1,0,0],bool),np.array([0,-1,-1]),10)
        result=audit_native_joins(*args,sheets=np.array([0,1,1]))
        self.assertEqual(result['boundary_pairs'],1)
        self.assertEqual(result['maximum_treated_boundary_correction_m'],2)
        self.assertEqual(result['maximum_signed_step_change_m'],2)


if __name__=='__main__': unittest.main()
