import unittest
import numpy as np

from river_flattening.joint import fit_joint,collapse_levels
from river_flattening.test_directions import network


def observations(levels):
    return dict(bank_owners=np.repeat(np.arange(len(levels)),5),
                bank_heights=np.repeat(np.asarray(levels,float),5),ocean_owners=np.array([],int))


class JointTests(unittest.TestCase):
    def test_flat_middle_branches_are_in_the_same_fit_as_adjoining_reaches(self):
        net=network([[0,4,1],[1,5,2],[1,6,7,2],[2,8,3]])
        native=observations([20,10,10,0,15,10,10,10,5])
        result,reports=fit_joint(native,net)
        self.assertEqual(len(reports),1)
        self.assertEqual(reports[0]['status'],'fitted')
        self.assertTrue(np.isfinite(result['fitted']).all())
        self.assertTrue(result['native_supported'].all())
        self.assertTrue(all(r['profile_applied'] for r in result['direction_diagnostics']))
        self.assertGreater(reports[0]['joint_direction_reaches'],0)

    def test_raw_water_heights_are_not_exact_boundary_targets(self):
        net=network([[0,1,2,3]])
        table=dict(owners=np.array([[1,2]]),heights=np.array([[500,-20]],dtype='float32'))
        result,reports=fit_joint(observations([10,7,4,1]),net,contacts=[table])
        self.assertEqual(reports[0]['status'],'fitted')
        self.assertFalse(result['boundary_acceptance']['raw_height_pins'])
        np.testing.assert_allclose(result['fitted'],[10,7,4,1],atol=.001)

    def test_two_outlets_still_share_the_same_upstream_unknown(self):
        net=network([[0,4,1],[1,5,2],[1,6,3]])
        native=observations([20,10,2,2,15,6,6]);native['ocean_owners']=np.array([2,3])
        result,reports=fit_joint(native,net)
        self.assertEqual(reports[0]['status'],'fitted')
        edge=result['supported_edges']
        self.assertTrue((result['fitted'][edge[:,0]]>=result['fitted'][edge[:,1]]).all())
        self.assertEqual(len(reports),1)

    def test_unmodelled_contact_retains_connected_group_without_raw_height_pins(self):
        net=network([[0,1,2]])
        net['ids']=np.arange(4);net['left']=np.array([0,1,2,-1]);net['right']=net['left'].copy();net['fraction']=np.zeros(4)
        table=dict(owners=np.array([[2,3]]),heights=np.array([[1,0]],dtype='float32'))
        result,reports=fit_joint(observations([7,4,1]),net,contacts=[table])
        self.assertTrue(np.isnan(result['fitted']).all())
        self.assertFalse(result['native_supported'].any())
        self.assertEqual(result['boundary_acceptance']['rejected_profile_groups'],[0])

    def test_exact_level_links_collapse_without_arbitrary_flow_direction(self):
        mapping,edges=collapse_levels(4,[(0,1),(2,3)],[[1,2]])
        self.assertEqual(mapping[1],mapping[2])
        self.assertEqual(len(edges),2)
        mapping,edges=collapse_levels(3,[(0,1),(1,2)],[[0,2]])
        self.assertTrue((mapping==mapping[0]).all())
        self.assertEqual(len(edges),0)


if __name__=='__main__':unittest.main()
