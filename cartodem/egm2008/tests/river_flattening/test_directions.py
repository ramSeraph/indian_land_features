import unittest
from unittest.mock import patch

import numpy as np

from river_flattening.directions import infer_directions
from river_flattening.routing import analysis_window
from river_flattening.vendor.rivgraph import dir_shortest_paths_links


def network(reaches):
    n=max(max(r) for r in reaches)+1
    return dict(reaches=reaches,distances=[np.arange(len(r))*180. for r in reaches],
                station_pixels=np.arange(n),ids=np.arange(n),left=np.arange(n),right=np.arange(n),
                fraction=np.zeros(n),unsupported_skeleton_pixels=0)


class DirectionTests(unittest.TestCase):
    def test_two_sea_outlets_share_one_upstream_network(self):
        net=network([[0,4,1],[1,5,2],[1,6,3]])
        reference=np.array([20,10,2,2,15,6,6.])
        direction,report=infer_directions(net,reference,np.array([2,3]))
        np.testing.assert_array_equal(direction,[1,1,1])
        self.assertTrue(all(r['status']=='accepted' for r in report))

    def test_unequal_braid_preserves_both_downhill_branches(self):
        net=network([[0,4,1],[1,5,2],[1,6,7,8,2],[2,9,3]])
        ref=np.array([30,20,10,0,25,15,18,15,12,5.])
        direction,_=infer_directions(net,ref,np.array([3]))
        np.testing.assert_array_equal(direction,[1,1,1,1])

    def test_geometry_alone_does_not_resolve_flat_river(self):
        net=network([[0,1,2,3]])
        a,report=infer_directions(net,np.ones(4)*7,np.array([],int))
        np.testing.assert_array_equal(a,[0])
        self.assertEqual(report[0]['status'],'unresolved')

    def test_bank_routing_conflict_remains_unknown(self):
        net=network([[0,1,2,3]])
        with patch('river_flattening.directions.reach_routing',return_value=[(-1,dict(agreement=1.))]):
            direction,report=infer_directions(net,np.array([15,10,5,0.]),np.array([],int),routing={},width=10)
        self.assertEqual(direction[0],0)
        self.assertIn('conflicts',report[0]['reason'])


    def test_equal_outlet_distances_do_not_break_tie_by_link_order(self):
        links=dict(id=[0,1,2],conn=[[0,1],[1,2],[0,3]],len=[10,10,10],guess=[[],[],[]],guess_alg=[[],[],[]])
        nodes=dict(id=[0,1,2,3],conn=[[0,2],[0,1],[1],[2]],outlets=[2,3])
        dir_shortest_paths_links(links,nodes)
        self.assertEqual(links['guess'][0],[])

    def test_analysis_fill_is_separate_from_original_holes_and_sea(self):
        height=np.full((80,80),12.,dtype=np.float32);classes=np.zeros((80,80),np.uint8)
        height[30:33,30:33]=-255;classes[:,:5]=1;height[:,:5]=-255
        original=height.copy()
        analysis,evidence,filled=analysis_window(height,classes)
        np.testing.assert_array_equal(height,original)
        self.assertTrue(np.all(analysis[30:33,30:33]==12))
        self.assertFalse(evidence[30:33,30:33].any())
        self.assertTrue(np.all(analysis[:,:5]==-255))
        self.assertFalse(filled[:,:5].any())


if __name__=='__main__':
    unittest.main()
