import json
import sqlite3
import tempfile
import unittest
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin

from build import atomic_json,read_json,sha256,setup
from river_flattening.prepare import geometry
from river_flattening.national import fit,apply,skeleton_components
from river_flattening.regional import accumulate,route_regions
from river_flattening.staging import geometry_sheet
from .helpers import attach_whitebox


class NationalTests(unittest.TestCase):
    def test_routing_excludes_incomplete_catchments_and_invalid_samples(self):
        result=accumulate(np.array([10,10,20,20,30]),np.array([8,8,4,4,2]),
                          np.array([10.,-10.,20.,20.,30.]),np.array([1,1,1,0,0],bool))
        np.testing.assert_array_equal(result['nodes'],[10,20])
        np.testing.assert_array_equal(result['count'],[1,1])
        np.testing.assert_allclose(result['vector'][:,0],0)
        self.assertGreater(result['vector'][0,1],0)
        self.assertLess(result['vector'][1,1],0)

    def test_component_graph_does_not_wrap_at_row_edges(self):
        count,labels=skeleton_components(np.array([9,10,50,51]),10)
        self.assertEqual(count,3)
        self.assertNotEqual(labels[0],labels[1])
        self.assertEqual(labels[2],labels[3])

    def test_whitebox_and_joint_fit_cross_sheets_and_preserve_holes(self):
        with tempfile.TemporaryDirectory() as directory:
            root=setup(Path(directory)/'work');sheets=[];masks=[]
            attach_whitebox(self,root)
            for tile in range(2):
                name=f'sheet{tile}.tif';mask_name=f'mask{tile}.tif'
                transform=from_origin(75+tile*128/3600,13,1/3600,1/3600)
                rows,cols=np.indices((96,128));values=(70-(cols+tile*128)*.16).astype('float32')
                classes=np.zeros(values.shape,'uint8');classes[42:51,:]=3
                classes[:2,:]=1;values[:2,:]=-255
                classes[10:14,20:26]=2  # Lake elevations remain outside the river treatment.
                classes[15:17,20:22]=3  # Separate tiny component remains unchanged.
                # A separate four-pixel closed loop has no station or outlet.
                classes[[24,25,25,26],[21,20,22,21]]=3
                if tile==0:values[44:48,110:114]=-255
                profile=dict(driver='GTiff',width=128,height=96,count=1,crs='EPSG:4326',transform=transform)
                path=root/'sheets/coastal'/name
                with rasterio.open(path,'w',**profile,dtype='float32',nodata=-255) as dst:dst.write(values,1)
                mask=root/'inputs/glo30-wbm'/mask_name
                with rasterio.open(mask,'w',**profile,dtype='uint8',nodata=0) as dst:dst.write(classes,1)
                masks.append(dict(name=mask_name,sha256=sha256(mask)))
                sheets.append(dict(name=name,mask=mask_name,width=128,height=96,transform=list(transform)[:6]))
                atomic_json(root/'reports/sheets'/(Path(name).stem+'.json'),dict(output_sha256=sha256(path)))
            atomic_json(root/'reports/inputs.json',dict(sheets=sheets,masks=masks))
            geometry(root,'partitioned',halo=32)
            self.assertEqual(read_json(root/'rivers/partitioned/geometry-recipe.json')['river_source'],'glo30')
            # The isolated worker must reproduce the same global IDs and bank
            # samples as the inline geometry path, including halo ownership.
            for sheet in sheets:
                path=root/'rivers/partitioned/geometry'/(Path(sheet['name']).stem+'.npz')
                previous=sha256(path)
                geometry_sheet(root,root/'rivers/partitioned',sheet)
                self.assertEqual(sha256(path),previous)
            route_regions(root,'partitioned',halo=64)
            report=fit(root,'partitioned',spacing=120,max_support_gap=600)
            self.assertTrue(report['partitioned']);self.assertGreater(report['fitted_stations'],10)
            with sqlite3.connect(root/'rivers/partitioned/components/models.sqlite') as db:
                loops=[(json.loads(text),model) for text,model in db.execute('SELECT report,data FROM models')
                       if json.loads(text)['stations']==0]
            self.assertEqual(len(loops),2)
            for loop,model in loops:
                self.assertEqual(loop['native_nodes'],4)
                self.assertEqual(loop['components'][0]['reason'],'No modelled river reach')
                self.assertIsNone(model)
            self.assertEqual(fit(root,'partitioned',spacing=120,max_support_gap=600),report)
            result=apply(root,'partitioned')
            self.assertEqual(result['river_source'],'glo30')
            self.assertEqual(result['treated_river_pixels'],2288)
            self.assertEqual(result['nodata_pixels'],528)
            outputs=[]
            for sheet in sheets:
                with rasterio.open(root/'output/partitioned/cogs'/sheet['name']) as src:
                    actual=src.read(1);self.assertEqual(src.nodata,-255)
                    self.assertEqual(src.tags(ns='IMAGE_STRUCTURE')['LAYOUT'],'COG')
                    self.assertEqual(src.tags()['RIVER_MASK'],'glo30')
                with rasterio.open(root/'sheets/coastal'/sheet['name']) as src:
                    before=src.read(1)
                with rasterio.open(root/'inputs/glo30-wbm'/sheet['mask']) as src:classes=src.read(1)
                preserve=(classes!=3)|(before==-255)
                np.testing.assert_array_equal(actual[preserve],before[preserve])
                np.testing.assert_array_equal(actual[24:27,20:23],before[24:27,20:23])
                self.assertLess(float(np.ptp(actual[42:51,70])),.01)
                outputs.append(actual)
            edge_drop=outputs[0][42:51,-1]-outputs[1][42:51,0]
            self.assertTrue(np.all((edge_drop>0)&(edge_drop<.5)))
            self.assertEqual(apply(root,'partitioned')['treated_river_pixels'],result['treated_river_pixels'])


if __name__=='__main__':unittest.main()
