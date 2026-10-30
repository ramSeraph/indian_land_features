import tempfile
import unittest
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin

from build import atomic_json,sha256,setup
from river_flattening.production import geometry,fit,apply
from river_flattening.audit import audit


class CrossSheetTests(unittest.TestCase):
    def test_shared_profile_crosses_sheet_boundary_and_exports_cogs(self):
        with tempfile.TemporaryDirectory() as directory:
            root=setup(Path(directory)/'work')
            sheets=[]; masks=[]; originals={}
            for tile in range(2):
                name=f'sheet{tile}.tif';mask_name=f'mask{tile}.tif'
                transform=from_origin(75+tile*128/3600,13,1/3600,1/3600)
                rows,cols=np.indices((96,128))
                values=(70-(cols+tile*128)*.16).astype(np.float32)
                classes=np.zeros(values.shape,dtype=np.uint8)
                classes[42:51,:]=3
                classes[12:16,20:26]=2
                if tile==0:
                    values[44:48,110:114]=-255
                originals[name]=values.copy()
                profile=dict(driver='GTiff',width=128,height=96,count=1,crs='EPSG:4326',transform=transform)
                path=root/'sheets/coastal'/name
                with rasterio.open(path,'w',**profile,dtype='float32',nodata=-255) as dst:
                    dst.write(values,1)
                mask=root/'inputs/glo30-wbm'/mask_name
                with rasterio.open(mask,'w',**profile,dtype='uint8',nodata=0) as dst:
                    dst.write(classes,1)
                masks.append(dict(name=mask_name,sha256=sha256(mask)))
                sheets.append(dict(name=name,mask=mask_name,width=128,height=96,transform=list(transform)[:6]))
                atomic_json(root/'reports/sheets'/(Path(name).stem+'.json'),dict(output_sha256=sha256(path)))
            atomic_json(root/'reports/inputs.json',dict(sheets=sheets,masks=masks))
            geometry(root,'check',halo=32,river_source='glo30',allocation_backend='scipy')
            report=fit(root,'check',spacing=120,max_support_gap=600,use_routing=False)
            self.assertGreater(report['fitted_stations'],10)
            result=apply(root,'check')
            self.assertEqual(len(result['sheets']),2)
            self.assertEqual(result['nodata_pixels'],16)
            outputs=[]
            for sheet in sheets:
                with rasterio.open(root/'output/check/cogs'/sheet['name']) as src:
                    array=src.read(1)
                    self.assertEqual(src.nodata,-255)
                    self.assertEqual(src.tags(ns='IMAGE_STRUCTURE')['LAYOUT'],'COG')
                before=originals[sheet['name']]
                preserve=np.ones(before.shape,dtype=bool); preserve[42:51,:]=False
                preserve|=before==-255
                np.testing.assert_array_equal(array[preserve],before[preserve])
                outputs.append(array)
            joined=np.concatenate(outputs,axis=1)
            # This check spans the actual sheet boundary at columns 127/128.
            self.assertTrue(np.all(np.diff(joined[46,120:136])<=.001))
            self.assertLess(abs(float(joined[46,127]-joined[46,128])),.3)
            checked=audit(root,'check')
            self.assertTrue(all(r['maximum_treated_boundary_correction_m']<=.001
                                for r in checked['native_treatment_joins']))


if __name__=='__main__':
    unittest.main()
