"""A changed river mask must never reuse the previous mask's cached geometry."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class ResumeTests(unittest.TestCase):
    def test_restart_rejects_incompatible_geometry_before_detaching(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            geometry=root/'rivers/final/geometry-recipe.json'
            geometry.parent.mkdir(parents=True)
            geometry.write_text(json.dumps(dict(river_source='unsupported')))
            result=subprocess.run([sys.executable,str(Path(__file__).resolve().parents[1]/'run_build.py'),
                                   '--workdir',str(root),'--resume-from','routing','--detach'],
                                  text=True,capture_output=True)
            self.assertEqual(result.returncode,2)
            self.assertIn('does not match cached geometry',result.stderr)
            self.assertFalse((root/'logs/build-final.log').exists())


if __name__=='__main__':unittest.main()
