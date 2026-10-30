"""Small temporary workspaces using the same Whitebox binary as the pipeline."""
import os
from pathlib import Path
import tempfile

from build import setup
from river_flattening.allocation import allocation_recipe, make_allocator


def attach_whitebox(test, root):
    directory = os.environ.get('CARTODEM_WHITEBOX_DIR')
    if not directory:
        test.skipTest('Set CARTODEM_WHITEBOX_DIR to the installed tools/whitebox-prebuilt directory')
    source = Path(directory).expanduser().resolve()
    if not (source/'build.json').is_file():
        raise ValueError('CARTODEM_WHITEBOX_DIR has no Whitebox installation record')
    (root/'tools').mkdir()
    (root/'tools/whitebox-prebuilt').symlink_to(source, target_is_directory=True)
    return root


def whitebox_workspace(test):
    temporary = tempfile.TemporaryDirectory()
    test.addCleanup(temporary.cleanup)
    return attach_whitebox(test, setup(Path(temporary.name)/'work'))


def allocator(test):
    root = whitebox_workspace(test)
    work = root/'rivers/check'
    work.mkdir(parents=True)
    return make_allocator(root, work, 'fixture', dict(allocation=allocation_recipe(root)))
