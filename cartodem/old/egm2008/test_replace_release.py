import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from build import atomic_json, read_json
from replace_release import publish


class FakeGitHub:
    base = 'repos/owner/repo'

    def __init__(self, release, files):
        self.release_data = dict(release)
        self.files = {f['name']: dict(f) for f in files}
        self.mutations = []
        self.corrupt_upload = False

    def release(self, tag):
        return dict(self.release_data)

    def inventory(self, release_id):
        return dict(self.files)

    def request(self, endpoint, method='GET', fields=None, path=None, mime=None):
        self.mutations.append((method, endpoint))
        if method == 'DELETE':
            asset_id = int(endpoint.rsplit('/', 1)[1])
            self.files = {k: v for k, v in self.files.items() if v['id'] != asset_id}
            return None
        assert method == 'POST'
        from urllib.parse import parse_qs, urlparse
        name = parse_qs(urlparse(endpoint).query)['name'][0]
        item = dict(id=999, name=name, size=path.stat().st_size, state='uploaded',
                    digest='sha256:'+hashlib.sha256(path.read_bytes()).hexdigest())
        if self.corrupt_upload:
            item['digest'] = 'sha256:wrong'
        self.files[name] = item
        return item


class ReplacementTests(unittest.TestCase):
    def fixture(self, root):
        folder = root/'publication'
        folder.mkdir()
        (root/'cog.tif').write_bytes(b'new raster')
        entry = dict(path='cog.tif', name='cog.tif', bytes=10,
                     sha256=hashlib.sha256(b'new raster').hexdigest())
        release = dict(id=12, tag_name='tag', name='same title', body='unchanged\r\nnotes',
                       draft=False, prerelease=False, target_commitish='main',
                       upload_url='https://uploads.github.com/repos/owner/repo/releases/12/assets{?name,label}')
        old = dict(id=34, name='cog.tif', size=3, digest='sha256:old', state='uploaded')
        atomic_json(folder/'release-before.json', release)
        atomic_json(folder/'assets-before.json', [old])
        plan = dict(entries=[entry], assets=1, bytes=10, tag='tag', release_id=12)
        return folder, plan, FakeGitHub(release, [old])

    def test_replace_verifies_upload_preserves_notes_and_resumes(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            folder, plan, gh = self.fixture(root)
            before = gh.release_data.copy()
            publish(root, folder, gh, plan)
            self.assertEqual([m[0] for m in gh.mutations], ['DELETE', 'POST'])
            self.assertEqual(gh.release_data, before)
            self.assertEqual(read_json(folder/'status.json')['state'], 'complete')
            gh.mutations.clear()
            publish(root, folder, gh, plan)
            self.assertEqual(gh.mutations, [])

    def test_corrupt_local_input_never_deletes_old_asset(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            folder, plan, gh = self.fixture(root)
            (root/'cog.tif').write_bytes(b'bad')
            with self.assertRaisesRegex(ValueError, 'Prepared asset changed'):
                publish(root, folder, gh, plan)
            self.assertEqual(gh.mutations, [])
            self.assertEqual(gh.files['cog.tif']['id'], 34)

    def test_external_replacement_and_note_edit_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            folder, plan, gh = self.fixture(root)
            gh.files['cog.tif']['id'] = 35
            with self.assertRaisesRegex(ValueError, 'changed outside'):
                publish(root, folder, gh, plan)
            self.assertEqual(gh.mutations, [])
            gh.files['cog.tif']['id'] = 34
            gh.release_data['body'] += '!'
            with self.assertRaisesRegex(ValueError, 'metadata changed'):
                publish(root, folder, gh, plan)
            self.assertEqual(gh.mutations, [])

    def test_bad_upload_digest_is_not_reported_as_success(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            folder, plan, gh = self.fixture(root)
            gh.corrupt_upload = True
            with self.assertRaisesRegex(ValueError, 'digest or size differs'):
                publish(root, folder, gh, plan)
            self.assertEqual(read_json(folder/'status.json')['state'], 'failed')


if __name__ == '__main__':
    unittest.main()
