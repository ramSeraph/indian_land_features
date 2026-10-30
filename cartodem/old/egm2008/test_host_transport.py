import runpy
import shlex
import unittest
from pathlib import Path

WRAPPER = runpy.run_path(str(Path(__file__).with_name('build-host')))


class HostTransportTests(unittest.TestCase):
    def setUp(self):
        self.config = dict(host='builder@example.test', key='/tmp/a key',
                           root='/build workspace', code='code', venv='.venv')

    def test_run_quotes_every_argument_and_preserves_stdin(self):
        args = ['python', '-c', 'print("$HOME; `false`")', 'a b', "quote'foo"]
        command = WRAPPER['run_command'](self.config, args)
        self.assertIn('-T', command)
        self.assertEqual(command[-2], self.config['host'])
        self.assertEqual(shlex.split(command[-1].split(' && exec ', 1)[1]), args)
        self.assertIn("cd '/build workspace/code'", command[-1])

    def test_tunnel_binds_only_loopback_and_fails_on_port_conflict(self):
        command=WRAPPER['tunnel_command'](self.config,8767,8768,True)
        self.assertIn('127.0.0.1:8767:127.0.0.1:8768',command)
        self.assertIn('ExitOnForwardFailure=yes',command)
        self.assertIn('-N',command);self.assertIn('-f',command)
        self.assertEqual(command[-1],self.config['host'])
        with self.assertRaises(ValueError):WRAPPER['tunnel_command'](self.config,0,8768)

    def test_transfers_are_rooted_and_sync_excludes_private_data(self):
        command = WRAPPER['transfer_command'](self.config, 'sync', [], '')
        self.assertIn('--exclude=.local/', command)
        self.assertIn('--exclude=*', command)
        self.assertNotIn('--delete', command)
        self.assertEqual(command[-1], "builder@example.test:'/build workspace/code/'")
        for path in ('../escape', '/absolute', 'code/../../escape'):
            with self.assertRaisesRegex(ValueError, 'configured build root'):
                WRAPPER['remote_path'](self.config, path)
        with self.assertRaisesRegex(ValueError, 'Private connection'):
            WRAPPER['transfer_command'](self.config, 'push',
                [str(Path(__file__).parent/'.local/build-host.json')], 'code/')


if __name__ == '__main__':unittest.main()
