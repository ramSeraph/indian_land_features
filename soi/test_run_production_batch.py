"""Check upload retry and stop behaviour without contacting GitHub."""

import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from run_production_batch import run_producer


class UploadRetryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.path = self.root / "state.json"
        self.tick = 1_000_000_000

    def checkpoint(self, *, error_type="CalledProcessError", sheet="44K_6"):
        self.path.write_text(
            json.dumps(
                {
                    "status": "error",
                    "error_type": error_type,
                    "current_sheet": sheet,
                    "sheets": {
                        sheet: {"status": "ready", "outcome": {"status": "built"}}
                    },
                }
            )
        )
        self.tick += 1_000_000_000
        os.utime(self.path, ns=(self.tick, self.tick))

    def test_upload_failure_resumes_then_succeeds(self):
        def execute(*args, **kwargs):
            if runner.call_count == 1:
                self.checkpoint()
                raise subprocess.CalledProcessError(1, ["producer"])

        with (
            patch("run_production_batch.subprocess.run", side_effect=execute) as runner,
            patch("run_production_batch.time.sleep") as sleep,
        ):
            run_producer(["producer"], self.root, self.root)
        self.assertEqual(runner.call_count, 2)
        sleep.assert_called_once_with(15)

    def test_retries_are_bounded(self):
        def execute(*args, **kwargs):
            self.checkpoint()
            raise subprocess.CalledProcessError(1, ["producer"])

        with (
            patch("run_production_batch.subprocess.run", side_effect=execute) as runner,
            patch("run_production_batch.time.sleep") as sleep,
        ):
            with self.assertRaises(subprocess.CalledProcessError):
                run_producer(["producer"], self.root, self.root)
        self.assertEqual(runner.call_count, 5)
        self.assertEqual([call.args[0] for call in sleep.call_args_list], [15, 30, 60, 120])

    def test_old_upload_error_does_not_hide_startup_failure(self):
        self.checkpoint()
        with (
            patch(
                "run_production_batch.subprocess.run",
                side_effect=subprocess.CalledProcessError(1, ["producer"]),
            ) as runner,
            patch("run_production_batch.time.sleep") as sleep,
        ):
            with self.assertRaises(subprocess.CalledProcessError):
                run_producer(["producer"], self.root, self.root)
        self.assertEqual(runner.call_count, 1)
        sleep.assert_not_called()

    def test_verification_mismatch_is_not_retried(self):
        def execute(*args, **kwargs):
            self.checkpoint(error_type="ValueError")
            raise subprocess.CalledProcessError(1, ["producer"])

        with (
            patch("run_production_batch.subprocess.run", side_effect=execute) as runner,
            patch("run_production_batch.time.sleep") as sleep,
        ):
            with self.assertRaises(subprocess.CalledProcessError):
                run_producer(["producer"], self.root, self.root)
        self.assertEqual(runner.call_count, 1)
        sleep.assert_not_called()


if __name__ == "__main__":
    unittest.main()
