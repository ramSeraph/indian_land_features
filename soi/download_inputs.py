"""Download pinned public input assets, resuming transfers and checking SHA256.

Usage: python3 download_inputs.py manifests/soi-contours.json data/contours/soi-raw
The manifest is immutable input. A different upstream asset fails verification.
"""

import argparse
import hashlib
import json
import shutil
import subprocess
from datetime import UTC, datetime
from pathlib import Path


def checksum(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def download(manifest_path, destination, min_free_gib=2):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(Path(manifest_path).read_text())
    names = set()
    for asset in manifest["assets"]:
        name = asset["name"]
        if Path(name).name != name or name in names or name in ("", ".", ".."):
            raise ValueError(f"Unsafe or duplicate asset name: {name}")
        names.add(name)
        if not asset["url"].startswith("https://") or len(asset["sha256"]) != 64:
            raise ValueError("Expected HTTPS URL and SHA256 for every asset")
    receipt = {"manifest_sha256": checksum(manifest_path), "assets": []}
    for asset in manifest["assets"]:
        target = destination / asset["name"]
        if target.exists():
            if (
                target.stat().st_size != asset["bytes"]
                or checksum(target) != asset["sha256"]
            ):
                raise ValueError(f"Existing input differs from pinned asset: {target}")
            print(f"Verified cached {target.name}", flush=True)
        else:
            staged = target.with_name(target.name + ".part")
            remaining = max(
                0, asset["bytes"] - (staged.stat().st_size if staged.exists() else 0)
            )
            if shutil.disk_usage(destination).free < remaining + min_free_gib * 1024**3:
                raise RuntimeError(
                    f"Insufficient disk reserve to download {target.name}"
                )
            subprocess.run(
                [
                    "curl",
                    "--fail",
                    "--location",
                    "--show-error",
                    "--retry",
                    "5",
                    "--retry-delay",
                    "2",
                    "--connect-timeout",
                    "30",
                    "--continue-at",
                    "-",
                    "--output",
                    str(staged),
                    asset["url"],
                ],
                check=True,
            )
            if (
                staged.stat().st_size != asset["bytes"]
                or checksum(staged) != asset["sha256"]
            ):
                raise ValueError(f"Download differs from pinned asset: {staged}")
            staged.replace(target)
        receipt["assets"].append(
            {
                **asset,
                "path": str(target.resolve()),
                "mtime_ns": target.stat().st_mtime_ns,
            }
        )
    receipt["verified_at"] = datetime.now(UTC).isoformat()
    receipt_path = destination / "inputs_verified.json"
    staged_receipt = receipt_path.with_suffix(".json.tmp")
    staged_receipt.write_text(json.dumps(receipt, indent=2) + "\n")
    staged_receipt.replace(receipt_path)
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--min-free-gib", type=float, default=2)
    args = parser.parse_args()
    download(args.manifest, args.destination, args.min_free_gib)
