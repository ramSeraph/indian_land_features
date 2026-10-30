# Recorded build provenance

These compact records describe the completed national COG build and the later
opaque TerrainRGB rebuild. Downloaded rasters and working databases are not
stored in Git.

- [inputs.lock.json](inputs.lock.json): URLs, sizes and checksums for every
  CartoDEM, geoid and GLO-30 WBM input.
- [cogs.json](cogs.json): the 543 final COG filenames, sizes, hashes, nodata and
  river-treatment counts. Pass this to `rebuild_tiles.py` when using those COGs.
- [river-results.json](river-results.json): actual fitting settings, treatment
  coverage, native checks, the maximum added adjacent-water step, and hashes of
  the full runtime reports from which this summary was extracted.
- [whitebox-installation.json](whitebox-installation.json): the pinned packages,
  binary/runtime hashes, settings and license information used for the build.
- [display-results.json](display-results.json): the opaque tile counts, four
  standard-tool partitions, and successful mosaic lookup check.
- [tests.json](tests.json): the original build's 65-test result in an isolated
  environment with the pinned processing dependencies. Tests for removed
  experiments are no longer part of the current suite.

The full per-sheet diagnostics are generated in the build workspace. Source
code and runtime provenance can affect container hashes even when native pixel
values agree; these checksums identify the recorded artifacts. The method and
its experimental limitations are documented in
[../river_flattening/README.md](../river_flattening/README.md).
