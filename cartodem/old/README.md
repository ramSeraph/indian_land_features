# Earlier CartoDEM scripts

These files preserve earlier workflows. The current entry point is
[../egm2008/README.md](../egm2008/README.md).

- `steps.sh` records the previous GDAL/rio-rgbify workflow and its zoom-dependent
  quantization. The original scraping, contour, naming and packaging scripts
  are retained beside it.
- `egm2008/` contains superseded conversion helpers, the low-zoom extension,
  an unused repacking approach, the evaluated rio-rgbify alternative, river
  experiments, the one-time release replacement helper, and historical method notes.

Historical instructions and artifact links describe their original workspace
layout and processing choices, including alpha 254 where explicitly recorded.
They are not the instructions for the current release. Shared numerical helpers
that the current pipeline imports remain in `egm2008/` at the repository's
CartoDEM level.
