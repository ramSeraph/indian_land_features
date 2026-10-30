# CartoDEM noise correction: source review

Reviewed 2026-10-23. This is a proposed experiment, not a correction already
applied to the CartoDEM build. No land elevations have been denoised.

## Evidence specific to CartoDEM

Goyal, Featherstone, Dikshit and Balasubramania, *Comparison and Validation of
Satellite-Derived Digital Surface/Elevation Models over India*,
[DOI 10.1007/s12524-020-01273-7](https://doi.org/10.1007/s12524-020-01273-7),
was read from the supplied PDF, including Figure 6. PDF page numbers below count
the title page as page 1.

Supplied PDF SHA256:
`8d52c54ee857db58eb5a2f78f480fe36471cc3b77242ee4dda851b878893f0f0`.

- Page 7 and Figure 6 on page 8 identify stripe effects in CartoDEM through
  comparisons with SRTM and ASTER. Their conclusions also identify CartoDEM
  striping. This is direct motivation to inspect our release for directional
  errors, not proof that every linear feature is an artefact.
- Table 9 (page 11): over 3,842 benchmarks, CartoDEM MAE/RMSE are 8.5/19.7 m;
  MERIT's are 7.1/17.3 m. Over 145 GCPs, the corresponding figures are 5.7/9.6 m
  and 4.5/7.1 m.
- Table 10 (page 12): below 500 m, CartoDEM MAE/RMSE are 6.7/14.5 m versus
  MERIT's 5.5/11.7 m. Ground control is concentrated in lower elevations and
  geographically limited areas; the authors explicitly qualify mountain results.
- CartoDEM is treated as a DSM. The authors add modelled canopy heights to
  ground truth for DSM comparisons, but not for MERIT. They do not independently
  validate the canopy model. These numbers therefore do not isolate the benefit
  of destriping or prove that a denoised 30 m CartoDEM will match MERIT accuracy.
- The paper harmonizes vertical references. Its evaluation does not establish
  that the exact CartoDEM v3r1 sheets in our build have identical errors.

## Transferable methods

[Yamazaki et al. (2017), section 2.2](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017GL072874)
separates stripe, absolute, vegetation and short-scale random errors. It uses
2-D Fourier filtering for stripes and Gallant's adaptive smoothing for random
noise. MERIT also processes optical stereo AW3D, so applicability is not confined
to radar SRTM. Its analysis grid is about 90 m; CartoDEM needs separate calibration.

[Gallant (2011)](https://www.geomorphometry.org/uploads/pdf/pdf2011/Gallant2011geomorphometry.pdf)
uses multiple scales and an estimated noise variance: stronger smoothing where
noise dominates local relief, weaker smoothing where terrain variation dominates.
Its optional missing-data filling must be disabled for published CartoDEM output.

## Proposed CartoDEM experiment

1. Preserve both original nodata masks, normalize to `-255`, correct to EGM2008,
   and exclude ocean and already treated water surfaces from noise estimation.
2. Diagnose directional spectra over overlapping multi-sheet windows. Identify
   stripe orientation, wavelength and spatially varying amplitude from our data.
   Do not assume SRTM's parameters or erase an entire frequency direction.
3. Estimate only the supported stripe component and subtract it. Keep a raster
   of removed elevations. Evaluate sharp scene/strip steps separately: the H44N
   discontinuity is not established to be periodic striping. Blurring it into a
   ramp would not establish a corrected surface.
4. Test adaptive noise smoothing on representative plains, coast, hills and
   mountains, preserving original nodata and class-2 treatment. Retaining a 30 m
   output grid does not guarantee retention of all 30 m terrain detail.
5. Validate hillshade, slopes, ridge/valley profiles, boundary continuity and
   independent height evidence. Use withheld reference points where available;
   agreement with a second noisy DSM alone cannot prove an improvement. Resolve
   datum, alignment, canopy and observation-date differences before comparison.

The available **SOI DEM covering at least the Gangetic plains** is a candidate
independent reference for the first pilot (user confirmed 2026-10-24). After
checking horizontal alignment and recording the SOI vertical reference, inspect
CartoDEM-minus-SOI residuals for coherent directional patterns and broad biases.
SOI's contour-derived surface has its own interpolation errors and limited
effective detail; fine residuals must not automatically be labelled CartoDEM
noise. **No automatic SOI vertical conversion is required for this pilot.** Its
sea-level-based heights must not receive the ellipsoidal CartoDEM `h - N`
subtraction. A meaningful reference offset matters for absolute height blending;
for stripe diagnostics, separate a broad offset from the spatial residual
patterns without zero-centering every sheet independently. The SOI workflow
explicitly retains source heights (see `../soi/PRODUCTION.md` in the parent
repository). Use separate areas for fitting a correction and
checking its effect. No new national correction is selected by this review.

If accepted, noise correction should precede river-direction inference and bank
profile fitting, so those stages use the improved terrain. EGM2008 conversion
changes the vertical reference; it does not remove observational noise or canopy
bias. Final COGs and TerrainRGB follow the accepted elevation processing.

The existing GLO-30 role remains ocean/river masks only. This review does not add
GLO-30 heights or MERIT heights as replacements for CartoDEM.
