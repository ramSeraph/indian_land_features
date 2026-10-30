"""Audit native river-treatment interfaces, including sheet seams."""
import numpy as np


def audit_native_joins(pixels,before,after,applied,groups,width,sheets=None):
    """Inspect written raster values at treatment interfaces, including seams."""
    order=np.argsort(pixels)
    pixels,before,after,applied,groups=(v[order] for v in (pixels,before,after,applied,groups))
    if sheets is not None:sheets=sheets[order]
    pairs=0;worst=0.;maximum_signed_step_change=0.
    for offset in (1,width):
        if not len(pixels):continue
        j=np.searchsorted(pixels,pixels+offset);safe=np.minimum(j,len(j)-1)
        good=(j<len(j))&(pixels[safe]==pixels+offset)
        if offset==1:good&=pixels%width<width-1
        if sheets is not None:good&=sheets!=sheets[safe]
        good&=(before!=-255)&(before[safe]!=-255)&np.isfinite(before)&np.isfinite(before[safe])
        boundary=(applied!=applied[safe])|(applied&applied[safe]&(groups!=groups[safe]))
        i=np.flatnonzero(good&boundary);jj=j[i]
        pairs+=len(i)
        for side in (i,jj):
            changed=side[applied[side]]
            worst=max(worst,float(np.abs(after[changed].astype(float)-before[changed]).max(initial=0)))
        step_change=(after[jj].astype(float)-after[i])-(before[jj].astype(float)-before[i])
        maximum_signed_step_change=max(maximum_signed_step_change,float(np.abs(step_change).max(initial=0)))
    return dict(boundary_pairs=pairs,maximum_treated_boundary_correction_m=worst,
                maximum_signed_step_change_m=maximum_signed_step_change)
