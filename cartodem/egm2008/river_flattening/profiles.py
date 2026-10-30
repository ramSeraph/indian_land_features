"""Fit a shared monotone river profile from CartoDEM observations.

This module does not infer flow direction. Edges must be supplied upstream to
downstream. The independently implemented objective follows the constrained
quantile-smoothing formulation described by Schwanghart & Scherler (2017):
https://doi.org/10.5194/esurf-5-821-2017

Coordinates and distances are in metres; observed elevations are EGM2008 metres.
Junctions are represented by a single station index shared by adjoining reaches.
"""
from collections import deque

import numpy as np
import clarabel
from scipy import sparse
from scipy.sparse.csgraph import connected_components


def topological_order(n, edges):
    """Reject cycles, including an accidentally reversed tributary/delta branch."""
    edges = np.asarray(edges, dtype=np.int64).reshape(-1, 2)
    if np.any((edges < 0) | (edges >= n)) or np.any(edges[:, 0] == edges[:, 1]):
        raise ValueError('Invalid directed river edge')
    if len(np.unique(edges, axis=0)) != len(edges):
        raise ValueError('Duplicate directed river edge')
    adjacency = sparse.csr_matrix((np.ones(len(edges)), (edges[:, 0], edges[:, 1])), shape=(n, n))
    indegree = np.bincount(edges[:, 1], minlength=n)
    pending = deque(np.flatnonzero(indegree == 0).tolist())
    order = []
    while pending:
        u = pending.popleft()
        order.append(u)
        for v in adjacency.indices[adjacency.indptr[u]:adjacency.indptr[u+1]]:
            indegree[v] -= 1
            if indegree[v] == 0:
                pending.append(int(v))
    if len(order) != n:
        raise ValueError('River directions contain a cycle; resolve it before fitting')
    return np.asarray(order), adjacency


def curvature_matrix(n, reaches, distances):
    """Second derivatives along each reach, with unequal station spacing."""
    rows, cols, values = [], [], []
    spacings = []
    row = 0
    if len(reaches) != len(distances):
        raise ValueError('Every reach requires station distances')
    for ids, x in zip(reaches, distances):
        ids, x = np.asarray(ids, dtype=np.int64), np.asarray(x, dtype=float)
        if (len(ids) != len(x) or len(ids) < 2 or np.any((ids < 0) | (ids >= n))
                or not np.isfinite(x).all() or np.any(np.diff(x) <= 0)):
            raise ValueError('Invalid reach stations or distances')
        spacings.extend(np.diff(x))
        for k in range(1, len(ids)-1):
            a, b = x[k]-x[k-1], x[k+1]-x[k]
            rows.extend([row]*3)
            cols.extend(ids[k-1:k+2])
            values.extend([2/(a*(a+b)), -2/(a*b), 2/(b*(a+b))])
            row += 1
    c = sparse.csc_matrix((values, (rows, cols)), shape=(row, n))
    spacing = float(np.median(spacings)) if spacings else 1.0
    return c, spacing


def fit_network(n, edges, observation_station, observations, *, reaches=(),
                distances=(), quantile=0.2, stiffness=2.0, weights=None,
                upper_bounds=None, max_iterations=200,
                unoriented_reaches=(), unoriented_distances=()):
    """Return Float32 station heights and numerical validation results.

    Missing observations must be excluded by the caller, using source validity.
    No interpolation is allowed for an entirely unobserved connected component.
    Optional upper bounds constrain station heights; they are not invented from
    a reference DEM. They must use the same datum as the observations.
    """
    if not isinstance(n, (int, np.integer)) or n < 1:
        raise ValueError('At least one river station is required')
    if not 0 < quantile < 1 or not np.isfinite(stiffness) or stiffness < 0:
        raise ValueError('Invalid quantile or smoothing stiffness')
    edges = np.asarray(edges, dtype=np.int64).reshape(-1, 2)
    ids = np.asarray(observation_station, dtype=np.int64)
    z = np.asarray(observations, dtype=np.float64)
    if (ids.ndim != 1 or z.ndim != 1 or len(ids) != len(z) or len(z) == 0
            or not np.isfinite(z).all() or np.any((ids < 0) | (ids >= n))):
        raise ValueError('Invalid bank observations; exclude nodata first')
    order, graph = topological_order(n, edges)
    geometry_edges=list(map(tuple,edges))
    for reach in unoriented_reaches:
        geometry_edges.extend(zip(reach[:-1],reach[1:]))
    geometry_edges=np.asarray(geometry_edges,dtype=int).reshape(-1,2)
    geometry_graph=sparse.csr_matrix((np.ones(len(geometry_edges)),
                                     (geometry_edges[:,0],geometry_edges[:,1])),shape=(n,n))
    _, labels = connected_components(geometry_graph, directed=False)
    if len(np.unique(labels[ids])) != len(np.unique(labels)):
        raise ValueError('A river component has no valid bank observations')
    m = len(z)
    if weights is None:
        # Equal total weight per sampled station, regardless of bank pixel count.
        weights = 1.0 / np.bincount(ids, minlength=n)[ids]
    weights = np.asarray(weights, dtype=float)
    if weights.shape != z.shape or not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError('Observation weights must be finite and positive')

    c, spacing = curvature_matrix(n, list(reaches)+list(unoriented_reaches),
                                  list(distances)+list(unoriented_distances))
    edge_set = set(map(tuple, edges.tolist()))
    for reach in reaches:
        if any((int(u), int(v)) not in edge_set for u, v in zip(reach[:-1], reach[1:])):
            raise ValueError('A reach disagrees with the directed network')
    scale = spacing**2 * stiffness * np.sqrt(n/c.shape[0]) if c.shape[0] else 0.0
    p_height = 2 * (scale*c).T @ (scale*c)
    # Unknowns: station heights, positive residuals, negative residuals.
    p = sparse.block_diag((p_height, sparse.csc_matrix((2*m, 2*m))), format='csc')
    q = np.concatenate([np.zeros(n), quantile*weights, (1-quantile)*weights])
    sample = sparse.csc_matrix((np.ones(m), (np.arange(m), ids)), shape=(m, n))
    identity = sparse.eye(m, format='csc')
    equality = sparse.hstack((sample, identity, -identity), format='csc')
    residual_bounds = sparse.hstack((sparse.csc_matrix((2*m, n)), sparse.eye(2*m)), format='csc')
    gradient = sparse.csc_matrix((np.tile([1., -1.], len(edges)),
                                  (np.repeat(np.arange(len(edges)), 2), edges.ravel())),
                                 shape=(len(edges), n))
    monotonic = sparse.hstack((gradient, sparse.csc_matrix((len(edges), 2*m))), format='csc')
    offset = float(np.median(z))
    matrices = [equality, -residual_bounds, -monotonic]
    bounds = [z-offset, np.zeros(2*m), np.zeros(len(edges))]
    cones = [clarabel.ZeroConeT(m), clarabel.NonnegativeConeT(2*m+len(edges))]
    if upper_bounds is not None:
        upper_bounds = np.asarray(upper_bounds, dtype=float)
        if upper_bounds.shape != (n,) or np.any(np.isnan(upper_bounds)) or np.any(upper_bounds == -np.inf):
            raise ValueError('Invalid bank upper bounds')
        capped=np.flatnonzero(np.isfinite(upper_bounds))
        selector=sparse.csc_matrix((np.ones(len(capped)),(np.arange(len(capped)),capped)),shape=(len(capped),n+2*m))
        matrices.append(selector)
        bounds.append(upper_bounds[capped]-offset)
        cones.append(clarabel.NonnegativeConeT(len(capped)))

    settings=clarabel.DefaultSettings()
    settings.verbose=False
    settings.max_iter=max_iterations
    settings.tol_gap_abs=1e-8
    settings.tol_gap_rel=1e-8
    settings.tol_feas=1e-9
    settings.max_threads=1
    solver=clarabel.DefaultSolver(sparse.triu(p,format='csc'),q,
                                  sparse.vstack(matrices,format='csc'),np.concatenate(bounds),cones,settings)
    solution=solver.solve()
    if solution.status != clarabel.SolverStatus.Solved:
        raise ValueError(f'River profile fit did not converge: {solution.status}')
    h = np.asarray(solution.x[:n]) + offset
    violation = float(np.maximum(0, h[edges[:, 1]]-h[edges[:, 0]]).max(initial=0))
    cap_violation = float(np.maximum(0, h-upper_bounds).max(initial=0)) if upper_bounds is not None else 0.0
    if not np.isfinite(h).all() or max(violation, cap_violation) > 0.001:
        raise ValueError('Fitted profile violates constraints beyond solver tolerance')
    original = h.copy()
    if upper_bounds is not None:
        h = np.minimum(h, upper_bounds)
    # Remove only numerical solver tolerance, so Float32 output is exactly monotone.
    for u in order:
        downstream = graph.indices[graph.indptr[u]:graph.indptr[u+1]]
        h[downstream] = np.minimum(h[downstream], h[u])
    numerical_adjustment = float(np.abs(original-h).max())
    if numerical_adjustment > 0.001:
        raise ValueError('Accumulated solver tolerance exceeds one millimetre')
    result = h.astype(np.float32)
    if np.any(result == -255) or not np.isfinite(result).all():
        raise ValueError('Fitted river height collides with output nodata or Float32 range')
    if np.any(result[edges[:, 1]] > result[edges[:, 0]]):
        raise ValueError('Float32 river profile is not monotone')
    residuals = z-result[ids]
    report = dict(stations=n, observations=m, edges=len(edges), quantile=quantile,
                  stiffness=stiffness, median_station_spacing_m=spacing,
                  iterations=int(solution.iterations), solver='Clarabel', solver_status='solved',
                  maximum_solver_monotonicity_error_m=violation,
                  maximum_numerical_adjustment_m=numerical_adjustment,
                  median_observation_minus_surface_m=float(np.median(residuals)),
                  p95_absolute_observation_residual_m=float(np.percentile(np.abs(residuals), 95)))
    return result, report
