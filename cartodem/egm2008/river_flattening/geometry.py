"""Native-grid river geometry and cross-channel assignment for a bounded window.

The water footprint is never enlarged. Skeleton pixels are mapped to the common
mosaic's pixel IDs, allowing geometry from overlapping windows to be joined.
Water pixels use the nearest skeleton point measured *through the river mask*;
a Euclidean shortcut across a meander or island is not allowed.
"""
import numpy as np
from scipy import ndimage, sparse
from scipy.sparse.csgraph import dijkstra, connected_components
from scipy.sparse.linalg import spsolve
from scipy.spatial import cKDTree
from skimage.morphology import skeletonize


def pixel_graph(ids, width, dx, dy, *, remove_corner_diagonals=False):
    """8-neighbour sparse graph on sorted unique linear pixel IDs."""
    ids = np.asarray(ids, dtype=np.int64)
    if len(ids) == 0 or np.any(np.diff(ids) <= 0):
        raise ValueError('Expected sorted unique pixel IDs')
    dx=np.asarray(dx)
    if dx.ndim>1 or (dx.ndim==1 and len(dx)!=len(ids)) or np.any(dx<=0) or dy<=0:
        raise ValueError('Invalid pixel ground distances')
    rows, cols, weights = [], [], []
    for dr, dc in ((0, 1), (1, -1), (1, 0), (1, 1)):
        candidate = ids + dr*width + dc
        j = np.searchsorted(ids, candidate)
        safe = np.minimum(j, len(ids)-1)
        good = (j < len(ids)) & (ids[safe] == candidate)
        good &= (ids % width + dc >= 0) & (ids % width + dc < width)
        if remove_corner_diagonals and dr and dc:
            # A redundant diagonal around an L corner creates a false tiny cycle.
            for offset in (dr*width, dc):
                k = np.searchsorted(ids, ids+offset)
                good &= ~((k < len(ids)) & (ids[np.minimum(k, len(ids)-1)] == ids+offset))
        i = np.flatnonzero(good)
        rows.extend([i, j[good]])
        cols.extend([j[good], i])
        local_dx=dx if dx.ndim==0 else (dx[i]+dx[j[good]])/2
        steps=np.broadcast_to(np.hypot(dc*local_dx,dr*dy),i.shape).copy()
        weights.extend([steps]*2)
    return sparse.csr_matrix((np.concatenate(weights),
                              (np.concatenate(rows), np.concatenate(cols))), shape=(len(ids), len(ids)))


def extract_window(heights, classes, *, dx, dy, grid_width, row_offset=0,
                   col_offset=0, core=None, nodata=-255.0, bank_pixels=2, allocator=None):
    """Extract geometry using a halo, returning only core-owned output/samples.

    Heights/classes must already share the original CartoDEM grid. `core` is
    (row_start, row_end, col_start, col_end) in this window. Nonfinite and declared
    nodata elevations cannot be bank observations. Lake class 2 is not sampled.
    An allocator returns local seed pixel IDs in row-major river-cell order.
    The default SciPy implementation is retained for small reference studies;
    production supplies the memory-bounded Whitebox raster allocator.
    """
    heights, classes = np.asarray(heights), np.asarray(classes)
    if heights.ndim != 2 or classes.shape != heights.shape or dx <= 0 or dy <= 0:
        raise ValueError('Invalid aligned river inputs')
    height, width = heights.shape
    if core is None:
        core = (0, height, 0, width)
    r0, r1, c0, c1 = core
    if not (0 <= r0 < r1 <= height and 0 <= c0 < c1 <= width):
        raise ValueError('Invalid core window')
    river = classes == 3
    flat = np.flatnonzero(river).astype(np.int64)
    if len(flat) == 0:
        return {k:np.empty(0, dtype=d) for k,d in
                [('nodes','int64'),('pixels','int64'),('owners','int64'),
                 ('bank_owners','int64'),('bank_heights','float32'),('ocean_owners','int64')]}
    skeleton = skeletonize(river)
    skeleton,pruned=prune_spurs(skeleton,river,dx,dy)
    skel_flat = np.flatnonzero(skeleton)
    if allocator is None:
        water_graph = pixel_graph(flat, width, dx, dy)
        sources = np.searchsorted(flat, skel_flat)
        distance, _, owner = dijkstra(water_graph, directed=False, indices=sources,
                                      min_only=True, return_predecessors=True)
        if not np.isfinite(distance).all() or np.any(owner < 0):
            raise ValueError('A river pixel has no connected skeleton owner')
        local_owner = flat[owner]
        del water_graph, distance, owner
    else:
        local_owner = allocator(river, skeleton, dx, dy)
        if local_owner.shape != flat.shape or not np.isin(local_owner, skel_flat).all():
            raise ValueError('Invalid river allocation')

    def global_ids(local):
        return (local//width + row_offset)*grid_width + local%width + col_offset

    def in_core(local):
        r, c = local//width, local%width
        return (r >= r0) & (r < r1) & (c >= c0) & (c < c1)

    global_owner = global_ids(local_owner)
    del local_owner
    core_water = in_core(flat)
    bank = ndimage.binary_dilation(river, iterations=bank_pixels)
    bank &= (classes == 0) & np.isfinite(heights) & (heights != nodata)
    bank_flat = np.flatnonzero(bank)
    bank_flat = bank_flat[in_core(bank_flat)]
    # Bank samples are outside water; enter via their nearest river-edge pixel,
    # then use that pixel's in-water assignment to the skeleton.
    points = np.column_stack((flat%width*dx, flat//width*dy))
    tree = cKDTree(points)
    _, nearest_water = tree.query(np.column_stack((bank_flat%width*dx, bank_flat//width*dy)))
    near_ocean = ndimage.binary_dilation(classes == 1)[river]
    return dict(nodes=global_ids(skel_flat[in_core(skel_flat)]),
                pixels=global_ids(flat[core_water]), owners=global_owner[core_water],
                bank_owners=global_owner[nearest_water], bank_heights=heights.ravel()[bank_flat].astype(np.float32),
                ocean_owners=np.unique(global_owner[near_ocean & core_water]),
                pruned_skeleton_pixels=pruned)


def prune_spurs(skeleton,river,dx,dy):
    """Remove short medial-axis branches toward bank corners, preserving ports.

    Only a terminal-to-junction branch can be pruned. Standalone channels and
    endpoints at the processing-window boundary are retained. The threshold is
    the greater of 90 m and twice the junction's distance to the river edge.
    This changes the centreline representation, never the GLO-30 footprint.
    """
    skeleton=skeleton.copy()
    radii=ndimage.distance_transform_edt(river,sampling=(dy,dx))
    height,width=river.shape
    original=int(skeleton.sum())
    for _ in range(4):
        ids=np.flatnonzero(skeleton)
        graph=pixel_graph(ids,width,dx,dy,remove_corner_diagonals=True)
        degree=np.diff(graph.indptr)
        remove=[]
        for leaf in np.flatnonzero(degree==1):
            pixel=ids[leaf]
            row,col=divmod(int(pixel),width)
            if min(row,col,height-1-row,width-1-col)<=1:
                continue
            path=[int(leaf)]
            previous=-1; current=int(leaf); length=0.
            while True:
                a,b=graph.indptr[current:current+2]
                options=graph.indices[a:b]
                options=options[options!=previous]
                if not len(options):
                    break
                nxt=int(options[0])
                index=a+np.searchsorted(graph.indices[a:b],nxt)
                length+=graph.data[index]
                if degree[nxt]!=2:
                    if degree[nxt]>=3 and length<max(90.,2*radii.ravel()[ids[nxt]]):
                        remove.extend(path)
                    break
                path.append(nxt)
                previous,current=current,nxt
        if not remove:
            break
        skeleton.ravel()[ids[np.asarray(remove)]]=False
    return skeleton,original-int(skeleton.sum())


def compress_network(ids, width, dx, dy, *, spacing=180.0):
    """Retain junctions, sample degree-two chains, and keep native projections.

    Returns shared stations/reaches plus left/right station IDs and interpolation
    fractions for every original skeleton point. Pure closed loops are reported
    as unsupported; no arbitrary flow direction is introduced for them.
    """
    ids = np.unique(np.asarray(ids, dtype=np.int64))
    graph = pixel_graph(ids, width, dx, dy, remove_corner_diagonals=True)
    degree = np.diff(graph.indptr)
    special = degree != 2
    high = np.flatnonzero(degree > 2)
    # Collapse neighbouring junction pixels into one water-level station.
    station = np.full(len(ids), -1, dtype=np.int32)
    representatives = []
    if len(high):
        ng, labels = connected_components(graph[high][:,high], directed=False)
        station[high] = labels
        for group in range(ng):
            representatives.append(int(high[np.flatnonzero(labels == group)[0]]))
    for i in np.flatnonzero(degree <= 1):
        station[i] = len(representatives)
        representatives.append(int(i))
    left, right = station.copy(), station.copy()
    fraction = np.zeros(len(ids), dtype=np.float32)
    seen = np.zeros(len(graph.indices), dtype=bool)
    reaches, distances = [], []

    def neighbours(u):
        a,b = graph.indptr[u:u+2]
        return graph.indices[a:b]

    def mark(u,v):
        a,b = graph.indptr[u:u+2]
        k = np.searchsorted(graph.indices[a:b], v)
        seen[a+k] = True
        return float(graph.data[a+k])

    for start in np.flatnonzero(special):
        for next_node in neighbours(start):
            a,b = graph.indptr[start:start+2]
            if seen[a+np.searchsorted(graph.indices[a:b], next_node)]:
                continue
            if special[next_node] and station[start] == station[next_node]:
                mark(start,next_node); mark(next_node,start)
                continue
            path, x = [int(start)], [0.0]
            previous, current = int(start), int(next_node)
            while True:
                step = mark(previous,current)
                mark(current,previous)
                path.append(current); x.append(x[-1]+step)
                if special[current]:
                    break
                options = neighbours(current)
                following = int(options[0] if options[1] == previous else options[1])
                previous,current = current,following
            if station[path[0]] == station[path[-1]]:
                # A branch leaves and returns to one junction. Its direction
                # needs additional evidence; retain its native values for now.
                continue
            picks = [0]
            for k in range(1,len(path)-1):
                if x[k]-x[picks[-1]] >= spacing:
                    picks.append(k)
            if len(path)>2 and len(picks)==1:
                # Keep parallel branches distinct in the reduced graph.
                picks.append(len(path)//2)
            picks.append(len(path)-1)
            for k in picks[1:-1]:
                station[path[k]] = len(representatives)
                representatives.append(path[k])
            reach = [int(station[path[k]]) for k in picks]
            distances.append(np.asarray([x[k] for k in picks]))
            reaches.append(reach)
            for a,b in zip(picks[:-1],picks[1:]):
                for k in range(a,b+1):
                    i = path[k]
                    left[i],right[i] = station[path[a]],station[path[b]]
                    fraction[i] = (x[k]-x[a])/(x[b]-x[a])
    # A pure closed loop has no representatives; its empty index must still
    # be integer so the caller can retain the component unchanged.
    return dict(ids=ids, station_pixels=ids[np.asarray(representatives,dtype=np.int64)],
                reaches=reaches, distances=distances, left=left, right=right,
                fraction=fraction, unsupported_skeleton_pixels=int(np.sum((left<0)|(right<0))))


def orient_component(n, reaches, distances, root):
    """Geometric direction candidate for a component with an established outlet.

    Interior harmonic potential is computed with outlet=0 and terminal inlets=1.
    Unlike shortest paths, this can orient both arms of a split/rejoined channel
    consistently. It is a hypothesis, not measured flow: CartoDEM observations
    must independently support it before the surface is applied.
    """
    edge_lengths = {}
    for reach, x in zip(reaches, distances):
        for u,v,length in zip(reach[:-1], reach[1:], np.diff(x)):
            key = tuple(sorted((int(u), int(v))))
            edge_lengths[key] = min(edge_lengths.get(key, np.inf),float(length))
    edges = np.asarray(list(edge_lengths),dtype=int).reshape(-1,2)
    if not len(edges) or root not in edges:
        raise ValueError('Outlet is not in this river component')
    conductance = np.array([1/edge_lengths[tuple(e)] for e in edges])
    adjacency = sparse.csr_matrix((np.tile(conductance,2),
                                   (np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(n,n))
    active = np.flatnonzero(np.diff(adjacency.indptr)>0)
    n_components, _ = connected_components(adjacency[active][:,active],directed=False)
    if n_components != 1:
        raise ValueError('Fit one connected river component at a time')
    degree = np.diff(adjacency.indptr)
    terminals = np.flatnonzero(degree==1)
    if root not in terminals or len(terminals)<2:
        raise ValueError('Outlet/inlets are not unambiguous terminal stations')
    fixed = terminals
    potential = np.full(n,np.nan)
    potential[fixed] = 1
    potential[root] = 0
    free = np.setdiff1d(active,fixed)
    laplacian = sparse.diags(np.asarray(adjacency.sum(axis=1)).ravel())-adjacency
    if len(free):
        potential[free] = spsolve(laplacian.tocsr()[free][:,free],
                                  -laplacian.tocsr()[free][:,fixed]@potential[fixed])
    directed_reaches, directed_distances = [], []
    for reach,x in zip(reaches,distances):
        delta = potential[reach[0]]-potential[reach[-1]]
        if not np.isfinite(delta) or abs(delta)<1e-9:
            raise ValueError('A channel direction is unresolved by outlet geometry')
        if delta>0:
            directed_reaches.append(list(reach)); directed_distances.append(np.asarray(x))
        else:
            directed_reaches.append(list(reversed(reach)))
            directed_distances.append(x[-1]-np.asarray(x)[::-1])
    directed_edges = np.unique(np.asarray([e for r in directed_reaches for e in zip(r[:-1],r[1:])]),axis=0)
    return directed_edges,directed_reaches,directed_distances,potential
