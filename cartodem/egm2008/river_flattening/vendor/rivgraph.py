"""Extracted RivGraph functions; see LICENSE-RivGraph.txt and provenance.json.
Modified extraction: explicit tie rejection and minimum parallel-edge weight.
Copyright 2019 Triad National Security, LLC; BSD-3-Clause.
"""


import numpy as np
import networkx as nx


def algmap(key=None):
    """
    Returns a numeric key corresponding to each algorithm used by RivGraph to
    set link directionality. These numbers are found in links['guess_alg'] and
    links['certain_alg'] and are useful for diagnosing issues in setting
    directionality.

    Parameters
    ----------
    key : str
        See the mapper below for the possible keywords.

    Returns
    -------
    algno : float or int
        Numeric representation of direction-setting algorithm.
    """

    mapper = {
              'sourcesinkfix' : -2,
              'manual_set' : -1,
              'inletoutlet' : 0,
              'continuity' : 1,
              'parallels' : 2,
              'artificials' : 2.1,
              'main_chans' : 4,
              'bridges' : 5,
              'known_fdr' : 6,
              'known_fdr_rs' : 6.1,
              'syn_dem' : 10,
              'syn_dem_med' : 10.1,
              'sym_dem_leftover' : 10.2,
              'sp_links' : 11,
              'sp_links_sink' : -11.1,
              'sp_nodes' : 12,
              'longest_steepest' : 13,
              'three_agree' : 15,
              'syn_dem_and_sp' : 16,
              'cl_dist_guess' : 20,
              'cl_ang_guess' : 21,
              'cl_dist_set' : 22,
              'cl_ang_set' : 23,
              'cl_ang_rs' : 23.1,
              'cl_dist_and_ang' : 24,
              'short_no_bktrck' : 25,
              'wid_pctdiff' : 26,
              'wid_continuity' : 27,
              'wid_continuity_violation' : 27.1,
              'darea_grad': 28,
              'random': 29,
              'n_agree': 30.,  # decimal indicates number of agreement
              'main_channel_darea_grad': 31,
    }

    return mapper[key] if key else mapper


def dir_shortest_paths_links(links, nodes, difthresh=0):
    """
    Guesses a link's directionality by which of its endpoint nodes is
    closer to the nearest outlet (or preoutlet). "difthresh" refers to the
    difference in distances between endpoint nodes; higher means directionality
    is less likely to be ascertained.

    Parameters
    ----------
    links : dict
        Network links and associated properties.
    nodes : dict
        Network nodes and associated properties.
    difthresh : int, optional
        Threshold for the length (in number of nodes) that must be surpassed
        between an unknown link's endpoints for directionality to be set.
        Lengths are measured from the link's endpoints to the nearest
        outlet (or preoutlet). The default is 0.

    Returns
    -------
    links : dict
        Network links and associated properties updated to set directionality
        according to this algorithm.
    nodes : dict
        Network nodes and associated properties updated to set directionality
        according to this algorithm.
    """
    # alg = 11 # identifier for diagnostics
    alg = algmap('sp_links')

    # Create networkX graph, adding weighted edges
    G = nx.Graph()
    G.add_nodes_from(nodes['id'])
    for lc, wt in zip(links['conn'], links['len']):
        G.add_edge(lc[0], lc[1], weight=min(wt, G.get_edge_data(lc[0], lc[1], {}).get('weight', np.inf)))

    # Get all "pre-outlet", i.e. nodes one link upstream of outlets. Use these so that decision of where to chop off outlet links doesn't play a role in shortest path.
    preoutlets = []
    for o in nodes['outlets']:
        linkconn = links['conn'][links['id'].index(nodes['conn'][nodes['id'].index(o)][0])]
        othernode = linkconn[:]
        othernode.remove(o)
        preoutlets.append(othernode[0])

    # NetworkX's multi-source dijkstra path length iterator
    msdpl = nx.multi_source_dijkstra_path_length(G, preoutlets)

    # METHOD 2: Direction ascertained based on the nearness of a link's endnodes
    # to the nearest outlet [closer is downstream]
    for il,lid in enumerate(links['id']):

        linkidx = links['id'].index(lid)
        lconn = links['conn'][linkidx][:]

        # Compute shortest distance from each end node to nearest pre-outlet
        node_min_len = []
        for endnode in lconn:
            node_min_len.append(msdpl[endnode])

        # Skip if the minima are too similar
        if abs(node_min_len[0] - node_min_len[1]) <= difthresh:
            continue
        else:  # DS node is the closer one
            dsnode = lconn[node_min_len.index(min(node_min_len))]
            usnode = list(set(lconn) - set([dsnode]))[0]

            # Update certainty
            links['guess'][linkidx].append(usnode)
            links['guess_alg'][linkidx].append(alg)

    return links, nodes


def guess_synthetic_slope(links, nodes, node_column="synthetic_elevation",
                          alg=algmap('syn_dem'), ascending=False, min_slope=0, max_slope=1):
    """Guess directions by the slope of a synthetic elevation, e.g. distance to outlets.

    Links are not flipped, only guessed upstream node is attached to links[guess].
    """
    slopes = np.zeros(len(links["id"]), dtype=float)
    for linkidx, lid in enumerate(links['id']):
        conn = links["conn"][linkidx]
        elev = [nodes[node_column][nodes["id"].index(i)] for i in conn]
        slope = (elev[-1] - elev[0]) / (links["len"][linkidx] / links["sinuosity"][linkidx])
        # Make sure slope is negative, else flip direction
        di = slope < 0 if ascending else slope > 0
        if np.abs(slope) > min_slope and np.abs(slope) <= max_slope:
            usnode = conn[int(di)]
            # Store guess
            links['guess'][linkidx].append(usnode)
            links['guess_alg'][linkidx].append(alg)
        # Store slope
        slopes[linkidx] = slope * (-1 if di else 1)

    links[node_column+'_slope'] = slopes
    return links, nodes
