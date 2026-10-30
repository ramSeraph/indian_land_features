"""Reproduce the RivGraph extract from its pinned, public source revision.

Run: python -m river_flattening.vendor.refresh
The selected functions remain upstream code except for the explicitly listed
compatibility/safety edits below. Full upstream package dependencies are avoided.
"""
import ast
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

REVISION = '0a0a7135752894f451344bbf3521cc0f85f6f426'
BASE = f'https://raw.githubusercontent.com/mwort/RivGraph/{REVISION}/'
SELECTION = {
    'rivgraph/directionality.py': ['algmap', 'dir_shortest_paths_links'],
    'rivgraph/networks/network_directionality.py': ['guess_synthetic_slope'],
}


def generate(sources, license_text):
    result = ['"""Extracted RivGraph functions; see LICENSE-RivGraph.txt and provenance.json.\n'
              'Modified extraction: explicit tie rejection and minimum parallel-edge weight.\n'
              'Copyright 2019 Triad National Security, LLC; BSD-3-Clause.\n"""\n',
              'import numpy as np\nimport networkx as nx\n']
    for path, names in SELECTION.items():
        source = sources[path]
        functions = {node.name: ast.get_source_segment(source, node)
                     for node in ast.parse(source).body if isinstance(node, ast.FunctionDef)}
        for name in names:
            code = functions[name]
            code = code.replace("alg=dy.algmap('syn_dem')", "alg=algmap('syn_dem')")
            # A MultiGraph is not needed for distances, but the last inserted
            # parallel branch must not replace a shorter branch nondeterministically.
            code = code.replace('G.add_edge(lc[0], lc[1], weight=wt)',
                "G.add_edge(lc[0], lc[1], weight=min(wt, G.get_edge_data(lc[0], lc[1], {}).get('weight', np.inf)))")
            code = code.replace('abs(node_min_len[0] - node_min_len[1]) < difthresh',
                                'abs(node_min_len[0] - node_min_len[1]) <= difthresh')
            result.append(code+'\n')
    directory = Path(__file__).parent
    (directory/'rivgraph.py').write_text('\n\n'.join(result))
    (directory/'LICENSE-RivGraph.txt').write_text(license_text)
    provenance = dict(repository='https://github.com/mwort/RivGraph', revision=REVISION,
                      license='BSD-3-Clause', functions=SELECTION,
                      normalized_source_sha256={p: hashlib.sha256(s.encode()).hexdigest() for p,s in sources.items()},
                      edits=['Remove dy namespace from default algmap argument',
                             'Reject equal outlet distances rather than choose an arbitrary endpoint',
                             'Keep minimum parallel-edge length in the distance graph'])
    (directory/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')


if __name__ == '__main__':
    sources = {path: urlopen(BASE+path).read().decode().replace('\r\n','\n') for path in SELECTION}
    generate(sources, urlopen(BASE+'LICENSE.txt').read().decode().replace('\r\n','\n'))
