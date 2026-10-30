import json

from pathlib import Path


with open('data/v3_r1_water_bodies.geojsonl', 'w') as of:
    for p in Path('data/water/v3_r1/').glob('*.geojsonl'):
        name = p.name
        tname = name.replace('.geojsonl', '')
        with open(p, 'r') as f:
            for line in f:
                feat = json.loads(line)
                feat['properties']['tile'] = tname
                of.write(json.dumps(feat))
                of.write('\n')
