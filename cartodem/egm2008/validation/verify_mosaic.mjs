// Verify every packaged tile is a candidate in the production mosaic resolver.
// Usage: node validation/verify_mosaic.mjs /path/to/pmtiles-mosaic/dist mosaic.json partition_info.json
import fs from 'node:fs';
import path from 'node:path';
import { pathToFileURL } from 'node:url';

const [moduleDir, mosaicPath, partitionPath] = process.argv.slice(2);
if (!moduleDir || !mosaicPath || !partitionPath) {
  throw new Error('Supply the pmtiles-mosaic dist directory, mosaic, and partition metadata');
}
const { MosaicResolver } = await import(pathToFileURL(path.resolve(moduleDir, 'index.js')));
const { createPartitionAwareIndex } = await import(
  pathToFileURL(path.resolve(moduleDir, 'indexers/partition-aware.js')));
const mosaic = JSON.parse(fs.readFileSync(mosaicPath, 'utf8'));
const partitions = JSON.parse(fs.readFileSync(partitionPath, 'utf8'));
const resolver = new MosaicResolver(mosaic, { createIndex: createPartitionAwareIndex });
let count = 0;
for (const [archive, part] of Object.entries(partitions)) {
  if (!(archive in mosaic.slices)) throw new Error(`Missing slice: ${archive}`);
  for (const xyz of Object.keys(part.tiles)) {
    const [z, x, y] = xyz.split(',').map(Number);
    if (!resolver.resolve(z, x, y).some(slice => slice.key === archive)) {
      throw new Error(`Mosaic cannot resolve ${xyz} to ${archive}`);
    }
    count++;
  }
}
console.log(JSON.stringify({ verified_tiles: count, archives: Object.keys(partitions).length,
  minzoom: mosaic.header.min_zoom, maxzoom: mosaic.header.max_zoom,
  mime_type: resolver.mimeType }));
