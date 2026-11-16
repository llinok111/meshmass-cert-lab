#!/usr/bin/env node
// bench/bench.mjs — timings + the machine-readable evidence file.
// Run: node bench/bench.mjs  ->  bench/result.json
// 1 warmup + 7 timed reps per case; cold start = very first call after import.

import { readFileSync, writeFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { massProperties } from '../src/core.mjs';
import { vertexMeanCentroid, aabbMassProperties, pointMassInertia, voxelMassProperties } from '../baseline/naive.mjs';

const dataPath = new URL('../data/shapes.json', import.meta.url);
const raw = readFileSync(dataPath);
const data = JSON.parse(raw.toString('utf8'));
const S = Object.fromEntries(data.shapes.map((s) => [s.name, s]));
const REPORTS = 7;

const stats = (xs) => {
  const s = [...xs].sort((a, b) => a - b);
  const n = s.length;
  const mean = s.reduce((a, b) => a + b, 0) / n;
  const sd = Math.sqrt(s.reduce((a, b) => a + (b - mean) ** 2, 0) / n);
  const q = (p) => s[Math.min(n - 1, Math.max(0, Math.ceil(p * n) - 1))];
  return { reps: n, min: s[0], median: q(0.5), p95: q(0.95), max: s[n - 1], mean, stddev: sd, relStddev: mean ? sd / mean : 0 };
};
const time = (fn) => {
  const t0 = process.hrtime.bigint();
  fn();
  return Number(process.hrtime.bigint() - t0) / 1e6;
};
const runCase = (label, fn) => {
  const cold = time(fn);
  fn(); // warmup
  const reps = [];
  for (let i = 0; i < REPORTS; i++) reps.push(time(fn));
  return { label, coldStartMs: cold, ...stats(reps) };
};

// ---- cold start: the first core call in this process
const coldProps = massProperties(S.sphere_uv.positions, S.sphere_uv.indices);

const timings = [];
timings.push(runCase('massProperties(sphere_uv, 960 tris)', () => massProperties(S.sphere_uv.positions, S.sphere_uv.indices)));
timings.push(runCase('massProperties(box_tjunction, 810 tris)', () => massProperties(S.box_tjunction.positions, S.box_tjunction.indices)));
timings.push(runCase('massProperties(all 9 shapes, 1920 tris)', () => {
  for (const s of data.shapes) massProperties(s.positions, s.indices);
}));
timings.push(runCase('baseline.voxelMassProperties(sphere_uv, n=32)', () => voxelMassProperties(S.sphere_uv.positions, S.sphere_uv.indices, 32)));

// ---- per-shape evidence
const shapes = {};
for (const s of data.shapes) {
  const p = massProperties(s.positions, s.indices);
  shapes[s.name] = {
    triangles: p.triangles, vertices: p.vertices, closed: p.closed, ok: p.ok, reason: p.reason,
    boundaryEdges: p.boundaryEdges, degenerateTriangles: p.degenerateTriangles,
    volume: p.volume, area: p.area, centroid: p.centroid, mass: p.mass,
    inertia: p.inertia, principalMoments: p.principalMoments, bboxDiag: p.bboxDiag,
  };
}

// ---- baseline errors on the same meshes (honest negative controls included)
const sphere = shapes.sphere_uv, cone = shapes.cone, hollow = shapes.hollow_box, box = shapes.box;
const vmCone = vertexMeanCentroid(S.cone.positions);
const centroidErr = Math.hypot(vmCone.centroid[0] - cone.centroid[0], vmCone.centroid[1] - cone.centroid[1], vmCone.centroid[2] - cone.centroid[2]);
const pm = pointMassInertia(S.sphere_uv.positions, sphere.mass);
const aabbHollow = aabbMassProperties(S.hollow_box.positions);
const aabbBox = aabbMassProperties(S.box.positions);
const voxel = [16, 24, 32].map((n) => {
  const v = voxelMassProperties(S.sphere_uv.positions, S.sphere_uv.indices, n);
  return { n, volume: v.volume, cells: v.cells, relErr: Math.abs(v.volume - sphere.volume) / sphere.volume };
});

const result = {
  generatedAt: new Date().toISOString(),
  node: process.version,
  sampleSha256: createHash('sha256').update(raw).digest('hex'),
  sampleBytes: raw.length,
  shapes,
  timings,
  baselines: {
    coneVertexMeanCentroid: {
      trueCentroidY: cone.centroid[1], vertexMeanY: vmCone.centroid[1],
      absErr: centroidErr, relErrOfTrueY: centroidErr / cone.centroid[1],
      pctOfHeight: centroidErr / S.cone.params.height,
    },
    spherePointMassInertia: {
      exactIxx: sphere.inertia[0][0], pointMassIxx: pm.inertia[0][0],
      relErr: Math.abs(pm.inertia[0][0] - sphere.inertia[0][0]) / sphere.inertia[0][0],
    },
    hollowAabbVolume: {
      exactVolume: hollow.volume, aabbVolume: aabbHollow.volume,
      relErr: Math.abs(aabbHollow.volume - hollow.volume) / hollow.volume,
    },
    solidBoxAabbVolume: {
      exactVolume: box.volume, aabbVolume: aabbBox.volume,
      relErr: Math.abs(aabbBox.volume - box.volume) / box.volume,
    },
    voxelSphere: voxel,
  },
};
writeFileSync(new URL('./result.json', import.meta.url), JSON.stringify(result, null, 2) + '\n');
console.log('coldStart(ms) ' + coldProps.volume.toFixed(3) + ' volume checked');
for (const t of timings) {
  console.log(`${t.label.padEnd(52)} median ${t.median.toFixed(3)} ms  min ${t.min.toFixed(3)}  p95 ${t.p95.toFixed(3)}  relSD ${(t.relStddev * 100).toFixed(1)}%  cold ${t.coldStartMs.toFixed(3)}`);
}
console.log('wrote bench/result.json');
