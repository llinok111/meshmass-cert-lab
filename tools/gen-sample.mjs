#!/usr/bin/env node
// tools/gen-sample.mjs — deterministic sample meshes for the mass-property
// certificate. Run: node tools/gen-sample.mjs  ->  data/shapes.json
// No randomness, no network, no date: two runs produce identical bytes.

import { mkdirSync, writeFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const here = dirname(fileURLToPath(import.meta.url));
const out = resolve(here, '..', 'data', 'shapes.json');

// Four corners of a box face, in CCW order seen from outside.
function box(w, h, d, center = [0, 0, 0]) {
  const hx = w / 2, hy = h / 2, hz = d / 2;
  const [cx, cy, cz] = center;
  const P = [
    [hx, -hy, -hz], [hx, hy, -hz], [hx, hy, hz], [hx, -hy, hz],       // +X
    [-hx, -hy, hz], [-hx, hy, hz], [-hx, hy, -hz], [-hx, -hy, -hz],   // -X
    [-hx, hy, hz], [hx, hy, hz], [hx, hy, -hz], [-hx, hy, -hz],       // +Y
    [-hx, -hy, -hz], [hx, -hy, -hz], [hx, -hy, hz], [-hx, -hy, hz],   // -Y
    [-hx, -hy, hz], [hx, -hy, hz], [hx, hy, hz], [-hx, hy, hz],       // +Z
    [-hx, hy, -hz], [hx, hy, -hz], [hx, -hy, -hz], [-hx, -hy, -hz],   // -Z
  ];
  const positions = [], indices = [];
  for (const p of P) positions.push(p[0] + cx, p[1] + cy, p[2] + cz);
  for (let f = 0; f < 6; f++) {
    const o = 4 * f;
    indices.push(o, o + 1, o + 2, o, o + 2, o + 3);
  }
  return { positions, indices };
}

function inverted(src) {
  const indices = src.indices.slice();
  for (let i = 0; i < indices.length; i += 3) {
    const t = indices[i + 1]; indices[i + 1] = indices[i + 2]; indices[i + 2] = t;
  }
  return { positions: src.positions.slice(), indices };
}

function merge(a, b) {
  const positions = a.positions.concat(b.positions);
  const off = a.positions.length / 3;
  return { positions, indices: a.indices.concat(b.indices.map((i) => i + off)) };
}

// Box whose +X face is tessellated 21x21 (non-uniform vertex density on
// purpose: this is what CAD/adaptive exporters actually produce).
function denseFaceBox(size, n) {
  const base = box(size, size, size);
  // Drop the two +X triangles (first face) and its four corners.
  const px = size / 2;
  const positions = [], indices = [];
  const keep = new Set();
  for (let i = 0; i < base.indices.length; i++) if (i >= 6) keep.add(base.indices[i]);
  const remap = new Map();
  for (const i of keep) {
    remap.set(i, positions.length / 3);
    positions.push(base.positions[3 * i], base.positions[3 * i + 1], base.positions[3 * i + 2]);
  }
  for (let i = 6; i < base.indices.length; i += 3) {
    indices.push(remap.get(base.indices[i]), remap.get(base.indices[i + 1]), remap.get(base.indices[i + 2]));
  }
  const base0 = positions.length / 3;
  for (let y = 0; y <= n; y++) for (let z = 0; z <= n; z++) {
    positions.push(px, -size / 2 + (size * y) / n, -size / 2 + (size * z) / n);
  }
  const at = (y, z) => base0 + y * (n + 1) + z;
  for (let y = 0; y < n; y++) for (let z = 0; z < n; z++) {
    indices.push(at(y, z), at(y + 1, z), at(y + 1, z + 1));
    indices.push(at(y, z), at(y + 1, z + 1), at(y, z + 1));
  }
  return { positions, indices };
}

function uvSphere(r, seg, ring) {
  const positions = [], indices = [];
  positions.push(0, r, 0);
  for (let i = 1; i < ring; i++) {
    const phi = (Math.PI * i) / ring;
    const y = r * Math.cos(phi), rad = r * Math.sin(phi);
    for (let j = 0; j < seg; j++) {
      const th = (2 * Math.PI * j) / seg;
      positions.push(rad * Math.sin(th), y, rad * Math.cos(th));
    }
  }
  positions.push(0, -r, 0);
  const top = 0, bottom = positions.length / 3 - 1;
  const row = (i, j) => 1 + (i - 1) * seg + ((j + seg) % seg);
  for (let j = 0; j < seg; j++) indices.push(top, row(1, j), row(1, j + 1));
  for (let i = 1; i < ring - 1; i++) for (let j = 0; j < seg; j++) {
    const a = row(i, j), b = row(i + 1, j), c = row(i + 1, j + 1), d = row(i, j + 1);
    indices.push(a, b, d);
    indices.push(b, c, d);
  }
  for (let j = 0; j < seg; j++) indices.push(bottom, row(ring - 1, j + 1), row(ring - 1, j));
  return { positions, indices };
}

// Closed polygonal cone: apex fan + base fan. Vertex density is extremely
// non-uniform (1 apex vertex vs n rim vertices), which is exactly what makes
// a vertex-mean centroid wrong.
function cone(R, H, n) {
  const positions = [0, H, 0, 0, 0, 0];
  for (let i = 0; i < n; i++) {
    const th = (2 * Math.PI * i) / n;
    positions.push(R * Math.cos(th), 0, R * Math.sin(th));
  }
  const ring = (i) => 2 + ((i + n) % n);
  const indices = [];
  for (let i = 0; i < n; i++) indices.push(0, ring(i + 1), ring(i)); // side, outward
  for (let i = 0; i < n; i++) indices.push(1, ring(i), ring(i + 1)); // base, outward -Y
  return { positions, indices };
}

function tetra(L) {
  const positions = [0, 0, 0, L, 0, 0, 0, L, 0, 0, 0, L];
  // outward-facing winding for the face opposite each vertex
  const indices = [0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];
  return { positions, indices };
}

const shapes = [];
const push = (name, params, mesh) => shapes.push({
  name, params,
  positions: Array.from(mesh.positions, (v) => Number(v.toFixed(12))),
  indices: mesh.indices,
  triangles: mesh.indices.length / 3,
});

const BW = 4, BH = 2, BD = 1.5;
push('box', { kind: 'box', size: [BW, BH, BD] }, box(BW, BH, BD));

// same box, rotated 0.7 rad about (1,2,3) and translated: rotation covariance
const rotBox = (() => {
  const R = (() => {
    const [x, y, z] = [1 / Math.sqrt(14), 2 / Math.sqrt(14), 3 / Math.sqrt(14)];
    const c = Math.cos(0.7), s = Math.sin(0.7), C = 1 - c;
    return [[c + x * x * C, x * y * C - z * s, x * z * C + y * s],
      [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
      [z * x * C - y * s, z * y * C + x * s, c + z * z * C]];
  })();
  const T = [10, -4, 2.5];
  const b = box(BW, BH, BD);
  const p = [];
  for (let i = 0; i < b.positions.length; i += 3) {
    const [x, y, z] = [b.positions[i], b.positions[i + 1], b.positions[i + 2]];
    p.push(R[0][0] * x + R[0][1] * y + R[0][2] * z + T[0],
      R[1][0] * x + R[1][1] * y + R[1][2] * z + T[1],
      R[2][0] * x + R[2][1] * y + R[2][2] * z + T[2]);
  }
  return { positions: p, indices: b.indices, R, T };
})();
push('box_rot', { kind: 'box', size: [BW, BH, BD], rotate: rotBox.R, translate: rotBox.T },
  { positions: rotBox.positions, indices: rotBox.indices });

const HO = 3, HI = 1.6;
push('hollow_box', { kind: 'hollow_box', outer: HO, inner: HI },
  merge(box(HO, HO, HO), inverted(box(HI, HI, HI))));

const TL = 2.4;
push('tetra', { kind: 'tetra', edge: TL }, tetra(TL));

const b1 = box(2, 1, 1.5, [-2.5, 0, 0]), b2 = box(2, 1, 1.5, [2.5, 0, 0]);
push('two_box', { kind: 'two_box', parts: [{ size: [2, 1, 1.5], center: [-2.5, 0, 0] }, { size: [2, 1, 1.5], center: [2.5, 0, 0] }] },
  merge(b1, b2));

// A +X face tessellated 20x20 against single-triangle neighbours: the classic
// T-junction crack. Volume still evaluates to 27 by coincidence, topology says
// 84 boundary edges -> must be refused.
push('box_tjunction', { kind: 'box', size: [3, 3, 3], faceGrid: 20 }, denseFaceBox(3, 20));

push('cone', { kind: 'cone', radius: 1.5, height: 4, segments: 32 }, cone(1.5, 4, 32));

// Box with the +Y face missing: must be refused, not silently measured.
const openBox = (() => {
  const b = box(2, 2, 2);
  const positions = [], indices = [];
  const dropped = new Set(b.indices.slice(12, 18));
  const remap = new Map();
  for (let i = 0; i < b.positions.length / 3; i++) {
    if (dropped.has(i)) continue;
    remap.set(i, positions.length / 3);
    positions.push(b.positions[3 * i], b.positions[3 * i + 1], b.positions[3 * i + 2]);
  }
  for (let f = 0; f < 6; f++) for (let k = 0; k < 6; k++) {
    const idx = 6 * f + k;
    if (idx >= 12 && idx < 18) continue;
    indices.push(remap.get(b.indices[idx]));
  }
  return { positions, indices };
})();
push('open_box', { kind: 'open_box', size: [2, 2, 2], missingFace: '+Y' }, openBox);

push('sphere_uv', { kind: 'sphere', radius: 1.5, segments: 32, rings: 16 }, uvSphere(1.5, 32, 16));

const payload = {
  meta: {
    generator: 'tools/gen-sample.mjs',
    deterministic: true,
    units: 'L (unitless); density 1 mass/L^3',
    shapeCount: shapes.length,
  },
  shapes,
};

mkdirSync(dirname(out), { recursive: true });
writeFileSync(out, JSON.stringify(payload));
console.log(`wrote ${out}: ${shapes.length} shapes, ${shapes.reduce((s, x) => s + x.triangles, 0)} triangles, ` +
  `${(JSON.stringify(payload).length / 1024).toFixed(1)} KiB`);
for (const s of shapes) console.log(`  ${s.name.padEnd(15)} tris=${String(s.triangles).padStart(5)} verts=${String(s.positions.length / 3).padStart(5)}`);
