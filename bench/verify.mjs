#!/usr/bin/env node
// bench/verify.mjs — assertions on the exact mass-property core.
// Every analytic oracle here is an independent closed form (box, hollow box,
// tetrahedron, compound of two boxes); no constant is shared with src/core.mjs.
// Run: node bench/verify.mjs   (exit 0 = all pass)

import { readFileSync } from 'node:fs';
import {
  massProperties, combineMassProperties, transformVertices, rotationAxisAngle, symEigen3,
} from '../src/core.mjs';
import {
  vertexMeanCentroid, aabbMassProperties, pointMassInertia, voxelMassProperties,
} from '../baseline/naive.mjs';

const data = JSON.parse(readFileSync(new URL('../data/shapes.json', import.meta.url), 'utf8'));
const S = Object.fromEntries(data.shapes.map((s) => [s.name, s]));
const props = (name, opts) => massProperties(S[name].positions, S[name].indices, opts);

let passed = 0;
const failures = [];
function check(id, name, cond, detail) {
  if (cond) { passed++; console.log(`PASS ${id} ${name.padEnd(34)} ${detail}`); }
  else { failures.push(`${id} ${name}`); console.log(`FAIL ${id} ${name.padEnd(34)} ${detail}`); }
}
const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-300);
const matRel = (A, B) => {
  let d = 0, s = 0;
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
    d = Math.max(d, Math.abs(A[i][j] - B[i][j]));
    s = Math.max(s, Math.abs(B[i][j]));
  }
  return d / Math.max(s, 1e-300);
};
const matVec = (A, v) => A.map((row) => row[0] * v[0] + row[1] * v[1] + row[2] * v[2]);
const diag = (a, b, c) => [[a, 0, 0], [0, b, 0], [0, 0, c]];
const mul3 = (A, B) => [0, 1, 2].map((i) => [0, 1, 2].map((j) => A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j]));
const T3 = (A) => [0, 1, 2].map((i) => [0, 1, 2].map((j) => A[j][i]));

// ---------------------------------------------------------------- analytic
const BW = 4, BH = 2, BD = 1.5, MB = BW * BH * BD;
const boxVol = MB;
const boxI = diag(MB * (BH ** 2 + BD ** 2) / 12, MB * (BW ** 2 + BD ** 2) / 12, MB * (BW ** 2 + BH ** 2) / 12);

// -- A1..A3: the box is a polyhedron whose exact analytic moments are known.
const box = props('box');
check('A1', 'box.volume', rel(box.volume, boxVol) < 1e-14,
  `got ${box.volume.toFixed(15)} want ${boxVol} relErr ${rel(box.volume, boxVol).toExponential(2)}`);
check('A2', 'box.centroid', Math.hypot(...box.centroid) < 1e-13 * box.bboxDiag,
  `got [${box.centroid.map((v) => v.toExponential(2))}] |c| ${Math.hypot(...box.centroid).toExponential(2)}`);
check('A3', 'box.inertia(diag)', matRel(box.inertia, boxI) < 1e-13,
  `got [${[0, 1, 2].map((i) => box.inertia[i][i].toFixed(12))}] want [${[0, 1, 2].map((i) => boxI[i][i])}] rel ${matRel(box.inertia, boxI).toExponential(2)}`);
check('A3b', 'box.area', rel(box.area, 2 * (BW * BH + BH * BD + BW * BD)) < 1e-14,
  `got ${box.area.toFixed(12)} want ${2 * (BW * BH + BH * BD + BW * BD)}`);

// -- A4: same box rotated 0.7 rad about (1,2,3) and translated.
const rot = S.box_rot.params.rotate, tr = S.box_rot.params.translate;
const wantRotI = mul3(mul3(rot, boxI), T3(rot));
const boxRot = props('box_rot');
check('A4', 'box_rot.inertia = R I R^T', matRel(boxRot.inertia, wantRotI) < 1e-12,
  `rel ${matRel(boxRot.inertia, wantRotI).toExponential(2)} across 9 entries`);
check('A4b', 'box_rot.centroid = T', Math.hypot(boxRot.centroid[0] - tr[0], boxRot.centroid[1] - tr[1], boxRot.centroid[2] - tr[2]) < 1e-12,
  `got [${boxRot.centroid.map((v) => v.toFixed(9))}]`);

// -- A5: hollow box = outer solid minus inner solid (signed tetrahedra).
const HO = 3, HI = 1.6, mO = HO ** 3, mI = HI ** 3;
const wantHollowV = mO - mI;
const wantHollowI = diag(
  (mO * (HO ** 2 + HO ** 2) - mI * (HI ** 2 + HI ** 2)) / 12,
  (mO * (HO ** 2 + HO ** 2) - mI * (HI ** 2 + HI ** 2)) / 12,
  (mO * (HO ** 2 + HO ** 2) - mI * (HI ** 2 + HI ** 2)) / 12,
);
const hollow = props('hollow_box');
check('A5', 'hollow_box.volume/I', rel(hollow.volume, wantHollowV) < 1e-13 && matRel(hollow.inertia, wantHollowI) < 1e-13,
  `V ${hollow.volume} vs ${wantHollowV.toFixed(12)} (rel ${rel(hollow.volume, wantHollowV).toExponential(2)}), Ixx ${hollow.inertia[0][0].toFixed(9)} vs ${wantHollowI[0][0].toFixed(9)}`);

// -- A6: right tetrahedron, hand-derived: V=L^3/6, Ixx_c=L^5/80, Ixy_c=L^5/480.
const TL = 2.4;
const tetra = props('tetra');
const wantTetraI = [[TL ** 5 / 80, TL ** 5 / 480, TL ** 5 / 480],
  [TL ** 5 / 480, TL ** 5 / 80, TL ** 5 / 480],
  [TL ** 5 / 480, TL ** 5 / 480, TL ** 5 / 80]];
check('A6', 'tetra.volume/I (hand-derived)', rel(tetra.volume, TL ** 3 / 6) < 1e-14 && matRel(tetra.inertia, wantTetraI) < 1e-13,
  `V ${tetra.volume.toExponential(6)} vs ${(TL ** 3 / 6).toFixed(12)}; Ixx ${tetra.inertia[0][0].toFixed(9)} vs L^5/80=${(TL ** 5 / 80).toFixed(9)}`);

// -- A7: two disjoint boxes: parallel-axis analytic.
const m2 = 2 * 1 * 1.5;
const wantTwoI = diag(
  2 * (m2 * (1 ** 2 + 1.5 ** 2) / 12),
  2 * (m2 * (2 ** 2 + 1.5 ** 2) / 12 + m2 * 2.5 ** 2),
  2 * (m2 * (2 ** 2 + 1 ** 2) / 12 + m2 * 2.5 ** 2),
);
const two = props('two_box');
check('A7', 'two_box.volume/I (parallel axis)', rel(two.volume, 2 * m2) < 1e-14 && matRel(two.inertia, wantTwoI) < 1e-13,
  `V ${two.volume} Iyy ${two.inertia[1][1].toFixed(9)} vs ${wantTwoI[1][1].toFixed(9)} rel ${matRel(two.inertia, wantTwoI).toExponential(2)}`);

// -- A8: combineMassProperties(two halves) must reproduce the merged mesh.
const halfA = massProperties(S.two_box.positions.slice(0, 72), S.two_box.indices.slice(0, 36));
const halfB = massProperties(S.two_box.positions.slice(72), S.two_box.indices.slice(36).map((i) => i - 24));
const combined = combineMassProperties([halfA, halfB]);
check('A8', 'combine(two parts) == merged mesh', matRel(combined.inertia, two.inertia) < 1e-13 && rel(combined.mass, two.mass) < 1e-14,
  `mass ${combined.mass} vs ${two.mass}; Iyy ${combined.inertia[1][1].toFixed(9)} vs ${two.inertia[1][1].toFixed(9)}`);

// -- A9: translation invariance at 1e6 (catches origin-dependent integration).
const far = props('box', {}); // warm
const farProps = massProperties(transformVertices(S.box.positions, { translate: [1e6, -1e6, 3e6] }), S.box.indices);
check('A9', 'translation 1e6: V/I invariant', rel(farProps.volume, box.volume) < 1e-14 && matRel(farProps.inertia, box.inertia) < 1e-13,
  `V ${farProps.volume} (rel ${rel(farProps.volume, box.volume).toExponential(2)}), Ixx ${farProps.inertia[0][0].toFixed(12)} vs ${box.inertia[0][0].toFixed(12)}`);
void far;

// -- A10: scale invariance 1e-6 ... 1e6 (volume s^3, inertia s^5, verdict stable).
const scaleRows = [];
let scaleOk = true, verdictOk = true;
for (const s of [1e-6, 1e-3, 1, 1e3, 1e6]) {
  const p = massProperties(transformVertices(S.box.positions, { scale: s }), S.box.indices);
  const t = massProperties(transformVertices(S.box_tjunction.positions, { scale: s }), S.box_tjunction.indices);
  const ev = rel(p.volume, boxVol * s ** 3);
  const ei = matRel(p.inertia, boxI.map((r) => r.map((v) => v * s ** 5)));
  const ec = Math.hypot(p.centroid[0] / s, p.centroid[1] / s, p.centroid[2] / s);
  scaleRows.push(`s=${s.toExponential(0)} volRel=${ev.toExponential(1)} Irel=${ei.toExponential(1)}`);
  if (!(ev < 1e-13 && ei < 1e-12 && ec < 1e-12)) scaleOk = false;
  if (p.closed !== true || t.closed !== false) verdictOk = false;
}
check('A10', 'scale 1e-6..1e6 invariance', scaleOk && verdictOk, scaleRows.join(' | ') + ` verdictsStable=${verdictOk}`);

// -- A11: rotation invariance of trace / principal moments (2nd rotation).
const R2 = rotationAxisAngle([1, -4, 2], 1.1);
const rotBox2 = massProperties(transformVertices(S.box.positions, { rotate: R2 }), S.box.indices);
const tr2 = (M) => M[0][0] + M[1][1] + M[2][2];
const pmRel = Math.max(...[0, 1, 2].map((i) => rel(rotBox2.principalMoments[i], box.principalMoments[i])));
check('A11', 'rotation covariance (trace, moments)', rel(tr2(rotBox2.inertia), tr2(box.inertia)) < 1e-13 && pmRel < 1e-12,
  `trace ${tr2(rotBox2.inertia).toFixed(12)} vs ${tr2(box.inertia).toFixed(12)}; principal moments rel ${pmRel.toExponential(2)}`);

// -- A12: eigenpair validity: I = V diag(λ) V^T and V orthonormal.
const Vv = props('box').principalAxes;
const recon = mul3(mul3(Vv, diag(...box.principalMoments)), T3(Vv));
let ortho = 0;
for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) ortho = Math.max(ortho, Math.abs(mul3(T3(Vv), Vv)[i][j] - (i === j ? 1 : 0)));
check('A12', 'eigen: I = Q diag(λ) Q^T', matRel(recon, box.inertia) < 1e-12 && ortho < 1e-14,
  `reconRel ${matRel(recon, box.inertia).toExponential(2)} orthonormality ${ortho.toExponential(2)} λ=[${box.principalMoments.map((v) => v.toFixed(6))}]`);

// -- A13: open / cracked meshes must be refused, not silently measured.
const open = props('open_box'), cracked = props('box_tjunction');
const openIfIgnored = rel(open.volume, 8); // the true closed box is 2x2x2
check('A13', 'open + T-junction refused', open.ok === false && cracked.ok === false
  && open.boundaryEdges === 4 && cracked.boundaryEdges === 84,
  `open: ok=${open.ok} boundaryEdges=${open.boundaryEdges} reason=${open.reason} (its volume ${open.volume.toFixed(6)} would be ${(openIfIgnored * 100).toFixed(2)}% off the true 8); crack: ok=${cracked.ok} boundaryEdges=${cracked.boundaryEdges}`);

// -- A14: degenerate triangle audit (zero-area triangle inserted).
const degPos = S.box.positions.slice();
const degIdx = S.box.indices.slice(0, 33);
degIdx.push(0, 0, 0);
const deg = massProperties(degPos, degIdx);
check('A14', 'degenerate triangle flagged', deg.degenerateTriangles === 1 && deg.reason.includes('degenerate-triangles'),
  `degenerateTriangles=${deg.degenerateTriangles} reason=${deg.reason}`);

// -- A15: honest baseline comparison on the same meshes.
const cone = props('cone');
const vm = vertexMeanCentroid(S.cone.positions);
const centroidErr = Math.hypot(vm.centroid[0] - cone.centroid[0], vm.centroid[1] - cone.centroid[1], vm.centroid[2] - cone.centroid[2]);
const sphere = props('sphere_uv');
const pm = pointMassInertia(S.sphere_uv.positions, sphere.mass);
const pmErr = rel(pm.inertia[0][0], sphere.inertia[0][0]);
const aabbHollow = aabbMassProperties(S.hollow_box.positions);
const aabbHollowErr = rel(aabbHollow.volume, hollow.volume);
const aabbBox = aabbMassProperties(S.box.positions);
const aabbBoxErr = rel(aabbBox.volume, box.volume);
check('A15', 'vertex-mean centroid (cone)', centroidErr / cone.centroid[1] > 0.5,
  `true COM y ${cone.centroid[1].toFixed(6)}, vertex mean y ${vm.centroid[1].toFixed(6)} -> off by ${(centroidErr / cone.centroid[1] * 100).toFixed(2)}% (${(centroidErr / 4 * 100).toFixed(1)}% of height)`);
check('A15b', 'point-mass inertia (sphere)', pmErr > 0.3,
  `point-mass Ixx ${pm.inertia[0][0].toFixed(6)} vs exact ${sphere.inertia[0][0].toFixed(6)} -> ${(pmErr * 100).toFixed(2)}%`);
check('A15c', 'AABB volume, hollow part', aabbHollowErr > 0.1,
  `AABB ${aabbHollow.volume} vs exact ${hollow.volume} -> +${(aabbHollowErr * 100).toFixed(2)}% (a hollow part billed as solid)`);
check('A15d', 'AABB volume, solid box (control)', aabbBoxErr < 1e-15,
  `AABB ${aabbBox.volume} vs exact ${box.volume} -> ${aabbBoxErr.toExponential(2)} (baseline is not a strawman everywhere)`);

// -- A16: independent sampling oracle vs exact, on the sphere.
const exactV = sphere.volume;
const v16 = voxelMassProperties(S.sphere_uv.positions, S.sphere_uv.indices, 16);
const v24 = voxelMassProperties(S.sphere_uv.positions, S.sphere_uv.indices, 24);
const e16 = rel(v16.volume, exactV), e24 = rel(v24.volume, exactV);
check('A16', 'voxel oracle converges to exact', Math.abs(e16) < 0.08 && Math.abs(e24) < 0.06 && Math.abs(e24) < Math.abs(e16),
  `mesh=${exactV.toFixed(6)} voxel16=${v16.volume.toFixed(6)} (${(e16 * 100).toFixed(2)}%) voxel24=${v24.volume.toFixed(6)} (${(e24 * 100).toFixed(2)}%)`);

// -- A17: exact values vs smooth analytic references (tessellation error only).
const coneSmoothV = (Math.PI * 1.5 ** 2 * 4) / 3;
const coneErr = rel(cone.volume, coneSmoothV);
const sphereV = (4 / 3) * Math.PI * 1.5 ** 3;
const sphereErr = rel(sphere.volume, sphereV);
const sphereI = 0.4 * sphereV * 1.5 ** 2;
const sphereIErr = rel(sphere.inertia[0][0], sphereI);
check('A17', 'faceted vs smooth reference', Math.abs(coneErr) < 0.02 && Math.abs(sphereErr) < 0.03 && Math.abs(sphereIErr) < 0.04,
  `cone V ${(coneErr * 100).toFixed(3)}%; sphere V ${(sphereErr * 100).toFixed(3)}%, I ${(sphereIErr * 100).toFixed(3)}% (32x16 facets, not an integration error)`);

// -- A18: determinism.
const a = JSON.stringify(massProperties(S.sphere_uv.positions, S.sphere_uv.indices));
const b = JSON.stringify(massProperties(S.sphere_uv.positions, S.sphere_uv.indices));
check('A18', 'deterministic output', a === b, `two runs byte-identical (${a.length} bytes)`);

// -- A19: sum of parts over a 3-part assembly equals the exact one-shot result.
const three = combineMassProperties([box, { ...hollow, centroid: [10, 0, 0] }, { ...tetra, centroid: [0, 10, 0] }]);
const threeTotal = box.mass + hollow.mass + tetra.mass;
check('A19', '3-part assembly mass/centroid', rel(three.mass, threeTotal) < 1e-14
  && Math.abs(three.centroid[0] - (hollow.mass * 10) / threeTotal) < 1e-12,
  `mass ${three.mass.toFixed(9)} == ${threeTotal.toFixed(9)}; cx ${three.centroid[0].toFixed(9)}`);

console.log(`\n${passed} passed / ${failures.length} failed`);
if (failures.length) console.log('failed: ' + failures.join(', '));
process.exit(failures.length ? 1 : 0);
