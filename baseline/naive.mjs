// baseline/naive.mjs — the three shortcuts real engines/exports take when
// they need "mass properties" without an exact integration. All are kept
// runnable so the numbers in README can be reproduced:
//   node --input-type=module -e "..."  or via bench/verify.mjs.
//
// 1. vertexMean: centroid = arithmetic mean of the vertex positions.
// 2. aabb:        the body is approximated by its bounding box.
// 3. pointMass:   every vertex is a point mass (m/N) at the vertex-mean centroid.
// 4. voxel:       uniform grid + ray-parity inside test (independent of core).

const mean = (positions) => {
  let x = 0, y = 0, z = 0;
  const n = positions.length / 3;
  for (let i = 0; i < positions.length; i += 3) { x += positions[i]; y += positions[i + 1]; z += positions[i + 2]; }
  return [x / n, y / n, z / n];
};

export function vertexMeanCentroid(positions) {
  const c = mean(positions);
  const off = Math.hypot(...c);
  return { centroid: c, offsetFromOrigin: off };
}

export function aabbMassProperties(positions) {
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < positions.length; i += 3) for (let k = 0; k < 3; k++) {
    lo[k] = Math.min(lo[k], positions[i + k]);
    hi[k] = Math.max(hi[k], positions[i + k]);
  }
  const s = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
  const m = s[0] * s[1] * s[2];
  const c = [(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2];
  const I = [
    [m * (s[1] ** 2 + s[2] ** 2) / 12, 0, 0],
    [0, m * (s[0] ** 2 + s[2] ** 2) / 12, 0],
    [0, 0, m * (s[0] ** 2 + s[1] ** 2) / 12],
  ];
  return { volume: m, centroid: c, mass: m, inertia: I };
}

// Point-mass inertia about the vertex-mean centroid. Pass the true mass to
// isolate the distribution error from the mass error; default is the AABB mass.
export function pointMassInertia(positions, mass) {
  const c = mean(positions);
  const aabb = aabbMassProperties(positions);
  const m = mass === undefined ? aabb.mass : mass;
  const n = positions.length / 3;
  const w = m / n;
  const I = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < positions.length; i += 3) {
    const d = [positions[i] - c[0], positions[i + 1] - c[1], positions[i + 2] - c[2]];
    const dd = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) I[a][b] += w * ((a === b ? dd : 0) - d[a] * d[b]);
  }
  return { centroid: c, mass: m, inertia: I };
}

// Ray parity along +X. Independent of src/core.mjs: no shared constants, no
// shared code path; only the triangle soup is shared.
export function rayParityInside(px, py, pz, positions, indices) {
  let hits = 0;
  for (let t = 0; t < indices.length; t += 3) {
    const i0 = 3 * indices[t], i1 = 3 * indices[t + 1], i2 = 3 * indices[t + 2];
    const ax = positions[i0], ay = positions[i0 + 1], az = positions[i0 + 2];
    const bx = positions[i1], by = positions[i1 + 1], bz = positions[i1 + 2];
    const cx = positions[i2], cy = positions[i2 + 1], cz = positions[i2 + 2];
    // barycentric in the YZ projection
    const e1y = by - ay, e1z = bz - az, e2y = cy - ay, e2z = cz - az;
    const wy = py - ay, wz = pz - az;
    const den = e1y * e2z - e1z * e2y;
    if (den === 0) continue;
    const u = (wy * e2z - wz * e2y) / den;
    const v = (e1y * wz - e1z * wy) / den;
    if (u < 0 || v < 0 || u + v > 1) continue;
    const x = ax + u * (bx - ax) + v * (cx - ax);
    if (x > px) hits++;
  }
  return (hits & 1) === 1;
}

// Uniform grid sampling: volume/centroid/inertia from inside cell centres.
export function voxelMassProperties(positions, indices, n) {
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < positions.length; i += 3) for (let k = 0; k < 3; k++) {
    lo[k] = Math.min(lo[k], positions[i + k]);
    hi[k] = Math.max(hi[k], positions[i + k]);
  }
  const s = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
  const cellVol = (s[0] / n) * (s[1] / n) * (s[2] / n);
  let inside = 0;
  const acc = [0, 0, 0];
  const centers = [];
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) for (let k = 0; k < n; k++) {
    const p = [lo[0] + (i + 0.5) * s[0] / n, lo[1] + (j + 0.5) * s[1] / n, lo[2] + (k + 0.5) * s[2] / n];
    if (!rayParityInside(p[0], p[1], p[2], positions, indices)) continue;
    inside++;
    centers.push(p);
    for (let d = 0; d < 3; d++) acc[d] += p[d];
  }
  const volume = inside * cellVol;
  const c = inside ? acc.map((v) => v / inside) : [0, 0, 0];
  const m = volume / Math.max(inside, 1);
  const I = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (const p of centers) {
    const d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
    const dd = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) I[a][b] += m * ((a === b ? dd : 0) - d[a] * d[b]);
  }
  return { volume, centroid: c, mass: volume, inertia: I, cells: n * n * n, insideCells: inside };
}
