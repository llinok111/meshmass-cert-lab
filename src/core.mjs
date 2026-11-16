// src/core.mjs — exact mass properties of closed triangle meshes.
//
// Method: signed tetrahedron decomposition with the origin as apex
// (divergence theorem). Each triangle (a,b,c) contributes a tetra (0,a,b,c)
// with signed volume det(a,b,c)/6, first moment V*(a+b+c)/4 and second
// moment (V/20)*(aa^T + bb^T + cc^T + ss^T), s = a+b+c.  Because the
// contributions are signed, the result is origin-independent for a closed
// mesh, and the parallel-axis shift to the centroid then gives the inertia
// tensor about the centre of mass.
//
// This module is pure: no THREE, no DOM, no I/O, no environment reads.
// Units: positions in L, density = 1 mass/L^3 by default, so
//   volume ~ L^3, centroid ~ L, inertia ~ L^5, principal moments ~ L^5.
//
// Tolerances (all relative, see README "Known limits"):
//   * degenerate triangle: 0.5*|u x v| <= areaEpsRel * bboxDiag^2
//     (default areaEpsRel = 1e-14).  A triangle is never dropped; it only
//     counts as degenerate for the audit.
//   * zero-volume: |volume| <= volEpsRel * bboxDiag^3 (default 1e-13).
//   * vertex weld for the topology audit: grid cell weldEpsRel * bboxDiag
//     (default 1e-12); the geometric integration never welds anything.
// No absolute epsilon is used anywhere, so verdicts are scale invariant.

export const DEFAULTS = Object.freeze({ areaEpsRel: 1e-14, volEpsRel: 1e-13, weldEpsRel: 1e-12 });

function bboxOf(positions) {
  let x0 = Infinity, y0 = Infinity, z0 = Infinity;
  let x1 = -Infinity, y1 = -Infinity, z1 = -Infinity;
  for (let i = 0; i < positions.length; i += 3) {
    const x = positions[i], y = positions[i + 1], z = positions[i + 2];
    if (x < x0) x0 = x; if (x > x1) x1 = x;
    if (y < y0) y0 = y; if (y > y1) y1 = y;
    if (z < z0) z0 = z; if (z > z1) z1 = z;
  }
  return { min: [x0, y0, z0], max: [x1, y1, z1], size: [x1 - x0, y1 - y0, z1 - z0] };
}

// Weld vertices by position so that unwelded exports (duplicated corners,
// the norm in CAD/glTF exports) still get a meaningful topology audit.
// The cell size is relative to the bounding-box diagonal, never absolute.
export function weldMap(positions, eps) {
  const map = new Map();
  const canon = new Int32Array(positions.length / 3);
  for (let i = 0; i < positions.length; i += 3) {
    const key = Math.round(positions[i] / eps) + ',' +
      Math.round(positions[i + 1] / eps) + ',' + Math.round(positions[i + 2] / eps);
    let id = map.get(key);
    if (id === undefined) { id = map.size; map.set(key, id); }
    canon[i / 3] = id;
  }
  return { canon, unique: map.size };
}

// Edge audit on welded ids. Purely combinatorial: no epsilon, so scale free.
export function edgeAudit(positions, indices, weldEps) {
  const { canon } = weldMap(positions, weldEps);
  const counts = new Map();
  const tri = (indices.length / 3) | 0;
  for (let t = 0; t < tri; t++) {
    const a = canon[indices[3 * t]], b = canon[indices[3 * t + 1]], c = canon[indices[3 * t + 2]];
    for (const [i, j] of [[a, b], [b, c], [c, a]]) {
      const lo = i < j ? i : j, hi = i < j ? j : i;
      const key = lo * 4294967296 + hi;
      counts.set(key, (counts.get(key) || 0) + 1);
    }
  }
  let boundary = 0, nonManifold = 0;
  const boundaryEdges = [];
  for (const [key, n] of counts) {
    if (n === 1) {
      boundary++;
      if (boundaryEdges.length < 32) {
        boundaryEdges.push([Math.floor(key / 4294967296), key % 4294967296]);
      }
    } else if (n > 2) nonManifold++;
  }
  return { edges: counts.size, boundary, nonManifold, boundaryEdges };
}

// Jacobi eigenvalue decomposition of a symmetric 3x3. Returns eigenvalues
// sorted descending plus orthonormal eigenvectors (columns of `vectors`).
export function symEigen3(A) {
  const a = [[A[0][0], A[0][1], A[0][2]], [A[1][0], A[1][1], A[1][2]], [A[2][0], A[2][1], A[2][2]]];
  let V = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  const mul = (X, Y) => {
    const Z = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
      let s = 0;
      for (let k = 0; k < 3; k++) s += X[i][k] * Y[k][j];
      Z[i][j] = s;
    }
    return Z;
  };
  for (let sweep = 0; sweep < 60; sweep++) {
    const off = Math.abs(a[0][1]) + Math.abs(a[0][2]) + Math.abs(a[1][2]);
    const diag = Math.abs(a[0][0]) + Math.abs(a[1][1]) + Math.abs(a[2][2]);
    if (off <= 1e-17 * (diag || 1)) break;
    for (const [p, q] of [[0, 1], [0, 2], [1, 2]]) {
      const apq = a[p][q];
      if (apq === 0) continue;
      const theta = (a[q][q] - a[p][p]) / (2 * apq);
      const t = (theta >= 0 ? 1 : -1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      const R = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
      R[p][p] = c; R[q][q] = c; R[p][q] = s; R[q][p] = -s;
      const Rt = [[R[0][0], R[1][0], R[2][0]], [R[0][1], R[1][1], R[2][1]], [R[0][2], R[1][2], R[2][2]]];
      const next = mul(mul(Rt, a), R);
      for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) a[i][j] = next[i][j];
      V = mul(V, R);
    }
  }
  const order = [0, 1, 2].sort((i, j) => a[j][j] - a[i][i]);
  return {
    values: order.map((i) => a[i][i]),
    vectors: [[V[0][order[0]], V[0][order[1]], V[0][order[2]]],
      [V[1][order[0]], V[1][order[1]], V[1][order[2]]],
      [V[2][order[0]], V[2][order[1]], V[2][order[2]]]],
    diagonalOfA: [a[0][0], a[1][1], a[2][2]],
  };
}

// Exact mass properties of one closed triangle mesh.
// Returns { ok, reason, volume, area, mass, centroid, inertia, principalMoments, ... }.
export function massProperties(positions, indices, opts = {}) {
  const density = opts.density === undefined ? 1 : opts.density;
  const areaEpsRel = opts.areaEpsRel === undefined ? DEFAULTS.areaEpsRel : opts.areaEpsRel;
  const weldEpsRel = opts.weldEpsRel === undefined ? DEFAULTS.weldEpsRel : opts.weldEpsRel;
  const volEpsRel = opts.volEpsRel === undefined ? DEFAULTS.volEpsRel : opts.volEpsRel;
  const bbox = bboxOf(positions);
  const diag = Math.hypot(bbox.size[0], bbox.size[1], bbox.size[2]);
  if (!(diag > 0)) throw new Error('massProperties: empty or zero-extent mesh');

  // Integrate in bbox-centred coordinates. The signed tetra decomposition is
  // translation invariant in exact arithmetic, but evaluating det(a,b,c) with
  // |a| ~ 1e6 and a result ~ 1e1 loses every digit. Recentring keeps the
  // relative error at ~1e-16 for any offset (asserted in bench/verify.mjs).
  const ox = 0.5 * (bbox.min[0] + bbox.max[0]);
  const oy = 0.5 * (bbox.min[1] + bbox.max[1]);
  const oz = 0.5 * (bbox.min[2] + bbox.max[2]);
  const P = new Float64Array(positions.length);
  for (let i = 0; i < positions.length; i += 3) {
    P[i] = positions[i] - ox; P[i + 1] = positions[i + 1] - oy; P[i + 2] = positions[i + 2] - oz;
  }

  const areaEps = areaEpsRel * diag * diag;
  const tri = (indices.length / 3) | 0;
  let volume = 0, area = 0, degenerate = 0;
  let mx = 0, my = 0, mz = 0;
  let cxx = 0, cyy = 0, czz = 0, cxy = 0, cxz = 0, cyz = 0;

  for (let t = 0; t < tri; t++) {
    const ia = 3 * indices[3 * t], ib = 3 * indices[3 * t + 1], ic = 3 * indices[3 * t + 2];
    const ax = P[ia], ay = P[ia + 1], az = P[ia + 2];
    const bx = P[ib], by = P[ib + 1], bz = P[ib + 2];
    const cx = P[ic], cy = P[ic + 1], cz = P[ic + 2];
    const ux = bx - ax, uy = by - ay, uz = bz - az;
    const vx = cx - ax, vy = cy - ay, vz = cz - az;
    const nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
    const a2 = Math.hypot(nx, ny, nz);
    if (0.5 * a2 <= areaEps) degenerate++;
    area += 0.5 * a2;

    const det = ax * nx + ay * ny + az * nz; // a . (u x v) = det(a,b,c)
    const V = det / 6;
    volume += V;
    mx += V * (ax + bx + cx) / 4;
    my += V * (ay + by + cy) / 4;
    mz += V * (az + bz + cz) / 4;
    // s = a+b+c (the tetra's far vertex from the origin-apex decomposition)
    const sx = ax + bx + cx, sy = ay + by + cy, sz = az + bz + cz;
    const k = V / 20;
    cxx += k * (ax * ax + bx * bx + cx * cx + sx * sx);
    cyy += k * (ay * ay + by * by + cy * cy + sy * sy);
    czz += k * (az * az + bz * bz + cz * cz + sz * sz);
    cxy += k * (ax * ay + bx * by + cx * cy + sx * sy);
    cxz += k * (ax * az + bx * bz + cx * cz + sx * sz);
    cyz += k * (ay * az + by * bz + cy * cz + sy * sz);
  }

  const audit = edgeAudit(positions, indices, weldEpsRel * diag);
  const closed = audit.boundary === 0 && audit.nonManifold === 0;
  const mass = density * volume;
  const cx = (volume !== 0 ? mx / volume : 0) + ox;
  const cy = (volume !== 0 ? my / volume : 0) + oy;
  const cz = (volume !== 0 ? mz / volume : 0) + oz;

  // Inertia about the origin, then shift to the centroid.
  const Io = [[cyy + czz, -cxy, -cxz], [-cxy, cxx + czz, -cyz], [-cxz, -cyz, cxx + cyy]];
  const d = density;
  const shift = (m, c) => [
    [m * (c[1] * c[1] + c[2] * c[2]), -m * c[0] * c[1], -m * c[0] * c[2]],
    [-m * c[0] * c[1], m * (c[0] * c[0] + c[2] * c[2]), -m * c[1] * c[2]],
    [-m * c[0] * c[2], -m * c[1] * c[2], m * (c[0] * c[0] + c[1] * c[1])],
  ];
  const S = shift(mass, [cx - ox, cy - oy, cz - oz]); // I_o is about the local origin
  const I = [0, 1, 2].map((i) => [0, 1, 2].map((j) => d * Io[i][j] - S[i][j]));
  for (let i = 0; i < 3; i++) for (let j = i + 1; j < 3; j++) {
    const sym = 0.5 * (I[i][j] + I[j][i]);
    I[i][j] = sym; I[j][i] = sym;
  }
  const eig = symEigen3(I);

  const reasons = [];
  if (!closed) reasons.push('mesh-not-closed');
  if (Math.abs(volume) <= volEpsRel * diag * diag * diag) reasons.push('zero-volume');
  if (degenerate > 0) reasons.push('degenerate-triangles');

  return {
    ok: reasons.length === 0,
    reason: reasons.length ? reasons.join('+') : null,
    closed,
    triangles: tri,
    vertices: positions.length / 3,
    degenerateTriangles: degenerate,
    areaEps,
    boundaryEdges: audit.boundary,
    nonManifoldEdges: audit.nonManifold,
    boundaryEdgeList: audit.boundaryEdges,
    edgeCount: audit.edges,
    volume, area, density, mass,
    centroid: [cx, cy, cz],
    inertia: I,
    principalMoments: eig.values,
    principalAxes: eig.vectors,
    bbox,
    bboxDiag: diag,
  };
}

// Sum of parts {mass, centroid, inertia} expressed in one frame
// (parallel-axis theorem). Real use: compound rigid bodies.
export function combineMassProperties(parts) {
  let m = 0;
  const c = [0, 0, 0];
  for (const p of parts) {
    m += p.mass;
    for (let i = 0; i < 3; i++) c[i] += p.mass * p.centroid[i];
  }
  if (m === 0) throw new Error('combineMassProperties: zero total mass');
  for (let i = 0; i < 3; i++) c[i] /= m;
  const I = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (const p of parts) {
    const d = [p.centroid[0] - c[0], p.centroid[1] - c[1], p.centroid[2] - c[2]];
    const dd = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
      I[i][j] += p.inertia[i][j] + p.mass * ((i === j ? dd : 0) - d[i] * d[j]);
    }
  }
  return { mass: m, volume: m, centroid: c, inertia: I, principalMoments: symEigen3(I).values };
}

// Affine transform of a vertex buffer. rotate: row-major 3x3 or null.
export function transformVertices(positions, { scale = 1, rotate = null, translate = null } = {}) {
  const out = new Float64Array(positions.length);
  const R = rotate || [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  const T = translate || [0, 0, 0];
  for (let i = 0; i < positions.length; i += 3) {
    const x = positions[i] * scale, y = positions[i + 1] * scale, z = positions[i + 2] * scale;
    out[i] = R[0][0] * x + R[0][1] * y + R[0][2] * z + T[0];
    out[i + 1] = R[1][0] * x + R[1][1] * y + R[1][2] * z + T[1];
    out[i + 2] = R[2][0] * x + R[2][1] * y + R[2][2] * z + T[2];
  }
  return out;
}

export function rotationAxisAngle(axis, angle) {
  const n = Math.hypot(axis[0], axis[1], axis[2]);
  const [x, y, z] = [axis[0] / n, axis[1] / n, axis[2] / n];
  const c = Math.cos(angle), s = Math.sin(angle), C = 1 - c;
  return [
    [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
    [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
    [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
  ];
}
