# meshmass-cert

Exact mass properties of closed triangle meshes — volume, area, centre of mass and the full
inertia tensor — with an analytic oracle, and a refusal path for meshes that would silently
return garbage. On a 4×2×1.5 box the twelve analytic values match to **3.55e-16** relative
error; on the same cone the vertex-mean baseline misses the centre of mass by **88.24 %**.

Upstream: none

## What it is

Given `positions` + `indices`, `src/core.mjs` integrates the signed tetrahedra `(0,a,b,c)`
(divergence theorem), sums their first and second moments, and shifts the result to the centre
of mass with the parallel-axis theorem. Signed contributions make it exact for closed meshes of
any genus, including hollow shells (outer surface + inverted inner surface).

It also refuses: every result carries `ok`, `closed`, `boundaryEdges`, `nonManifoldEdges` and
`degenerateTriangles`. An open box is reported `ok:false, boundaryEdges:4`; a T-junction crack
box is `ok:false, boundaryEdges:84`. Ignoring `ok` on the open box gives a volume **16.67 %**
off the true value — which is exactly the bug this flag exists to prevent.

Real users: physics/rigid-body setup (needs `I` about the COM), 3D printing and mass estimation,
asset QA on glTF/CAD exports.

## Run it

```bash
node tools/gen-sample.mjs     # deterministic sample meshes -> data/shapes.json
node bench/verify.mjs         # 24 assertions, exit 0
node bench/bench.mjs          # timings -> bench/result.json
python3 -m http.server 4173   # then open http://localhost:4173 (three@0.169.0 via importmap)
```

No dependency to install. The browser imports the *same* `src/core.mjs` and cross-checks every
live number against `bench/result.json` on screen ("browser == node").

## Evidence (all from `bench/result.json` and the verify run)

| case | oracle | result |
| --- | --- | --- |
| box 4×2×1.5 | analytic box, `I=m(w²+d²)/12` | volume 12 (rel 0), `I` diag `[6.25, 18.25, 20]` rel **3.55e-16** |
| hollow box 3³ − 1.6³ | analytic outer − inner | V 22.904 rel **7.76e-16**, `Ixx` 38.752373333 |
| right tetra L=2.4 | hand-derived `Ixx=L⁵/80`, `Ixy=L⁵/480` | V 2.304, `Ixx` 0.995328000 exactly |
| two boxes at x=±2.5 | parallel-axis analytic | `Iyy` 40.625 rel 5.25e-16 |
| box rotated 0.7 rad + translated | `I' = R I Rᵀ`, `c' = T` | rel **4.17e-13**, centroid exact |
| box at +1e6 offset | translation invariance | volume rel 0, `Ixx` 6.25 identical |
| box × 1e-6 … 1e6 | `V∝s³`, `I∝s⁵` | rel ≤ 2.8e-16, closed/refused verdicts unchanged |

Baselines (`baseline/naive.mjs`, same meshes, all runnable):

| baseline | mesh | error |
| --- | --- | --- |
| vertex-mean centroid | cone (32-gon, 34 verts) | COM y 0.117647 vs 1.0 → **88.24 %** (22.1 % of height) |
| vertices as point masses | UV sphere 32×16 | `Ixx` 22.988698 vs 12.400301 → **+85.39 %** |
| bounding box | hollow box | 27 vs 22.904 → **+17.88 %** (a hollow part billed solid) |
| bounding box | solid box | 0.00 — the baseline is not a strawman everywhere |
| voxel n=16 / 24 / 32 | sphere volume | +3.11 % / +0.49 % / +0.169 %, 215 ms at n=32 |

Speed (`bench/result.json`, 1 warmup + 7 reps): exact sphere (960 tris) **median 2.654 ms**,
cold start 3.197 ms; the voxel oracle costs **214.876 ms** at n=32 and is still 0.169 % off.
Relative std-dev is 19–137 % because the box is shared with other builds: p95 spikes (13.54 ms on
a 1.49 ms job) are scheduler preemption, not algorithm variance — output is byte-identical across
runs (assertion A18), only wall-clock jitters.

`data/shapes.json` sha256 `dd87819d76160602f4fed1261ff4285af7b3799b53abc131a6312d5f14d856e6`.

## 本轮创新

第一次把「网格质量属性」做成带拒绝语义的证书，而不是一个静默返回数字的工具函数：解析闭式
oracle（盒/空心盒/手推四面体/双盒平行轴）零共享常量地校验到 1e-16，同时用拓扑审计把
非水密与 T-junction 网格判为 `ok:false`（开口盒若被强行使用会错 16.67%）；并把
尺度 1e-6…1e6、平移 1e6 的不变量写成断言。与台账里的几何证书的区别：那些判「有没有缺陷」，
这个判「一个物理量算得对不对」，且给出可执行的拒绝判据。

## Known limits

- **Absolute thresholds: none.** All three tolerances are relative to the bbox diagonal:
  degenerate triangle `0.5|u×v| ≤ 1e-14·d²`, zero volume `|V| ≤ 1e-13·d³`, topology weld
  `1e-12·d`. At scale 1e-6 the weld cell is 1e-18, at 1e6 it is 1e-6 — verdicts are unchanged
  across 1e-6…1e6 (assertion A10). Below float64 subnormals (~1e-308) nothing is guaranteed.
- Signed tetra integration needs a **closed** mesh; open meshes are refused rather than corrected.
  Self-intersecting meshes are integrated as their signed volume — the number is well defined but
  is not "the mass of a solid".
- Density is uniform; non-uniform density needs per-part `combineMassProperties`.
- Faceted approximations of smooth bodies differ from the smooth closed form by tessellation only:
  cone 0.641 %, sphere 1.596 % volume / 2.540 % inertia at 32×16.
- Timing on this 4-core VM has 19–137 % relative std-dev (see above); medians are the stable read.

## Layout

```
src/core.mjs          pure ES module: massProperties, combineMassProperties, edgeAudit, transformVertices
main.js               three.js rendering + drag rotate + live COM markers (green exact / red vertex-mean)
index.html            importmap (three@0.169.0), single entry point
tools/gen-sample.mjs  9 deterministic shapes -> data/shapes.json
bench/verify.mjs      24 assertions (A1–A19), analytic oracles, no shared constants
bench/bench.mjs       1 warmup + 7 reps, cold start, -> bench/result.json
bench/result.json     every number quoted above
baseline/naive.mjs    vertex-mean / AABB / point-mass / voxel-ray-parity baselines
```

License: MIT (see `LICENSE`).
