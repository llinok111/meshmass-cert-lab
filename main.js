// main.js — rendering + interaction only. All numbers come from src/core.mjs,
// the same module bench/verify.mjs and bench/bench.mjs run in node.
import * as THREE from 'three';
import { massProperties } from './src/core.mjs';
import { vertexMeanCentroid } from './baseline/naive.mjs';

const canvas = document.getElementById('c');
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x0d1117);
const camera = new THREE.PerspectiveCamera(45, 1, 0.01, 20000);
const group = new THREE.Group();
scene.add(group);
scene.add(new THREE.HemisphereLight(0xbfd4ff, 0x0b0f14, 1.4));
const key = new THREE.DirectionalLight(0xffffff, 2.2);
key.position.set(4, 7, 6);
scene.add(key);

const meshMat = new THREE.MeshStandardMaterial({ color: 0x3d8bfd, metalness: 0.15, roughness: 0.55, flatShading: true, side: THREE.DoubleSide });
const comMat = new THREE.MeshBasicMaterial({ color: 0x3fb950 });
const naiveMat = new THREE.MeshBasicMaterial({ color: 0xf85149 });

const data = await (await fetch('./data/shapes.json')).json();
const ref = await (await fetch('./bench/result.json')).json();
const sphere = new THREE.SphereGeometry(1, 16, 12);

let current = null;
const list = document.getElementById('list');
const hud = document.getElementById('hud');
const matchEl = document.getElementById('match');
const verdict = document.getElementById('verdict');

function clearGroup() {
  for (const child of [...group.children]) {
    group.remove(child);
    if (child.geometry && child.geometry !== sphere) child.geometry.dispose();
  }
}

function select(shape) {
  current = shape;
  clearGroup();
  const geom = new THREE.BufferGeometry();
  geom.setAttribute('position', new THREE.Float32BufferAttribute(shape.positions, 3));
  geom.setIndex(shape.indices);
  geom.computeVertexNormals();
  group.add(new THREE.Mesh(geom, meshMat));
  group.add(new THREE.LineSegments(new THREE.WireframeGeometry(geom), new THREE.LineBasicMaterial({ color: 0x9fb6d4, transparent: true, opacity: 0.18 })));

  const p = massProperties(shape.positions, shape.indices);
  const vm = vertexMeanCentroid(shape.positions);
  const size = Math.max(p.bboxDiag * 0.012, 1e-9);
  const m1 = new THREE.Mesh(sphere, comMat); m1.position.set(...p.centroid); m1.scale.setScalar(size);
  const m2 = new THREE.Mesh(sphere, naiveMat); m2.position.set(...vm.centroid); m2.scale.setScalar(size);
  group.add(m1, m2);

  // frame
  const dist = p.bboxDiag * 1.6;
  camera.position.set(dist * 0.55, dist * 0.35, dist * 0.75);
  camera.near = dist / 500; camera.far = dist * 60; camera.updateProjectionMatrix();
  camera.lookAt(new THREE.Vector3(...p.centroid));

  const cell = (k, v) => `<tr><td class="dim">${k}</td><td>${v}</td></tr>`;
  const f = (x, d = 4) => (Math.abs(x) < 1e-12 ? '0' : x.toFixed(d));
  hud.innerHTML = [
    cell('triangles', p.triangles),
    cell('volume', f(p.volume, 6)),
    cell('area', f(p.area, 4)),
    cell('COM x,y,z', p.centroid.map((v) => f(v, 5)).join(', ')),
    cell('Ixx / Iyy / Izz', [0, 1, 2].map((i) => f(p.inertia[i][i], 5)).join(' / ')),
    cell('principal', p.principalMoments.map((v) => f(v, 5)).join(' / ')),
    cell('closed', p.closed ? '<span class="ok">yes</span>' : `<span class="bad">no (${p.boundaryEdges} boundary edges)</span>`),
    cell('certified', p.ok ? '<span class="ok">ok</span>' : `<span class="bad">${p.reason}</span>`),
    cell('vertex-mean err', `${(Math.hypot(vm.centroid[0] - p.centroid[0], vm.centroid[1] - p.centroid[1], vm.centroid[2] - p.centroid[2]) / Math.max(p.bboxDiag, 1e-30) * 100).toFixed(2)}% of bbox`),
  ].join('');

  const node = ref.shapes[shape.name];
  const dev = Math.max(
    Math.abs(p.volume - node.volume) / Math.max(Math.abs(node.volume), 1e-300),
    Math.abs(p.inertia[0][0] - node.inertia[0][0]) / Math.max(Math.abs(node.inertia[0][0]), 1e-300),
  );
  matchEl.innerHTML = dev < 1e-12
    ? `browser == node: <span class="ok">yes</span> (max rel Δ ${dev.toExponential(1)})`
    : `browser != node: <span class="bad">${dev.toExponential(2)}</span>`;
  verdict.innerHTML = p.ok
    ? `<span class="ok">closed mesh, exact moments</span><br /><span class="dim">same core as bench</span>`
    : `<span class="bad">refused: ${p.reason}</span><br /><span class="dim">numbers shown for diagnosis only</span>`;
}

for (const shape of data.shapes) {
  const b = document.createElement('button');
  b.innerHTML = `${shape.name} <span class="dim">${shape.triangles} tris</span>`;
  b.onclick = () => {
    for (const other of list.children) other.classList.remove('on');
    b.classList.add('on');
    select(shape);
  };
  list.appendChild(b);
}
list.firstChild.classList.add('on');
select(data.shapes[0]);

let drag = null;
canvas.addEventListener('pointerdown', (e) => { drag = { x: e.clientX, y: e.clientY }; canvas.setPointerCapture(e.pointerId); });
canvas.addEventListener('pointerup', () => { drag = null; });
canvas.addEventListener('pointermove', (e) => {
  if (!drag) return;
  group.rotation.y += (e.clientX - drag.x) * 0.008;
  group.rotation.x += (e.clientY - drag.y) * 0.008;
  drag = { x: e.clientX, y: e.clientY };
});

function resize() {
  const w = canvas.clientWidth, h = canvas.clientHeight;
  if (canvas.width !== w || canvas.height !== h) renderer.setSize(w, h, false);
  camera.aspect = w / Math.max(h, 1);
  camera.updateProjectionMatrix();
  renderer.render(scene, camera);
  requestAnimationFrame(resize);
}
resize();
