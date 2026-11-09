// Plays an animated SPZ clip (SPZ_NIANTIC_temporal) in PlayCanvas.
//
// Each frame, spz.TemporalPlayer evaluates in WASM only the Gaussians the clip's time index says may
// be visible, and the visible ones are pushed into a PlayCanvas GSplatContainer sized for the most
// any frame shows (the path PlayCanvas's own animated-splat (.gsaf) player uses), so per-frame work
// scales with what's on screen rather than the whole clip. Containers
// draw SH degree 0 only, so view-dependent color is dropped. WebGPU sorts every frame; on WebGL the
// sort runs on the CPU and can lag the motion by a frame.
//
// Serve the repository root (e.g. `python3 -m http.server`) and open
// http://localhost:8000/tools/4dgs/viewer/?clip=<url of an .spz>, or drop a file on the page.

import * as pc from "playcanvas";
import { CameraControls } from "playcanvas/scripts/esm/camera-controls.mjs";
import createSpzModule from "../../../dist/spz.js";

const SH_C0 = 0.28209479177387814;
const TEXEL_FLOATS = 4;
// Gaussians fainter than this don't count towards the bounds used to frame the camera.
const BOUNDS_MIN_OPACITY = 0.05;
// Clip times sampled for the bounds, since the poster frame alone can include far-off Gaussians
// that are invisible then.
const BOUNDS_SAMPLES = 16;
const CAMERA_DISTANCE_FACTOR = 1.5;

const toHalf = pc.FloatPacking.float2Half;
const sigmoid = (v) => 1 / (1 + Math.exp(-v));

const canvas = document.getElementById("canvas");
const statusEl = document.getElementById("status");
const bar = document.getElementById("bar");
const playButton = document.getElementById("play");
const timeline = document.getElementById("timeline");
const timeLabel = document.getElementById("time");
const loopBox = document.getElementById("loop");

const spz = await createSpzModule();
const device = await pc.createGraphicsDevice(canvas, { deviceTypes: ["webgpu", "webgl2"], antialias: false });
const app = new pc.Application(canvas, {
  graphicsDevice: device,
  mouse: new pc.Mouse(canvas),
  touch: new pc.TouchDevice(canvas),
  keyboard: new pc.Keyboard(window),
});
app.setCanvasFillMode(pc.FILLMODE_FILL_WINDOW);
app.setCanvasResolution(pc.RESOLUTION_AUTO);
window.addEventListener("resize", () => app.resizeCanvas());
app.start();

const camera = new pc.Entity("camera");
camera.addComponent("camera", { clearColor: new pc.Color(0.07, 0.07, 0.07), fov: 60 });
camera.addComponent("script");
app.root.addChild(camera);

let clip = null;  // { player, entity, setFrame, duration }
let playing = true;
let playbackSeconds = 0;

function showStatus(text) {
  statusEl.textContent = text;
}

// Bounds of the visible Gaussians at least BOUNDS_MIN_OPACITY opaque, after player.evaluateVisible.
function visibleBounds(player, count, bounds) {
  const positions = player.visiblePositions();
  const alphas = player.visibleAlphas();
  for (let k = 0; k < count; k++) {
    if (sigmoid(alphas[k]) < BOUNDS_MIN_OPACITY) continue;
    for (let c = 0; c < 3; c++) {
      bounds.min[c] = Math.min(bounds.min[c], positions[k * 3 + c]);
      bounds.max[c] = Math.max(bounds.max[c], positions[k * 3 + c]);
    }
  }
  return bounds;
}

function emptyBounds() {
  return { min: [Infinity, Infinity, Infinity], max: [-Infinity, -Infinity, -Infinity] };
}

function finiteOr(bounds, fallback) {
  return Number.isFinite(bounds.min[0]) ? bounds : fallback;
}

// Bounds of the Gaussians that are visible at some point in the clip (just the poster frame for a
// static cloud).
function clipBounds(player) {
  const bounds = emptyBounds();
  const start = player.clipStart();
  const span = player.clipEnd() - start;
  const samples = player.isAnimated() ? BOUNDS_SAMPLES : 1;
  for (let s = 0; s < samples; s++) {
    visibleBounds(player, player.evaluateVisible(start + (span * s) / samples), bounds);
  }
  return finiteOr(bounds, { min: [-1, -1, -1], max: [1, 1, 1] });
}

// A PlayCanvas splat with room for the most Gaussians a frame shows (player.maxVisible()), not the
// whole clip, so per-frame uploads and sorting scale with what's visible. Each frame packs the
// player's visible Gaussians into the first slots; scales and colors, constant per Gaussian, are
// converted once and gathered through the visible indices.
function createSplat(player, bounds) {
  const slots = Math.max(1, player.maxVisible());
  const container = new pc.GSplatContainer(device, slots, pc.GSplatFormat.createDefaultFormat(device));
  const texels = container.textureDimensions.x * container.textureDimensions.y;
  container.aabb.setMinMax(new pc.Vec3(...bounds.min), new pc.Vec3(...bounds.max));

  const logScales = player.scales();
  const dcColors = player.colors();
  const n = player.numPoints();
  const scaleHalves = new Uint16Array(n * 3);
  const colorHalves = new Uint16Array(n * 3);
  for (let i = 0; i < n * 3; i++) {
    scaleHalves[i] = toHalf(Math.exp(logScales[i]));
    colorHalves[i] = toHalf(0.5 + dcColors[i] * SH_C0);
  }
  const center = new Float32Array(texels * TEXEL_FLOATS);
  const rotation = new Uint16Array(texels * TEXEL_FLOATS);
  const scale = new Uint16Array(texels * TEXEL_FLOATS);
  const color = new Uint16Array(texels * TEXEL_FLOATS);
  const centers = new Float32Array(slots * 3);

  const entity = new pc.Entity("clip");
  entity.addComponent("gsplat", { resource: container });
  app.root.addChild(entity);

  // Uploads the count Gaussians player.evaluateVisible just produced. Its arrays are views into
  // WASM memory, so they're read here, not kept.
  function setFrame(count) {
    const indices = player.visibleIndices();
    const positions = player.visiblePositions();
    const quats = player.visibleRotations();
    const alphas = player.visibleAlphas();
    for (let k = 0; k < count; k++) {
      const i3 = indices[k] * 3;
      const k3 = k * 3;
      const k4 = k * 4;
      center[k4] = centers[k3] = positions[k3];
      center[k4 + 1] = centers[k3 + 1] = positions[k3 + 1];
      center[k4 + 2] = centers[k3 + 2] = positions[k3 + 2];
      const x = quats[k4], y = quats[k4 + 1], z = quats[k4 + 2], w = quats[k4 + 3];
      const inv = 1 / (Math.hypot(x, y, z, w) || 1);
      // PlayCanvas containers store rotations as (w, x, y, z).
      rotation[k4] = toHalf(w * inv);
      rotation[k4 + 1] = toHalf(x * inv);
      rotation[k4 + 2] = toHalf(y * inv);
      rotation[k4 + 3] = toHalf(z * inv);
      for (let c = 0; c < 3; c++) {
        scale[k4 + c] = scaleHalves[i3 + c];
        color[k4 + c] = colorHalves[i3 + c];
      }
      color[k4 + 3] = toHalf(sigmoid(alphas[k]));
    }
    for (const [name, data] of [["dataCenter", center], ["dataRotation", rotation],
                                ["dataScale", scale], ["dataColor", color]]) {
      const texture = container.getTexture(name);
      texture.lock().set(data);
      texture.unlock();
    }
    container.centers.set(centers.subarray(0, count * 3));  // used by the CPU sort on WebGL
    container.update(count);
    entity.gsplat.workBufferUpdate = pc.WORKBUFFER_UPDATE_ONCE;
    app.renderNextFrame = true;
  }
  return { entity, setFrame };
}

function frameCamera(bounds) {
  const center = bounds.min.map((v, k) => (v + bounds.max[k]) / 2);
  const radius = Math.hypot(...bounds.max.map((v, k) => v - bounds.min[k])) / 2;
  camera.setLocalPosition(center[0], center[1], center[2] + radius * CAMERA_DISTANCE_FACTOR * 2);
  camera.script.destroy("cameraControls");
  camera.script.create(CameraControls, {
    properties: { focusPoint: new pc.Vec3(...center), enableFly: false },
  });
}

// Decodes the file once into a TemporalPlayer, shows its poster frame, and plays from there.
async function load(bytes, name) {
  if (clip) {
    clip.entity.destroy();
    clip.player.delete();
    clip = null;
  }
  const player = new spz.TemporalPlayer(bytes, spz.CoordinateSystem.RUB);
  const bounds = clipBounds(player);
  const splat = createSplat(player, bounds);
  frameCamera(bounds);
  const duration = player.clipEnd() - player.clipStart();
  clip = { ...splat, player, duration };
  document.getElementById("drop").hidden = true;
  const backend = device.isWebGPU ? "WebGPU" : "WebGL2 (sorting may lag motion)";
  const summary = `${name}\n${player.numPoints().toLocaleString()} Gaussians`;
  if (!player.isAnimated()) {
    splat.setFrame(player.evaluateVisible(0));
    showStatus(`${summary}, static\n${backend}`);
    return;
  }
  loopBox.checked = player.loop() !== spz.SpzTemporalLoop.ONCE;
  playbackSeconds = player.posterTime() - player.clipStart();
  bar.hidden = false;
  showStatus(`${summary}, ${duration.toFixed(2)} s clip, up to ` +
             `${player.maxVisible().toLocaleString()} visible at once\n${backend}`);
  showFrame();
}

function currentClipTime() {
  const { player, duration } = clip;
  if (loopBox.checked) return player.clipTime(playbackSeconds);
  return player.clipStart() + Math.min(Math.max(playbackSeconds, 0), duration * (1 - Number.EPSILON));
}

function showFrame() {
  const t = currentClipTime();
  const { player } = clip;
  clip.setFrame(player.evaluateVisible(t));
  const progress = (t - clip.player.clipStart()) / clip.duration;
  timeline.value = String(progress);
  timeLabel.textContent = `${(t - clip.player.clipStart()).toFixed(2)} s`;
}

app.on("update", (dt) => {
  if (!clip?.player.isAnimated() || !playing) return;
  playbackSeconds += dt;
  showFrame();
});

playButton.addEventListener("click", () => {
  playing = !playing;
  playButton.textContent = playing ? "Pause" : "Play";
});
timeline.addEventListener("input", () => {
  if (!clip?.player.isAnimated()) return;
  playing = false;
  playButton.textContent = "Play";
  playbackSeconds = Number(timeline.value) * clip.duration;
  showFrame();
});

window.addEventListener("dragover", (e) => e.preventDefault());
window.addEventListener("drop", async (e) => {
  e.preventDefault();
  const file = e.dataTransfer.files[0];
  if (!file) return;
  await load(new Uint8Array(await file.arrayBuffer()), file.name);
});

const url = new URLSearchParams(location.search).get("clip");
if (url) {
  showStatus(`loading ${url}…`);
  const response = await fetch(url);
  if (!response.ok) {
    showStatus(`couldn't load ${url}: HTTP ${response.status}`);
  } else {
    await load(new Uint8Array(await response.arrayBuffer()), url);
  }
}
