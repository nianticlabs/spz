# 4DGS tools

Scripts for making, converting and inspecting animated SPZ clips (`SPZ_NIANTIC_temporal`, see
[extensions/temporal-design.md](../../extensions/temporal-design.md)).

They need the `spz` Python package built with extension support:

```sh
pip install . -Ccmake.define.SPZ_BUILD_EXTENSIONS=ON   # from the repository root
```

| Script | What it does | Needs |
| --- | --- | --- |
| `make_synthetic_clip.py` | Writes two small test clips: a polynomial/window clip and a skinned clip | `spz` |
| `spz_to_ply_frames.py` | Exports a clip as a numbered PLY sequence (plays in SuperSplat) | `spz` |
| `temporal_ply_to_spz.py` | Converts a PLY in the temporal PLY convention (STG/FreeTimeGS-style) | `spz` |
| `ply_sequence_to_spz.py` | Converts per-frame 3DGS PLYs (untracked) into one clip, frame-exact via visibility intervals | `spz` |
| `sog_sequence_to_spz.py` | The same for PlayCanvas SOG frames (`.sog` files or folders), with `--first/--count/--step/--sh-degree` | `spz`, Pillow (WebP) |
| `sog.py` | Reads SOG v2 splats (PlayCanvas splat-transform's format) | Pillow (WebP) |
| `fudan4dgs_to_spz.py` | Converts a [Fudan 4DGS](https://github.com/fudan-zvg/4d-gaussian-splatting) checkpoint (exact geometry) | `spz`, `torch` (CPU is fine) |
| `scgs_export.py` | Exports a trained [SC-GS](https://github.com/yihua7/SC-GS) model to `.npz`; run inside the SC-GS repo | SC-GS environment (CUDA) |
| `scgs_to_spz.py` | Converts that export into a skinned clip and reports the error against SC-GS | `spz` |
| `spz4d.py` | Shared helpers, including object ids (`set_object_ids`, `rigid_objects_clip`) | |
| `conformance.py`, `check_conformance.mjs` | Regenerate the conformance corpus; check the WASM reader against it | `spz`; Node |
| `spz_layout.py` | Independent parser of the file layout, for tests and the corpus | |

The CLI tools also understand clips: `spz_info clip.spz` prints the clip, its channels with their
encodings and declared maximum error, and its drivers; `spz_to_ply clip.spz frame.ply --time 1.5`
exports one frame.

**Precision.** The converters store motion as fixed-step channels with declared error bounds rather
than bit depths:

| Flag | Bounds | Default |
| --- | --- | --- |
| `--position-error` | each position polynomial term | 0.5 mm |
| `--rotation-error` | each rotation term component | 0.001 |
| `--time-error` | window centers and intervals | 50 µs |
| `--log-sigma-error` | ln(window sigma) | 0.001 |
| `--opacity-error` | window peak | 1/512 |

Terms are displacements across each Gaussian's visible span, so these bounds hold for short- and
long-lived Gaussians alike. `scgs_to_spz.py` takes `--weight-error` for skin weights (indices are
exact).

**Checking readers.** `python -m pytest tests/python/test_conformance.py` checks the C++ reader,
and `node tools/4dgs/check_conformance.mjs dist/spz.js tests/conformance` the WASM build.

## Testing on D-NeRF

[D-NeRF](https://github.com/albertpumarola/D-NeRF) is a set of images with camera poses and
timestamps, not trained splats, so a 4DGS model has to be trained first, on a machine with a CUDA
GPU. Download the dataset from the link in the D-NeRF README and unpack it so each scene is at
`data/dnerf/<scene>` (e.g. `jumpingjacks`, `mutant`, `standup`, `trex`, `hellwarrior`, `hook`,
`lego`, `bouncingballs`).

D-NeRF's model time runs over [0, 1]; `--duration` sets the clip length in seconds. Its scenes are
Blender worlds (x right, y forward, z up), so the converters default to `--from-coord RFU`, which
makes them upright in Y-up viewers.

### Per-Gaussian motion: Fudan 4DGS

Each 4D Gaussian sliced at time t is a 3D Gaussian moving linearly with a Gaussian visibility
window, which the clip stores exactly. Time-varying color is reduced to plain SH at each Gaussian's
peak time.

```sh
git clone https://github.com/fudan-zvg/4d-gaussian-splatting && cd 4d-gaussian-splatting
# set up the environment as in its README; put the data in data/dnerf/<scene>
python train.py --config configs/dnerf/jumpingjacks.yaml
# -> output/dnerf/jumpingjacks/chkpnt30000.pth (bouncingballs: chkpnt20000.pth)

python /path/to/spz/tools/4dgs/fudan4dgs_to_spz.py \
    output/dnerf/jumpingjacks/chkpnt30000.pth jumpingjacks_4dgs.spz --duration 5
```

`chkpnt_best.pth` (best test PSNR) works too.

### Skinned motion: SC-GS

SC-GS moves Gaussians with a few hundred control nodes, which become the clip's rigid drivers: a
Gaussian stores only its node indices and weights. Positions are exact. Orientations are an
approximation: SC-GS rotates Gaussians with a separately learned quaternion residual, while SPZ
skinning rotates them with their nodes (as glTF and USD skins do); `scgs_to_spz.py` prints the
position, rotation and scale error against SC-GS's own output, so you can judge it per scene.
Train with `--local_frame`, or nodes only translate.

```sh
git clone https://github.com/yihua7/SC-GS && cd SC-GS
# set up the environment as in its README
python train_gui.py --source_path data/dnerf/jumpingjacks --model_path outputs/jumpingjacks \
    --deform_type node --node_num 512 --hyper_dim 8 --is_blender --eval \
    --gt_alpha_mask_as_scene_mask --local_frame --resolution 2 --W 800 --H 800
# -> outputs/jumpingjacks_node/

python /path/to/spz/tools/4dgs/scgs_export.py -m outputs/jumpingjacks --out jumpingjacks_scgs.npz
# then, anywhere with the spz package:
python /path/to/spz/tools/4dgs/scgs_to_spz.py jumpingjacks_scgs.npz jumpingjacks_skinned.spz --duration 5
```

### Looking at the result

- **Real-time playback in PlayCanvas** (`viewer/`): build the WASM module with extensions (it writes
  `dist/spz.js`), serve the repository root, and open the viewer with a clip:

  ```sh
  emcmake cmake -S . -B build-wasm -DSPZ_BUILD_EXTENSIONS=ON && cmake --build build-wasm
  python3 -m http.server 8000
  # http://localhost:8000/tools/4dgs/viewer/?clip=/path/under/repo/clip.spz  (or drop a file on the page)
  ```

  Frames are evaluated in WASM (`spz.TemporalPlayer`) and pushed into a PlayCanvas
  `GSplatContainer`. Containers draw base color only (no higher-order SH). WebGPU re-sorts every
  frame; on WebGL2 the CPU sort can lag the motion by a frame.

- `spz_info jumpingjacks_4dgs.spz` checks that the clip binds.
- `python spz_to_ply_frames.py jumpingjacks_4dgs.spz frames/` writes a PLY per frame; drag the
  folder into the [SuperSplat editor](https://superspl.at/editor) and use its timeline.
- Any SPZ viewer shows the poster frame (the middle of the clip for 4DGS, the bind pose for SC-GS).

### Frame sequences are big

Untracked sequences (independently reconstructed frames, like most per-frame exports) can't share
Gaussians between frames, so a clip holds every frame's Gaussians: N frames of M Gaussians is N·M
Gaussians, each visible only in its frame. Playback stays fast: the clip's time index gives one bin
per frame, so the viewer evaluates and uploads only the current frame's M Gaussians (about 5 ms per
frame for M = 100k). But the file and the memory holding the clip scale with N·M. For example, 164 frames of 100k Gaussians at SH
degree 0 make a 243 MB clip and needed about 10 GB of memory to convert. Keep such clips to a few
seconds, or convert excerpts. Tracked or trained 4DGS models (Fudan 4DGS, SC-GS) are far smaller.

