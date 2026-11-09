# Design: dynamic Gaussians (4DGS) in SPZ — `SPZ_NIANTIC_temporal`

Status: implemented on `lreyna/4dgs-extensino` (phases 0–4, plus skinning and the JS player); the
PlayCanvas viewer is in `tools/4dgs/viewer`. The normative spec is in
[README.md](README.md#built-in-extensions); this document explains the design.

## 1. Decision

An SPZ file can hold an **animated clip**: one set of Gaussians whose state is an explicit function
of time (the "parametric per-Gaussian" family, *A*). Each Gaussian can carry:

- a polynomial position curve;
- a polynomial rotation curve;
- a Gaussian-shaped opacity window that says when it is visible.

One self-contained record, `SPZ_NIANTIC_temporal`, describes the **clip** (start, end, poster time)
and its **channels**. Each channel animates one target property (position, rotation, opacity) with
one channel type, and holds its own per-point data in the record's extension frames, encoded like
custom attributes. This mirrors glTF animation, where a channel targets a property and a sampler
says how it varies.

An earlier draft kept channel data in separately named custom attributes. Binding by name needed a
dependency mechanism in the decoder, validation of dangling or mistyped references on load and
save, and an attribute semantic used only by this record. Owning the data removes all of that: a
clip is one record, valid or not as a whole.

v1 also supports **grouped motion through skinning** (§12a): a few hundred rigid drivers, keyframed
in the record, and per-Gaussian driver indices and weights (`TRANSFORM` / `SKIN`), as in SC-GS and
skeletal rigs. Keyframed Gaussians (family *B*) are not in v1. They arrive later as a new channel
type, `KEYFRAMES`, in the same record (see [§12](#12-future-keyframed-gaussians-family-b)).

Why A first:

- **Bakes down losslessly.** A can be sampled into keyframes or per-frame clouds (MPEG GSC's
  I-3DGS sequences) without loss; the reverse is a lossy fit. Storing A keeps every later target
  open.
- **Covers the field.** Fudan 4DGS (exactly, see §9), FreeTimeGS, `.4dgs` gaussian-birth, and
  Spacetime Gaussians' (STG) geometry all reduce to it. v1 implements the most general form
  (cubic position, rotation curve, opacity window); simpler methods use a subset.
- **Cheapest to render.** A few multiply-adds per Gaussian in a shader; no extra fetches.
- **Compact.** Size scales with the Gaussian count, not the frame count (§10).

## 2. Goals and non-goals

Goals:

- Old readers (any SPZ v4 reader, builds without extensions) render a correct still frame.
- Exact, normative evaluation, so every renderer draws the same frame for the same time.
- No core format change: reuse extension frames and the custom-attribute value encoding.
- v1 is **single-file**: one SPZ holds one complete clip, typically up to 5–10 s.
- One generic import path, a documented temporal PLY convention (§9), that any training pipeline
  can target.

Non-goals for v1:

- Long content, segments, and manifests ([§8](#8-future-long-content-segments-and-manifest)). The
  clip fields already work for segments, so adding them later changes nothing in the file.
- Time-varying color / SH (needs an SH-rotation attribute semantic; research has not converged).
- Time-varying scale.
- Neural deformation fields (MLPs, hexplanes): bake them to curves on export.
- Audio, cameras, object tracks.

## 3. Evaluation (normative)

Let `t` be time in seconds. For each Gaussian with core position `p`, rotation `q`, alpha `α₀`, and
the record's channels, with `u = (t − o) / s` its normalized time (below):

```
position   p_l = p + Σ_{k=1..d} d_k · u^k                    (POSITION / POLYNOMIAL, else p)
rotation   q_l = normalize(q + Σ_{k=1..e} e_k · u^k)         (ROTATION / POLYNOMIAL, else q)
skin       p(t) = Σ_j ŵ_j (R_j(t) p_l + T_j(t))              (TRANSFORM / SKIN, else p_l)
           q(t) = normalize(Σ_j ŵ_j s_j r_j(t)) ⊗ q_l         (else q_l)
opacity    α(t) = peak · exp(−½ · ((t − c) / σ)²)            (OPACITY / WINDOW, else α₀)
           σ = exp(lnSigma); 0 outside [lo, hi) for windows with an interval
all        α(t) = 0 when t ∉ [clipStart, clipEnd)
scale, SH  unchanged
```

Per-point channels apply first, in the drivers' bind-pose frame; skinning then moves the result.
ŵ are the weights divided by their sum (points whose weights sum to 0 aren't skinned), and the
signs s_j = ±1 align each driver quaternion with the first one, so q and −q don't cancel.

A channel with a point range (§5) applies only to points in its range. Points outside it behave
as if the channel were absent.

A renderer MAY skip Gaussians with `α(t) < opacityCutoff` (record field) before sorting. It MUST
NOT skip others.

### Time frames: where polynomials are measured from

Each point's polynomials use its own origin `o` and span `s`, derived from data the file already
holds:

- The point's **visible span** `[lo, hi]` is the clip, intersected with `[c − Kσ, c + Kσ]` (where
  its window stays at or above the cutoff, `K = sqrt(−2 ln opacityCutoff)`) and with its interval.
- `o = clamp(posterTime, lo, hi)`, `s = hi − lo` (an empty span gives `posterTime` and the clip's
  duration).

Why not measure every polynomial from `posterTime`, as an earlier draft did? Re-centering is exact
algebra, but a stored coefficient's quantization error ε then shows up as `ε·|t − posterTime|^k`. A
short-lived Gaussian seen 2 s from the poster time would amplify its error 2^k-fold. Its core
position would also be a far extrapolation, invisible (α₀ ≈ 0) but inflating the poster frame's
bounds. With per-point frames:

- **Error is bounded by construction.** `|u| ≤ 1` wherever the point is visible, so its position
  error is at most the sum of its terms' errors plus the core's, for short- and long-lived Gaussians
  alike. This is the per-Gaussian precision `.4dgs` derives from lifetimes (§12b), here exact.
- **The poster frame is unchanged.** A point visible at `posterTime` has `o = posterTime`.
- **Bounds stay honest.** An invisible point is stored where it's next seen.

The frame needs the window's stored (quantized) values, so writers derive it from what they store;
the builder quantizes windows first.

### Core attributes are the poster frame

The core position, rotation, and alpha hold the state at each point's origin, which is `posterTime`
for every point visible then, and `α₀ = α(posterTime)`. Consequences:

- Readers that ignore the record render a clean frame at `posterTime`, not every instant
  superimposed. Gaussians invisible then get `α₀ ≈ 0` and sit on their visible path.
- 4D readers ignore `α₀` and use `peak`. This is why `peak` lives in the opacity channel: it can't
  be recovered from a near-zero `α₀`.
- Core data streams before extension frames, so a viewer can show the poster frame as soon as the
  core arrives and start animating once the record's data lands (`tools/4dgs/viewer` does this).

### Quaternion sign

`−q` is the same rotation as `q`, and the packer canonicalizes it: v3+ negates `q` so its largest
component is positive (`packQuaternionSmallestThree`, `src/cc/load-spz.cc`). `q + Σ e_k u^k` is
only correct if `e_k` share the stored sign of `q`. Therefore:

- `e_k` are defined relative to the core rotation **as stored** (after canonicalization).
- Writers MUST negate every `e_k` of a Gaussian whenever they negate its `q`. The library does this
  at pack time when a rotation channel is bound (§7). Writers outside the library must canonicalize
  `q` the way the packer does before deriving `b_k`.

## 4. Record format

`SPZ_NIANTIC_temporal`, type `0x4E530003`. A framed record: descriptors in the plaintext zone, data
in extension frames. It is keyed: at most one per file, and duplicates are dropped on load like other
keyed records.

```
f64 clipStart          // seconds; normally 0 (non-zero for future segments, §8)
f64 clipEnd            // seconds, exclusive; clipStart < clipEnd. Duration = clipEnd - clipStart
f64 posterTime         // time the core attributes show; clipStart <= posterTime < clipEnd
u8  loopMode           // playback hint: 0 ONCE, 1 LOOP, 2 PING_PONG
u8  reserved           // 0
f32 nominalFps         // playback hint only, 0 = unspecified; evaluation is continuous
f32 opacityCutoff      // culling threshold for α(t), in [0, 1); writers default to 1/255
u8  channelCount
channelCount x {
  u32 descriptorBytes  // length of the rest, so readers skip parameters they don't know
  u8  target           // 1 POSITION, 2 ROTATION, 3 OPACITY, 6 TRANSFORM; 4 SCALE, 5 COLOR reserved
  u8  type             // 1 POLYNOMIAL, 2 WINDOW, 4 SKIN; 3 KEYFRAMES reserved (§12)
  u8  flags            // bit 0 REQUIRED; others reserved, written as 0
  u16 dim              // floats per covered point
  u32 pointBegin, u32 pointEnd          // pointEnd 0: every point
  u64 dataBytes        // encoded size, so readers skip data they can't decode
  u8  encoding         // 0 FLOAT32; 1 RANGE: u8 bits, f32 min[dim], f32 max[dim]; 2 STEP: f32 step[dim]
}
u32 driverCount, u32 keyframeCount      // both 0, or both > 0 (§12a)
f64 keyframeTimes[keyframeCount]        // strictly increasing, shared by every driver
u32 nameCount, names                    // 0 or driverCount names (§12a)
frame table                             // driver transforms, then each channel's values
u32 binCount, f64 boundaries[binCount + 1], binCount x {u32 pointBegin, pointEnd}
                                        // time index: points that may be visible per bin (§11)
```

Two explicit lengths make the record extensible in both directions: `descriptorBytes` lets a
reader step over encoding parameters it doesn't know, and `dataBytes` over the data. A reader can
therefore skip an optional channel in a future encoding and still decode the rest; a required one
makes the file static. Readers MUST ignore bytes after the frame table (later versions may append
fields, e.g. the per-Gaussian keyframe time table in §12).

### Validation and fallback

The record is all-or-nothing. A reader **treats the file as static** (renders the poster frame,
logs why) when any of these hold:

- `clipStart`, `clipEnd`, or `posterTime` is non-finite, `clipStart >= clipEnd`, or
  `posterTime ∉ [clipStart, clipEnd)`.
- Two channels have the same target.
- A channel with `REQUIRED` set has a target or type the reader doesn't know, or a combination it
  doesn't support (v1: POSITION/POLYNOMIAL, ROTATION/POLYNOMIAL, OPACITY/WINDOW).
- A supported channel's `dim` doesn't fit its layout (§5), or its encoding parameters, point range
  or `dataBytes` are invalid, or a required channel's encoding is unknown.
- A skin index is outside the record's drivers, or a weight is negative.

Channels without `REQUIRED` whose target or type is unknown are kept (re-saving preserves them) and
otherwise ignored. v1 writers set `REQUIRED`
on every channel; the flag exists for future hint-only channels (e.g. velocity for motion blur).

## 5. Channels

| Target / type | `dim` | Builder encoding (default bound) |
|---|---|---|
| POSITION / POLYNOMIAL | `3·d`, `d` ∈ 1..3 | STEP, 0.5 mm per term |
| ROTATION / POLYNOMIAL | `4·e`, `e` ∈ 1..3 | STEP, 0.001 per term component |
| OPACITY / WINDOW | 3, or 5 with an interval | STEP: times 50 µs, ln σ 0.001, peak 1/512 |
| TRANSFORM / SKIN | `2·J`, `J` ∈ 1..8 | STEP: indices exact (step 1), weights 1/4096 |
| any / KEYFRAMES | reserved for B (§12) | |

Component layouts, per point:

- POSITION / POLYNOMIAL: `[d1x d1y d1z  d2x d2y d2z  d3x d3y d3z]`, displacements per power of `u`.
  The degree is `dim / 3`.
- ROTATION / POLYNOMIAL: `[e1x e1y e1z e1w  e2…]`, quaternion component order matching the core
  (x, y, z, w). The degree is `dim / 4`.
- OPACITY / WINDOW: `[c  lnSigma  peak (lo hi)]`. `c` is the time of peak visibility in seconds on
  the clip's timeline; `lnSigma` is the natural log of σ in seconds; `peak` is in [0, 1]. The optional
  hard interval `[lo, hi)` makes per-frame content exact (§12b).
- TRANSFORM / SKIN: `[i_1 … i_J  w_1 … w_J]`: driver indices (rounded to the nearest integer), then
  non-negative weights.

**Encodings.** `FLOAT32`, `RANGE` (bits between per-component min and max, as custom attributes
use) or `STEP` (integer multiples of a per-component step, zigzag-coded byte planes). STEP's error is
`step / 2` regardless of the data; RANGE's depends on the spread, so one fast Gaussian coarsens every
other. STEP's high byte planes are almost all zero, which zstd removes. `channelMaxError` reports a
channel's bound, and `spz_info` prints it.

**Point ranges.** A channel holds values for every point, or only `[pointBegin, pointEnd)`. Writers
order points `[moving | static]` and give motion channels the moving range, so static Gaussians cost
nothing. (Custom attributes gained the same optional range.)

## 6. Coordinate conversion

Coordinate conversion runs on load (and save) through `CoordinateConverter`, including
cross-family ±90° rotations about X.

- POSITION / POLYNOMIAL: each xyz triple converts with `convertDirection`, like positions: the
  coefficients transform like positions minus the origin.
- ROTATION / POLYNOMIAL: each 4-tuple goes through `convertQuaternion` (`rotFlipQFunc` or `flipQ`),
  exactly like core rotations. This is correct because the core conversion is linear in `q`:
  `r ⊗ (q + Σ e_k u^k) = r ⊗ q + Σ (r ⊗ e_k) u^k`. Treating the terms as vectors would be wrong:
  they're quaternions, and vector conversion mis-signs pseudovectors under mirroring.
- STEP channels: conversions permute and mix a vector's components, so each vector's components
  take its smallest step. The declared bound never loosens; re-quantizing mixed values on save can
  add up to that bound again.
- OPACITY / WINDOW, SKIN and time values: unaffected.
- Drivers: a rigid transform (R, T) becomes (C R Cᵀ, C T), where C is the conversion's orthogonal
  matrix (axis flips and/or ±90° about X). This keeps skinning consistent with how core rotations
  convert, for proper rotations and reflections alike; tests check that converting then evaluating
  equals evaluating then converting, for every kind of conversion.

After conversion the packer may re-canonicalize `q`'s sign; §3's sign rule applies at that point
too.

## 7. Library behavior

### Loading

- `UnpackOptions::temporal` (`TemporalLoad::Skip` by default, or `Include`). With `Skip` the record's
  extension frames are never decompressed and the cloud loads as its static poster frame; with
  `Include` the record loads whole, or not at all if it fails validation (§4). `customAttributes =
  std::nullopt` (decode everything, for re-saving) implies `Include`.

### Saving

- `saveSpz` validates the record against the cloud's point count (§4 and §5) and **fails the save**
  (logged) if it doesn't fit. An inconsistent in-memory cloud is a caller bug.
- At pack time, any Gaussian whose `q` the packer normalizes or negates has its rotation terms
  multiplied by the same factor (§3). This is a general capability (`SpzRotationCoupledExtension`):
  the quaternion packers report the factor they applied, and `packGaussians` hands every point's
  factor to the records that need it.

### API (C++; Python and JS mirror it)

```cpp
// Record: SpzExtensionTemporalNiantic { clipStart, clipEnd, posterTime, loop, nominalFps,
//   opacityCutoff, channels, driverCount, driverKeyframeTimes, driverTransforms, driverNames }
// Channels: SpzTemporalChannel { target, type, flags, dim, encoding, bits, steps, pointRange,
//   values }

// The record validated against the cloud, with its channels looked up; std::nullopt (logged) means
// static.
std::optional<SpzTemporalModel> bindTemporalModel(const GaussianCloud& cloud);

// A static cloud at time t, without the temporal record. For tests, thumbnails,
// PLY frame export, and as the golden reference for shader implementations.
GaussianCloud evaluateTemporal(const GaussianCloud& cloud, const SpzTemporalModel& model, double t);

// Allocation-free playback: writes only positions, rotations and alphas into a reused frame.
void evaluateTemporalInto(const GaussianCloud& cloud, const SpzTemporalModel& model, double t,
                          SpzTemporalFrame* frame);

// Playback seconds to clip time, by loop mode.
double clipTime(const SpzExtensionTemporalNiantic& clip, double playbackSeconds);

// A point's origin and span (§3), and a channel's declared error bound per component.
SpzTemporalTimeFrame temporalTimeFrame(const SpzTemporalModel& model, size_t point);
std::vector<float> channelMaxError(const SpzTemporalChannel& channel);

// Builds an animated cloud from the generic source form (§9): re-expresses the polynomials in each
// point's time frame, writes each point's state at its origin into the core, orders points
// animated-first, and creates STEP channels at source.precision (error bounds) with point ranges.
std::optional<GaussianCloud> buildTemporalCloud(const SpzTemporalSource& source);
```

Python mirrors these (`evaluate_temporal`, `build_temporal_cloud`, `clip_time`, `is_animated`,
`find_temporal`, `temporal_time_frame`, `channel_max_error`). JavaScript has the record bindings and a `TemporalPlayer` that keeps a clip in the
WASM heap and evaluates frames into views the renderer uploads directly.

## 8. Future: long content, segments and manifest

Not in v1. A v1 file is one clip, normally `[0, duration)`. Past about 5–10 s, long content becomes a
sequence of independent segment files plus a manifest. Each segment uses the same record, with
`clipStart` and `clipEnd` on the global timeline:

- **Why not one file.**
  - Extension frames sit in the last core stream, so every reader downloads all of them.
  - A file has a single `numPoints`, and A's Gaussian count grows with duration.
  - Proper multi-cloud files with a time→byte-range index would need a format v5.
- **Segments.** Each is a self-contained SPZ v4 file, recommended 1–3 s long, with its poster time
  at its midpoint. A Gaussian crossing a boundary is written into both segments. The
  `t ∉ [clipStart, clipEnd)` rule clips it, so switching segments needs no crossfade. Seeking means
  fetching one file.
- **Manifest** (`*.spzseq.json`): duration, loop mode, and an ordered segment list (`uri`,
  `clipStart`, `clipEnd`, `byteLength`, optional hash). It also allows optional *renditions* (point
  budget, SH degree, bit depth) for adaptive streaming, in line with MPEG GSC use case 29 (dynamic
  splat track over HTTP adaptive streaming).

## 9. Import: the temporal PLY convention

Every training codebase already writes PLY, so the generic import path is a documented set of PLY
properties. Per-pipeline adapters are small scripts that rename and reshape; the library does the
re-centering, poster-frame computation, sign canonicalization, and point ordering.

Each Gaussian is described **at its own reference time `t_center`**, the convention STG, FreeTimeGS,
and `.4dgs` already use:

| PLY properties | Meaning | Required |
|---|---|---|
| standard 3DGS (`x y z`, `rot_0..3`, `scale_0..2`, `opacity`, `f_dc_*`, `f_rest_*`) | State at `t_center`; `opacity` (logit, as in 3DGS) is the window's `peak` | yes |
| `t_center` | Time of peak visibility, seconds | yes |
| `t_scale` | ln σ of the opacity window, seconds (log space, like `scale_*`); omit for always-visible | no |
| `motion_0..8` | Position polynomial about `t_center`: linear xyz, quadratic xyz, cubic xyz. Fewer properties mean lower degree | no |
| `omega_0..3` (`omega_4..11` for higher degrees) | Rotation polynomial about `t_center`, quaternion order matching `rot_*` | no |

Clip timing (`clipStart`, `clipEnd`, `posterTime`) comes from the PLY header comments or the import
options.

The conversion into each point's time frame (§3) is exact: with origin `o` and span `s`, re-centering
gives `a'_k = Σ_{j≥k} C(j,k) · a_j · (o − t_center)^{j−k}`, the stored terms are `d_k = a'_k · s^k`, and
the core state is the curve at `o`.

Adapters, in priority order:

| Source | Mapping to the convention | Loss |
|---|---|---|
| **Fudan 4DGS** (Yang et al.; also 4DGS-1K) — first | 4D mean `(μ_xyz, μ_t)` and 4D covariance `Σ` from the two rotor quaternions and 4D scale. Conditional on `t`: `motion = Σ_xyz,t / Σ_tt` (linear), 3D cov `Σ_xyz − Σ_xyz,t Σ_t,xyz / Σ_tt` (eigendecompose → `scale`, `rot`), `t_center = μ_t`, `t_scale = ln √Σ_tt` | Geometry: none. Color: 4D spherindrical harmonics → static SH evaluated at `t_center` (loses time-varying color) |
| **FreeTimeGS** | `μ(t) = μ + v·(t − c)` with Gaussian temporal opacity: `motion_0..2 = v`, `t_center = c`, `t_scale = ln s` | None (verify against the reference code) |
| **STG** (Spacetime Gaussians) | Native names already close: cubic `motion`, linear `omega`, temporal RBF (check its σ convention, e.g. `exp(−(Δt/s)²)` means `σ = s/√2`) | Geometry: none. **Appearance: STG splats neural features decoded by an MLP**, so SH must be fitted (distillation). Training-side tool |
| **Wu 4DGaussians / Deformable-3DGS** (neural) | Sample the deformation at frame times, least-squares fit a cubic per Gaussian | Fit error. B handles this better (§12) |

`tools/4dgs/temporal_ply_to_spz.py` implements this convention in Python.

## 10. Size estimate

Per Gaussian, after zstd, at the builder's default bounds:

- **Core at SH degree 3:** about 65 bytes before zstd (position 9, alpha 1, color 3, scale 3,
  rotation 4, SH 45).
- **STEP channels** store 4 byte planes per component, but levels within a few thousand steps leave
  the top two planes almost entirely zero, so they compress to about the size of 16-bit levels.

| Channels | Bytes per moving Gaussian (≈) | Overhead vs core |
|---|---|---|
| Linear position + window (FreeTimeGS / Fudan) | 6 + 6 = 12 | ≈ +18% |
| Cubic position + linear rotation + window (STG) | 18 + 8 + 6 = 32 | ≈ +49% |

With point ranges, static Gaussians add nothing to the motion channels.

## 11. Rendering guidance (informative)

First target: **PlayCanvas / SuperSplat on WebGPU** (`tools/4dgs/viewer`).

Per-frame work must scale with what's visible, not with the clip. Clips of short-lived Gaussians
(STG-style windows, and especially frame sequences, where each frame's Gaussians are new) can hold
tens of times more Gaussians than any frame shows.

- **Time index (§4).** The record stores, per time bin, the range of points that may be visible.
  - Writers sort points by when they become visible, and bins follow intervals when the clip has
    them, so a frame sequence gets exactly one bin per frame.
  - A renderer evaluates only `range(bin(t))`: `evaluateVisibleTemporal` in C++,
    `TemporalPlayer.evaluateVisible` in JS.
  - It sizes GPU buffers for the largest bin (`maxVisible`) instead of the whole clip.
- **Measured:** on a 164-frame, 16.4M-Gaussian frame sequence (100k per frame), the WASM player
  takes about 3.5 ms per frame to evaluate the visible set and about 1.7 ms to pack it. Before, it
  evaluated all 16.4M points (112 ms) and re-uploaded 656 MB of textures; now it uploads about 4 MB.
- **Sorting:** re-sort every frame, since positions change continuously. WebGPU sorts on the GPU;
  WebGL's CPU sort can lag a frame (PlayCanvas's `.gsaf` player animates only on WebGPU for this
  reason).
- **Next step: GPU evaluation.** Upload the core attributes and channel data once as textures, and
  evaluate the §3 equations in a shader from a time uniform; the time index still picks the bin's
  range. CPU work per frame then drops to setting a uniform. The conformance corpus is the reference
  for shader implementations.

## 12a. Grouped motion: skinning (implemented)

Most real motion is shared: a body part moves as a unit. Instead of every Gaussian carrying its own
curve, a `TRANSFORM` / `SKIN` channel lets Gaussians follow a few hundred **rigid drivers**:

- **Drivers** live in the temporal record: `driverCount` transforms per keyframe, at keyframe times
  shared by every driver, each a translation and rotation from the bind pose. Translation
  interpolates linearly, rotation by shortest-path slerp; times outside the keyframes hold the
  nearest one. Driver data is small (drivers × keyframes × 7 floats) and stored in the record's own
  extension frames.
- **Per Gaussian**, the skin channel holds J ≤ 8 driver indices and weights: about 16 bytes at
  16 bits for J = 4, against 32 for STG-style cubic motion. Per-point channels stack on top, in the
  bind-pose frame, for local detail.
- **Blending** is linear blend skinning (§3), the model glTF skins, USD Skel and SC-GS positions use.
  Rigid groups are the J = 1 case; with points of a group stored together, their identical index and
  weight columns compress to almost nothing.
- **Core attributes are the bind pose.** Blends don't compose, so rebinding to an arbitrary frame
  isn't exact; writers should make the bind pose a real frame (the drivers' identity) and set
  `posterTime` to it, so old readers still see a sensible still.

- **Names.** Drivers may carry labels (bone or object names) for editing, retargeting and picking.
- **Objects.** A rigid object is a driver its Gaussians skin to with weight 1. Membership, in
  static and animated clouds alike, is the `object_id` custom attribute (integers, 0 for
  background), so it needs no new mechanism. `tools/4dgs/spz4d.py` has `set_object_ids` and
  `rigid_objects_clip`.

Importers: `tools/4dgs/scgs_export.py` and `scgs_to_spz.py` turn an SC-GS model into a skinned clip.
SC-GS positions are exactly LBS of its control nodes. Its Gaussian orientations are not: it adds a
separately learned per-node quaternion residual to the stored quaternion (`normalize(q + Σ w δ)`),
and adds a small time-varying offset to linear scales. The converter rotates orientations with the
nodes instead and folds the time-averaged scale offset into the bind pose, then reports position,
rotation and scale error against SC-GS's own output per scene. If that error proves visible, the
format could add an additive-residual rotation mode for skinned channels; measure first.

## 12b. Ideas taken from `.4dgs` (avala-ai)

[avala-ai/4dgs](https://github.com/avala-ai/4dgs) defines a standalone 4DGS container. It shares the
per-Gaussian parametric approach, and several of its ideas carried over:

| `.4dgs` idea | Here |
|---|---|
| Fixed-step quantization with declared error bounds (pitch = 2 × bound) | `STEP` encoding: zigzag-coded byte planes, bound `step / 2`, reported by `channelMaxError` and `spz_info` |
| Velocity precision set per Gaussian from its lifetime | Per-point time frames (§3): terms are displacements over the visible span, so one step bounds every Gaussian's error exactly; also fixes far-extrapolated poster-frame positions |
| Validity windows `[lo, hi)` per Gaussian | Optional interval on opacity windows (dim 5); `ply_sequence_to_spz.py` stores untracked frame sequences exactly |
| Object table, object ids, object tracks | Driver names; `object_id` custom attribute convention; rigid objects as J = 1 skinning |
| Conformance corpus, Kaitai grammar, named refusals | `tests/conformance` (C++ and WASM readers checked), `extensions/spz-v4.ksy` (checked against the layout parser), static-fallback variants |
| Explicit lengths so readers skip what they don't know | Length-prefixed channel descriptors and per-channel `dataBytes` |

Not adopted:
- deflate as a default codec (SPZ v4 is zstd throughout);
- a new container (SPZ files must stay readable by every SPZ reader);
- profiles;
- cameras, audio and sensor rigs (separate extensions if needed).

Their keyframe-delta design (stable Gaussian ids; deaths, updates, births; exact integer deltas so
errors don't accumulate; bounded chain depth) is the template for family B and segments (§8, §12).

## 12. Future: keyframed Gaussians (family B)

A struggles with complex, non-smooth motion such as cloth, hands, fluids, and long non-repeating
movement. It compensates by spawning more short-lived Gaussians. Keyframes (Ex4DGS; tracked
captures such as DualGS and HiFi4G) represent arbitrary trajectories, and piecewise Hermite
keyframes include a single cubic as a special case. Capture-stage pipelines produce B natively, and
it is the natural target for baking neural deformation fields.

### B is a channel type

- **Keyframe time table.** A new optional trailing record field (readers already ignore trailing
  bytes): `u32 keyframeCount, f64 keyframeTimes[keyframeCount]`, or `(start, fps)` for uniform
  spacing. All keyframe channels share it.
- **Channel type `KEYFRAMES = 3`,** valid for every target, with an interpolation mode added to the
  channel (linear, Catmull-Rom/Hermite, slerp).
  - POSITION / KEYFRAMES: `dim = 3K`, converted like positions. Each keyframe is a position offset from the poster
    frame, so static content stays all-zero.
  - ROTATION / KEYFRAMES: `dim = 4K`, converted like rotations, slerp, with the same sign rule as §3.
  - OPACITY / KEYFRAMES: `dim = K`.
- **One target, several channels.** A target may then have one channel per type, as long as their
  point ranges don't overlap: e.g. POSITION / KEYFRAMES on points `[0, k)` and
  POSITION / POLYNOMIAL on `[k, m)`. The "one channel per target" rule in §4 relaxes to "one
  channel per target per point". OPACITY / WINDOW can still apply on top as a visibility envelope.
- **Old v1 readers.** New channels are `REQUIRED`, so v1 readers fall back to the poster frame
  rather than mis-render. That is the reason for the flag.

### Why point ranges matter for B

Ex4DGS's efficiency comes from keyframing only the dynamic subset. Point ranges (decided for
phase 0) give exactly that. Points are ordered `[keyframed | polynomial | static]`, and each
channel covers its range. Renderers then branch on point range instead of per-Gaussian
flags, which stays uniform on GPUs.

### Size, and why segments become mandatory

B grows with the keyframe count `K`. At 30 fps over 2 s with 16-bit position offsets, that's 60
keyframes × 6 bytes = 360 bytes per keyframed Gaussian, about 5× the core. Practical B content
therefore needs:

- **Sparse keyframes** (Ex4DGS-style adaptive spacing, a few per second), not one per frame.
- **Short segments** (§8), which become required rather than optional.
- **Per-segment keyframe tables**, each starting fresh so segments stay independently decodable.
  Cross-segment delta coding is the step after that, and the point where SPZ would converge with
  the MPEG GSC dynamic codecs (keyframe plus residual, as in V³, 4DGCPro, and ResStreamGS). It is
  also where adopting their payloads inside segments could be considered.

## 13. Implementation status

| Phase | Work | Status |
|---|---|---|
| 0. Prerequisites | 16-bit quantization; point range on custom attributes; shared component-major encoding; pack-time rotation factors (`SpzRotationCoupledExtension`) | Done. Custom-attribute PLY round trip replaced by the Python temporal PLY importer |
| 1. Record | Self-contained `SpzExtensionTemporalNiantic` (channels own their data): parse/write/validate, keyed, `TemporalLoad`, coordinate conversion of drivers and channels | Done (`extensions/cc/temporal-niantic.*`) |
| 2. Evaluation | `bindTemporalModel`, `evaluateTemporal(Into)`, `clipTime`, `buildTemporalCloud`; skinning | Done |
| 3. Bindings | Python and WASM/TS; JS `TemporalPlayer` | Done |
| 4. Import and tools | `spz_info`, `spz_to_ply --time`; `tools/4dgs`: synthetic clips, PLY frame export, temporal PLY, frame sequences, objects, Fudan 4DGS, SC-GS | Done; D-NeRF runs pending (need a CUDA machine) |
| 4b. Precision and robustness (§12b) | STEP encoding, per-point time frames, intervals, driver names, object ids, length-prefixed descriptors, conformance corpus, Kaitai grammar | Done |
| 5. PlayCanvas playback | `tools/4dgs/viewer`: `TemporalPlayer` frames pushed into a PlayCanvas `GSplatContainer` (base color only; WebGPU sorts per frame) | Written; not yet checked in a browser |
| 5b. Render cost (§11) | Time index in the record; time-sorted builder; visible-only evaluation (C++, Python, JS); viewer sized for the visible set and decoding once | Done; per-frame cost independent of clip length |
| 6. Segments (after v1) | Manifest spec, reader/writer helpers, segment prefetch | Not started |

Tests:
- `tests/python/test_custom_attribute.py`: point ranges.
- `tests/python/test_temporal.py`: equations, poster frame, fallbacks, skinning, conversion
  commuting with evaluation, encodings and their bounds, time frames, intervals, names, unknown
  encodings.
- `tests/python/test_4dgs_tools.py`: converters against reference math, frame sequences, objects.
- `tests/python/test_conformance.py` and `tools/4dgs/check_conformance.mjs`: the corpus, against the
  C++ and WASM readers.
- `tests/python/test_grammar.py`: the Kaitai grammar against the layout parser.

## 14. Open questions

1. **SC-GS orientations:** is LBS rotation close enough to SC-GS's learned residual on D-NeRF
   scenes (the converter reports it), or does the format need an additive-residual rotation mode?
2. **Manifest:** JSON as proposed, or an SPZ-native index file? (Post-v1.)
3. **Adapters:** in this repo, or a separate tools repo, given training-side dependencies (PyTorch
   for STG distillation)?
4. **Test content:** a real Fudan 4DGS checkpoint to use as a fixture (license permitting), in
   addition to synthetic clips.
