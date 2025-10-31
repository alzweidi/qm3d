# project TODO

backlog groups potential features and improvements. each item includes a short explanation and feasibility notes. no implementation plan is included yet.

## top 5 to start next

1. ~~profiler HUD (high)~~
   - ~~what: small HUD with frame-time histogram, min/avg/max FPS, dropped-frame counter, toggle via UI.~~
   - ~~feasibility: uses existing RAF timestamps; lightweight UI in React. no rendering-path changes.~~
2. capture UI (high)
   - what: integrate `capture.js` into Controls for WebM/MP4 export (duration, FPS, bitrate), progress, cancel.
   - feasibility: `capture.js` is present; add UI wiring and state. Browser encoder constraints apply.
3. shareable URL state (high)
   - what: encode current params in URL (query/hash) and add a Copy Link button.
   - feasibility: pure client-side serialisation/deserialisation; guarded for backward compatibility.
4. physics workerization (high)
   - what: move `timeStep` loop to a Web Worker; consider OffscreenCanvas for rendering on worker where supported.
   - feasibility: straightforward messaging (postMessage/transferables). feature-flag rollout; fallback to main thread.
5. slice plane tool (medium)
   - what: interactive X/Y/Z cross-sections of the scalar field; toggle between slice and full 3D points.
   - feasibility: reuse existing data buffers; extra shader/material or CPU sampling path.

---

## performance & rendering

- ~~profiler HUD (high)~~
  - ~~what: HUD with frame-time histogram, min/avg/max FPS, dropped-frame counter, quick toggle.~~
  - ~~feasibility: build on existing FPS logic; minimal overhead if sampling window is bounded.~~

- slice plane tool (medium)
  - what: view planar slices along X/Y/Z; adjust plane index/position.
  - feasibility: implement as an alternate geometry/material path or CPU-extracted slice; compatible with current shaders.

- axes/grid overlays (low)
  - what: tiny XYZ axes helper and optional grid overlays to improve orientation.
  - feasibility: Three.js helpers or custom lines; low-risk overlay independent of physics.

## capture, export, and sharing

- capture UI (high)
  - what: add Controls for duration/FPS/bitrate; show progress and allow cancel.
  - feasibility: leverage existing `capture.js`; ensure backpressure handling and UI responsiveness.

- snapshot gallery (medium)
  - what: one-click PNG capture, thumbnail strip, download/share.
  - feasibility: use renderer DOM canvas; memory-bound thumbnail list; optional lazy cleanup.

- save/load states (medium)
  - what: export/import psi and parameters as JSON; include RNG seed for reproducibility.
  - feasibility: serialize Float32Arrays (base64 or binary blob); versioned schema guards.

## UX & accessibility

- FPS settings overlay (medium)
  - what: expose FPS EMA alpha and UI update rate; hide/show and compact HUD modes.
  - feasibility: pure UI/state changes; no physics impact.

- onboarding tour (low)
  - what: guided walkthrough for controls and concepts; contextual tips.
  - feasibility: client-side tour library; non-invasive.

## physics & features

- new potentials (medium)
  - what: double-well, random disorder, tilted plane; pluggable preset interface.
  - feasibility: extend potential arrays and presets; compatible with current split-step scheme.

- camera path recorder (medium)
  - what: record/play camera orbits and parameter timelines; export/import JSON.
  - feasibility: sample controls over time; small timeline player; no engine changes.

- potential editor (low)
  - what: voxel "paint" interface to author V(x,y,z); basic brushes and smoothing.
  - feasibility: CPU-side edits to potential arrays; feature-flagged to mitigate perf risks.

## experiments (feature-flagged)

- physics workerisation (high)
  - what: move `timeStep` loop to a Web Worker; optional OffscreenCanvas.
  - feasibility: message passing for buffers; fallback maintained; progressive rollout.

- WebGPU path (low)
  - what: prototype WebGPU renderer with graceful fallback to WebGL.
  - feasibility: guard behind feature flag; keep Three.js path intact.

- isosurface/volume mode (low)
  - what: optional density isosurfaces or volume raymarching.
  - feasibility: separate render mode with its own shaders; perf-tuned and flag-gated.
