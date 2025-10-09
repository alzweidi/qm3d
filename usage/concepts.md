# Concepts & Reference

This document explains what you’re seeing in the renderer and what each
control does. It’s intentionally plain-language where useful, with short
references for the curious.

## Viewport and camera

Mouse/touch navigation (OrbitControls): left‑drag orbits, wheel/trackpad
zooms (dolly), right‑drag pans. This only moves the camera, not the physics.

In plain words: Drag to look around, scroll to zoom, right‑drag to slide the
view. It doesn’t affect the simulation itself.

Sources: Three.js OrbitControls

## What you’re looking at

Points encode probability density and phase: point size ∝ |ψ| (the square
root of probability density |ψ|²). Colour hue maps the phase arg(ψ) from
−π → π around the colour wheel.

In plain words: Bigger points mean more likely. Regions can look brighter
mostly because many points overlap (additive blending). Colour shows the local
wave’s angle.

Sources:
- MIT OCW: Probability density and current
- MIT OCW Quantum I (lecture notes)

## Simulation grid and time

Grid N (32/64): Number of samples along x, y, z (N³ total). Larger N resolves
finer structure but needs more memory/compute.

In plain words: Higher N = sharper detail, heavier on your GPU/CPU.

Domain L: Physical side length of the cube. Spatial step is dx = L/N. Increasing
L makes the same packet look smaller relative to the box (and vice‑versa).

In plain words: L is the box size. Bigger box → features look smaller on screen.

dt scale: Sets the timestep via dt = dt_scale × dx². Smaller dt reduces
splitting error (Strang split‑step is 2nd order in time).

In plain words: Lower dt scale = more accurate motion but slower.

Steps per frame: How many physics steps are integrated between renders (a
visual smoothness/performance knob).

In plain words: More steps/frame = smoother motion and more simulated time
advanced per frame.

Important (periodic boundaries): The FFT method is periodic. Waves exiting one
face re‑enter at the opposite face unless absorbed (see CAP).

In plain words: Think “wrap‑around world”. Turn on the sponge (CAP) to stop
echoes.

Sources:
- Bao (2003): Time‑splitting spectral discretisations (SIAM)
- MIT 18.336: Spectral methods (periodic domains)

## Packet (initial wave)

σ (sigma): Spatial width of the Gaussian packet. Smaller σ → narrower in space
but broader in momentum (Fourier/uncertainty trade‑off).

In plain words: Tighter blob spreads faster and contains more momenta.

Centre x/y/z: Where you drop the packet inside the box (clamped to ±L/2).

In plain words: Slide these to place the blob.

k0x/k0y/k0z: Central wave‑vector components. For a free particle in these
units, the group velocity ≈ k.

In plain words: Set the direction and speed the blob travels.

Amplitude: Linear scale of the packet you add (superposes with what’s already
there). The norm isn’t force‑normalised after each add. The CAP can drain total
probability.

In plain words: How strong the new blob is. Adding many blobs stacks them. The
sponge can slowly fade them.

Sources:
- MIT Chemistry 5.73: Wavepackets and spreading
- Group velocity primer (UT Austin)

## Absorbing boundaries (CAP)

Width: Thickness of the absorbing “sponge” near the box faces (as a fraction of
the grid). Reduces wrap‑around reflections.

In plain words: How deep the edge‑sponge is.

Strength: How strongly the CAP damps outgoing flux each step (a complex
absorbing potential). Too high can nibble the interior.

In plain words: How aggressive the sponge is. Crank it if you still see echoes,
but don’t overdo it.

Sources:
- Muga et al., Complex absorbing potentials (Phys. Rep. 2004)
- CAP method (arXiv overview)

## Potentials (presets)

Free: V = 0 everywhere (pure transport/dispersion).

In plain words: Nothing pulling on the particle. It just moves and spreads.

Plane barrier: Step barrier on the +x side (V = V0 for x > 0). Shows partial
reflection/transmission.

In plain words: A wall on one side. Some of the wave bounces, some tunnels.

Box well: Cubic region with negative V (can trap bound‑state‑like density).

In plain words: A pit that can hold the blob.

Spherical well: Spherical negative well (isotropic trap).

In plain words: A round pit. Traps from all directions.

Harmonic: Isotropic harmonic V(r) = ½ ω² r². Classic oscillator behaviour
(breathing modes, etc.).

In plain words: Like a spring pulling to the centre. The packet oscillates.
Note: the UI uses a fixed ω = 1 (the function accepts ω but the UI doesn’t
expose it).

Sources: Quantum harmonic oscillator notes (UW/BU‑style)

## Rendering

Density scale: Multiplies point size after auto‑normalisation so faint regions
are still visible. Doesn’t change physics.

In plain words: Visual brightness knob for the blobs.

Phase hue on/off: Toggles colour‑by‑phase. “Off” uses a fixed hue.

In plain words: Switch the rainbow encoding of wave angle on or off.

WebGL point‑size limits: gl_PointSize is clamped by the GPU/driver. The app
queries ALIASED_POINT_SIZE_RANGE and clamps, so very dense spots stop growing
past the hardware cap.

In plain words: If points won’t get any bigger, your GPU hit its limit.

Sources:
- MDN: getParameter and constants (ALIASED_POINT_SIZE_RANGE)
- MDN: WebGL constants

## Runtime controls

Run/Pause: Toggles the Strang split‑step evolution (V half‑step, K full‑step,
V half‑step). Pausing freezes state.

In plain words: Start/stop time without changing the wave.

Reset ψ: Clears the wavefunction and the renderer’s auto‑scaling.

In plain words: Wipe the slate clean.

Add packet: Injects a Gaussian with the current packet settings (adds to
existing ψ).

In plain words: Drop another blob with the sliders you chose.

Presets: Apply one of the potentials above (rebuilds the V‑propagator).

In plain words: Switch the scene’s force field.

Sources: Accuracy and splitting order (SIAM)

## Tips (numerics and physics)

Accuracy: Smaller dt and/or larger N generally improves time‑splitting accuracy
(Strang is 2nd order).

In plain words: For cleanest results: tiny steps, finer grid.

Periodicity and wrap‑around: FFT methods are periodic. Use the CAP to absorb
outgoing waves and avoid re‑entry artefacts.

In plain words: The world loops. The sponge stops ghosts.

Wavepacket intuition: Narrower σ spreads faster and carries a wider momentum
spread. Larger |k0| moves faster (free case).

In plain words: Tight blobs blur quickly. Bigger k makes them travel faster.

Sources:
- MIT OCW: Wavepackets
- Bao (2003) — SIAM
- MIT 18.336 — Spectral (periodic)

## Further reading

- Bao, Jin, Markowich: Time‑splitting spectral methods for (non)linear
  Schrödinger (SIAM 2003; JCP 2002)
- Steven G. Johnson (MIT): Notes on FFT‑based differentiation (spectral/periodic)
- Muga et al. (2004): Complex absorbing potentials, Phys. Rep. 395
- MDN: WebGL constants and limits
