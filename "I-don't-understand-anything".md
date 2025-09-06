# Viewport & Camera

Mouse / touch navigation (orbit controls): left‑drag orbits, wheel/trackpad zooms (dolly), right‑drag pans. this only moves the camera, not the physics.

*In plain words:* drag to look around, scroll to zoom, right‑drag to slide the view; it doesn't affect the simulation itself.

**Sources:** three.js orbitcontrols

# What You're Looking At

Points encode probability density and phase: point size \( \propto |\psi| \) (that's √ of probability density \( |\psi|^2 \)); colour hue maps the phase \( \arg(\psi) \) from \( -\pi \to \pi \) around the colour wheel.

*In plain words:* bigger points mean more likely; regions may look brighter mainly because many points overlap (additive blending). Colour shows the local wave's angle.

**Sources:**
- MIT OCW: probability density & current
- MIT OCW Quantum I — lecture notes

# Simulation Grid & Time

**Grid n (32/64):** number of samples along x, y, z (\( n^3 \) total). larger n resolves finer structure but needs more memory/compute.

*In plain words:* higher n = sharper detail, but heavier on your GPU/CPU.

**Domain l:** physical side length of the cube; spatial step is \( dx = l/n \). increasing l makes the same packet look smaller relative to the box (and vice-versa).

*In plain words:* l is the box size; bigger box → features look smaller on screen.

**dt scale:** sets the time step via \( dt = dt_scale \times dx^2 \). smaller dt reduces splitting error (strang split-step is 2nd order in time).

*In plain words:* lower dt scale = more accurate motion but slower.

**Steps per frame:** how many physics steps you integrate between renders (visual smoothness/perf knob).

*In plain words:* more steps/frame = smoother motion and more simulated time advanced between renders.

**Important — Periodic boundaries:** the FFT method is periodic; waves exiting one face re-enter at the opposite face unless absorbed (see cap).

*In plain words:* think "wrap-around world"; turn on the sponge (cap) to stop echoes.

**Sources:**
- Bao 2003 — time-splitting spectral discretisations (SIAM)
- MIT 18.336 — spectral methods (periodic domains)

# Packet (Initial Wave)

**σ (sigma):** spatial width of the gaussian packet. smaller σ → narrower in space but broader in momentum (Fourier/uncertainty trade-off).

*In plain words:* tighter blob spreads faster and contains more momenta.

**Centre x / y / z:** where you drop the packet inside the box (clamped to \( \pm l/2 \)).

*In plain words:* slide these to place the blob.

**k0x / k0y / k0z:** central wave-vector components; for a free particle in these units, the group velocity \( \approx \mathbf{k} \).

*In plain words:* set the direction and speed the blob travels.

**Amplitude:** linear scale of the packet you add (superposes with what's already there). norm isn't force-normalised after each add; the cap can drain total probability.

*In plain words:* how strong the new blob is; adding many blobs stacks them; the sponge can slowly fade them.

**Sources:**
- MIT Chemistry 5.73 — wavepackets & spreading
- Group velocity primer (UT Austin)

# Absorbing Boundaries (Cap)

**Width:** thickness of the absorbing "sponge" near the box faces (as a fraction of the grid). reduces wrap-around reflections.

*In plain words:* how deep the edge-sponge is.

**Strength:** how strongly the cap damps outgoing flux each step (a complex absorbing potential). too high can nibble the interior.

*In plain words:* how aggressive the sponge is; crank it if you still see echoes, but don't overdo it.

**Sources:**
- Muga et al., complex absorbing potentials (Phys. Rep. 2004)
- CAP method (arxiv overview)

# Potentials (Presets)

**Free:** \( v = 0 \) everywhere (pure transport/dispersion).

*In plain words:* nothing pulling on the particle; it just moves and spreads.

**Plane barrier:** step barrier on the +x side (\( v = v_0 \) for x>0). shows partial reflection/transmission.

*In plain words:* a wall on one side — some of the wave bounces, some tunnels.

**Box well:** cubic region with negative v (can trap bound-state-like density).

*In plain words:* a pit that can hold the blob.

**Spherical well:** spherical negative well (isotropic trap).

*In plain words:* a round pit; traps from all directions.

**Harmonic:** isotropic harmonic \( v(r) = \frac{1}{2} \omega^2 r^2 \). classic oscillator behaviour (breathing modes etc.).

*In plain words:* like a spring pulling to the centre; the blob jiggles. Harmonic uses a fixed \( \omega = 1 \) in the UI (the function accepts \( \omega \) but the UI doesn't expose it).

**Sources:** quantum harmonic oscillator notes (UW / BU-style)

# Rendering

**Density scale:** multiplies point size after auto-normalisation so faint regions are still visible; doesn't change physics.

*In plain words:* visual brightness knob for the blobs.

**Phase hue on/off:** toggles colour-by-phase; off uses a fixed hue.

*In plain words:* switch the rainbow encoding of wave angle on or off.

**WebGL point-size limits:** gl_PointSize is clamped by the GPU/driver; the app queries ALIASED_POINT_SIZE_RANGE and clamps, so very dense spots stop growing past the hardware cap.

*In plain words:* if points won't get any bigger, your GPU has hit its limit.

**Sources:**
- MDN: getParameter & constants (ALIASED_POINT_SIZE_RANGE)
- MDN: WebGL constants

# Runtime Buttons

**Run / Pause:** toggles the Strang split-step evolution (v half-step → k full-step → v half-step). pausing freezes state.

*In plain words:* start/stop time without changing the wave.

**Reset ψ:** clears the wavefunction and the renderer's auto-scaling.

*In plain words:* wipe the slate clean.

**Add packet:** injects a gaussian with the current packet settings (adds to existing ψ).

*In plain words:* drop another blob with the sliders you chose.

**Presets:** apply one of the potentials above (rebuilds the v-propagator).

*In plain words:* switch the scene's force field.

**Sources:** accuracy & splitting order — SIAM

# Tips (Numerics & Physics)

**Accuracy:** smaller dt and/or larger n → better time-splitting accuracy (Strang is 2nd order).

*In plain words:* for cleanest results: tiny steps, finer grid.

**Periodicity & wrap-around:** FFT methods are periodic; use the cap to absorb outgoing waves and avoid re-entry artefacts.

*In plain words:* the world loops — the sponge stops ghosts.

**Wavepacket intuition:** narrower σ spreads faster and carries a wider momentum spread; larger \( |\mathbf{k}_0| \) moves faster (free case).

*In plain words:* tight blobs blur quickly; bigger k makes them travel faster.

**Sources:**
- MIT OCW wavepackets
- Bao 2003 — SIAM
- MIT 18.336 spectral (periodic)

# Further Reading (Reputable)

Bao, Jin, & Markowich: time-splitting spectral methods for (non)linear Schrödinger — SIAM & JCP overviews: SIAM 2003, JCP 2002

Steven G. Johnson (MIT): notes on FFT-based differentiation (spectral/periodic context): PDF

CAP background: Muga et al., complex absorbing potentials, Phys. Rep. 395 (2004): abstract; recent overview: arXiv:1310.1872

WebGL point-size and limits: MDN getParameter, MDN constants
