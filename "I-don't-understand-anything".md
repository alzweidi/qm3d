# viewport and camera

mouse/touch navigation (orbit controls): left-drag orbits, wheel/trackpad zooms (dolly), right-drag pans. this only moves the camera, not the physics.

*in plain words:* drag to look around, scroll to zoom, right-drag to slide the view. it doesn't affect the simulation itself.

**sources:** three.js orbitcontrols

## what you're looking at

points encode probability density and phase: point size \( \propto |\psi| \) (that's √ of probability density \( |\psi|^2 \)). colour hue maps the phase \( \arg(\psi) \) from \( -\pi \to \pi \) around the colour wheel.

*in plain words:* bigger points mean more likely. regions may look brighter mainly because many points overlap (additive blending). colour shows the local wave's angle.

**sources:**

- mit ocw: probability density and current
- mit ocw quantum i (lecture notes)

## simulation grid and time

**grid n (32/64):** number of samples along x, y, z (\( n^3 \) total). larger n resolves finer structure but needs more memory/compute.

*in plain words:* higher n = sharper detail, but heavier on your gpu/cpu.

**domain l:** physical side length of the cube. spatial step is \( dx = l/n \). increasing l makes the same packet look smaller relative to the box (and vice-versa).

*in plain words:* l is the box size. bigger box means features look smaller on screen.

**dt scale:** sets the time step via \( dt = dt_{scale} \times dx^2 \). smaller dt reduces splitting error (strang split-step is 2nd order in time).

*in plain words:* lower dt scale = more accurate motion but slower.

**steps per frame:** how many physics steps you integrate between renders (visual smoothness/perf knob).

*in plain words:* more steps/frame = smoother motion and more simulated time advanced between renders.

**important (periodic boundaries):** the fft method is periodic. waves exiting one face re-enter at the opposite face unless absorbed (see cap).

*in plain words:* think "wrap-around world". turn on the sponge (cap) to stop echoes.

**sources:**

- bao 2003: time-splitting spectral discretisations (siam)
- mit 18.336: spectral methods (periodic domains)

## packet (initial wave)

**σ (sigma):** spatial width of the gaussian packet. smaller σ means narrower in space but broader in momentum (fourier/uncertainty trade-off).

*in plain words:* tighter blob spreads faster and contains more momenta.

**centre x/y/z:** where you drop the packet inside the box (clamped to \( \pm l/2 \)).

*in plain words:* slide these to place the blob.

**k0x/k0y/k0z:** central wave-vector components. for a free particle in these units, the group velocity \( \approx \mathbf{k} \).

*in plain words:* set the direction and speed the blob travels.

**amplitude:** linear scale of the packet you add (superposes with what's already there). norm isn't force-normalised after each add. the cap can drain total probability.

*in plain words:* how strong the new blob is. adding many blobs stacks them. the sponge can slowly fade them.

**sources:**

- mit chemistry 5.73: wavepackets and spreading
- group velocity primer (ut austin)

## absorbing boundaries (cap)

**width:** thickness of the absorbing "sponge" near the box faces (as a fraction of the grid). reduces wrap-around reflections.

*in plain words:* how deep the edge-sponge is.

**strength:** how strongly the cap damps outgoing flux each step (a complex absorbing potential). too high can nibble the interior.

*in plain words:* how aggressive the sponge is. crank it if you still see echoes, but don't overdo it.

**sources:**

- muga et al., complex absorbing potentials (phys. rep. 2004)
- cap method (arxiv overview)

## potentials (presets)

**free:** \( v = 0 \) everywhere (pure transport/dispersion).

*in plain words:* nothing pulling on the particle. it just moves and spreads.

**plane barrier:** step barrier on the +x side (\( v = v_0 \) for x>0). shows partial reflection/transmission.

*in plain words:* a wall on one side. some of the wave bounces, some tunnels.

**box well:** cubic region with negative v (can trap bound-state-like density).

*in plain words:* a pit that can hold the blob.

**spherical well:** spherical negative well (isotropic trap).

*in plain words:* a round pit. traps from all directions.

**harmonic:** isotropic harmonic \( v(r) = \frac{1}{2} \omega^2 r^2 \). classic oscillator behaviour (breathing modes etc.).

*in plain words:* like a spring pulling to the centre. the blob jiggles. harmonic uses a fixed \( \omega = 1 \) in the ui (the function accepts \( \omega \) but the ui doesn't expose it).

**sources:** quantum harmonic oscillator notes (uw/bu-style)

## rendering

**density scale:** multiplies point size after auto-normalisation so faint regions are still visible. doesn't change physics.

*in plain words:* visual brightness knob for the blobs.

**phase hue on/off:** toggles colour-by-phase. off uses a fixed hue.

*in plain words:* switch the rainbow encoding of wave angle on or off.

**webgl point-size limits:** gl_pointsize is clamped by the gpu/driver. the app queries aliased_point_size_range and clamps, so very dense spots stop growing past the hardware cap.

*in plain words:* if points won't get any bigger, your gpu has hit its limit.

**sources:**

- mdn: getparameter and constants (aliased_point_size_range)
- mdn: webgl constants

## runtime buttons

**run/pause:** toggles the strang split-step evolution (v half-step, k full-step, v half-step). pausing freezes state.

*in plain words:* start/stop time without changing the wave.

**reset ψ:** clears the wavefunction and the renderer's auto-scaling.

*in plain words:* wipe the slate clean.

**add packet:** injects a gaussian with the current packet settings (adds to existing ψ).

*in plain words:* drop another blob with the sliders you chose.

**presets:** apply one of the potentials above (rebuilds the v-propagator).

*in plain words:* switch the scene's force field.

**sources:** accuracy and splitting order (siam)

## tips (numerics and physics)

**accuracy:** smaller dt and/or larger n means better time-splitting accuracy (strang is 2nd order).

*in plain words:* for cleanest results: tiny steps, finer grid.

**periodicity and wrap-around:** fft methods are periodic. use the cap to absorb outgoing waves and avoid re-entry artefacts.

*in plain words:* the world loops. the sponge stops ghosts.

**wavepacket intuition:** narrower σ spreads faster and carries a wider momentum spread. larger \( |\mathbf{k}_0| \) moves faster (free case).

*in plain words:* tight blobs blur quickly. bigger k makes them travel faster.

**sources:**

- mit ocw wavepackets
- bao 2003 (siam)
- mit 18.336 spectral (periodic)

## further reading (reputable)

bao, jin, and markowich: time-splitting spectral methods for (non)linear schrödinger (siam and jcp overviews: siam 2003, jcp 2002)

steven g. johnson (mit): notes on fft-based differentiation (spectral/periodic context): pdf

cap background: muga et al., complex absorbing potentials, phys. rep. 395 (2004): abstract. recent overview: arxiv:1310.1872

webgl point-size and limits: mdn getparameter, mdn constants
