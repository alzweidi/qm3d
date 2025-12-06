# qm3d: interactive 3d quantum wave engine

&#x20;  &#x20;

small, fast, and fairly pretty. a split‑step fourier (strang) tdse solver in the browser with a point‑cloud renderer. built to poke quantum intuition and stress‑test my brain.

> colour encodes phase arg(ψ); point size encodes |ψ|. orbit with the mouse, add packets, flip potentials, watch it breathe.

---

## demo

- live: [https://qm3d.netlify.app/](https://qm3d.netlify.app/)

---

## what this is

- 3d gaussian packet injection with steerable centre, width, k‑vector, and amplitude
- presets: free space, plane barrier, box well, spherical well, harmonic
- absorbing boundaries (cap) to tame fft periodic wrap‑around
- real‑time view: |ψ|² as density; phase as hue; orbit/zoom/pan via orbitcontrols
- responsive ui with sliders + sensible clamps (k respects nyquist)

*note:* imports and paths are plain js/esm; see `engine.jsx` for usage.

---

## why i built this

summer boredom + wanting to relearn some physics. i saw a few open‑source 2d tdse toys, built my own 2d version in another repo (more art than classroom), then went hunting for a good open‑source 3d one. didn’t find one, so i made one and i’m pushing it towards “proper” while keeping it tiny.

---

## install & run

prereqs: node 18+ recommended.

```bash
npm i
npm run dev
```

open the printed localhost url. vite will hot‑reload when you tweak code.

---

## controls (cheat sheet)

- **grid n (32/64):** samples per axis (n³ total). bigger n ⇒ more detail, more gpu/ram.
- **domain l:** cube side length. spatial step is `dx = l / n`.
- **dt scale:** time step via `dt = dt_scale × dx²` (stable + grid‑independent feel).
- **steps/frame:** extra physics steps between renders.
- **σ:** gaussian width. smaller σ spreads faster (uncertainty tax).
- **centre (x/y/z):** packet drop point (clamped to ±l/2).
- **k0x/y/z:** initial k. clamped inside safe band.
- **amplitude:** linear scale before normalisation.
- **cap width/strength:** damping sponge near faces. too strong chews interior.
- **phase hue:** toggle phase colouring.

---

## under the hood (short + honest)

- **solver:** split‑step spectral (strang): `e^{-iVΔt/2} → fft → e^{-iKΔt} → ifft → e^{-iVΔt/2}`
- **k‑grid:** centred using fft ordering; `k = 2π m / L` with `m ∈ [−n/2,…,n/2−1]`.
- **kinetic:** `θ_k = −½ k² dt` baked into a complex exp table per cell.
- **potential:** half‑step phase with optional cap decay `exp(−½·strength·s²·dt)`.
- **cap:** smooth polynomial in a shell of width `absorb_frac × n` on each face (set width or strength to zero to disable).
- **normalisation:** renormalise after packet add; thereafter norm naturally drifts only where cap damps outgoing flux.
- **boundaries:** periodic because fft; the cap exists to hide that.
- **units:** nondimensional (ħ = m = 1). `x` is in `[−l/2,l/2]`; energies match your chosen `V` scale.

### algorithm complexity

| operation | complexity | notes |
|-----------|------------|-------|
| 3d fft | O(N³ log N) | radix‑2 cooley‑tukey per axis |
| time step | O(N³ log N) | 2× fft + O(N³) pointwise ops |
| potential build | O(N³) | one‑time precomputation |
| memory | O(N³) | ~8 float32 arrays of size N³ |

for N=32: ~32K cells, ~2ms/step on modern CPUs. for N=64: ~262K cells, ~20ms/step.

---

## numerics defaults (good starting point)

- `n = 32`, `l = 10`, `σ ≈ 0.6`, `dt_scale ≈ 0.08` ⇒ `dt ≈ 0.08 dx²`
- safe k: target \~8 points/λ → `|k| ≲ min(0.6·π/dx, 2π/(8·dx))`
- cap: `width ≈ 15%` of box; `strength ≈ 3`

if you see sparkle or aliasing, drop `dt_scale`, or lower `k0*`, or bump `n`.

---

## potentials (built‑ins)

- **free:** `v = 0`
- **plane barrier:** potential step for `x > 0`
- **box well:** cubic well at centre (depth configurable in code)
- **spherical well:** radius \~ `0.25 l`
- **harmonic:** `v = ½ ω² r²` (ω=1 default)

\*add your own in \**`potentials.js`*

---

## performance notes

- points render with additive sprites; auto‑exposure on density keeps things legible.
- point size scales with distance and dpr; we cap at the gpu’s `ALIASED_POINT_SIZE_RANGE`.
- try n=32 for smoothness; n=64 is doable on decent gpus, but heavier.

---

## troubleshooting

- **energy blow‑up or ringing:** reduce `dt_scale`; check `k0*` within safe band.
- **wrap‑around ghosts:** widen the cap or increase strength a touch.
- **everything looks tiny:** your `l` is large; remember `dx = l/n`.
- **phase looks static:** toggle phase hue; stationary states do that.

---

## roadmap (rough)

- i'll keep going and see where it leads.

---

## dev notes

- react + three.js + vite + tailwind.
- hand‑rolled 1d/3d fft with round‑trip self‑tests. k‑space arrays precomputed; exponentials cached and rebuilt when `dt`, `v`, or cap changes.
- shaders are tiny: vertex sets sprite size from amplitude; fragment maps phase → hue; additive blending.

---

---

## acknowledgements

textbook split‑step method; shout‑out to all the fft giants we stand on. any bugs are mine.

---

## references

- feit, fleck, & steiger (1982): split‑operator methods for the time‑dependent schrödinger equation.
- tanguy, leforestier, & kosloff (various): spectral/fft approaches to tdse in confined domains.
- good fft primers: frigo & johnson (fftw), press et al. *numerical recipes*.

---

## want to cite?

> “abedalaziz alzweidi, qm3d: a tiny 3d quantum wave engine in the browser (2025), [https://github.com/alzweidi/qm3d](https://github.com/alzweidi/qm3d).”

---

## licence

educational use licence (eul) — non‑commercial educational and research use only. see [`LICENSE`](./LICENSE) for the full terms and contact info.

---

## notes for reviewers&#x20;

if you’re peeking at the numerics: the readme keeps it light; see code comments for exact k‑indexing, dt scaling, and cap definition.

---

## and don't forget to star if you like it
