// potential energy functions for quantum wave simulation
// various preset potentials: free space, barriers, wells, harmonic oscillator
//
// units: dimensionless (ℏ = m = 1). energy scale set by V₀ values.
// spatial coordinates: x, y, z ∈ [-L/2, L/2]

/**
 * set potential to zero everywhere (free space).
 * formula: V(r) = 0
 * @param {Float32Array} V - potential array to modify (dimensionless)
 */
export function presetFree(V) {
  V.fill(0);
}

/**
 * create a plane barrier potential (step function at x = 0).
 * formula: V(x,y,z) = V₀ if x > 0, else 0
 * physics: models quantum tunneling and reflection at a potential step.
 * @param {Float32Array} V - potential array to modify (dimensionless)
 * @param {Float32Array} coord - coordinate array
 * @param {number} N - grid size
 * @param {number} V0 - barrier height (default: 4, dimensionless energy units)
 */
export function presetPlaneBarrier(V, coord, N, V0 = 4) {
  for (let z = 0; z < N; z++) {
    const zOff = N*N*z;
    for (let y = 0; y < N; y++) {
      const yOff = zOff + N*y;
      for (let x = 0; x < N; x++) {
        V[yOff + x] = (coord[x] > 0 ? V0 : 0);
      }
    }
  }
}

/**
 * create a cubic potential well centered at origin.
 * formula: V(x,y,z) = V₀ if |x|,|y|,|z| < w/2, else 0  (where w = 0.4L)
 * physics: 3D infinite square well approximation for bound states.
 * @param {Float32Array} V - potential array to modify (dimensionless)
 * @param {Float32Array} coord - coordinate array
 * @param {number} N - grid size
 * @param {number} L - domain size
 * @param {number} V0 - well depth (default: -6, negative = attractive)
 */
export function presetBoxWell(V, coord, N, L, V0 = -6) {
  const w = L * 0.4;
  for (let z = 0; z < N; z++) {
    const zOff = N*N*z;
    for (let y = 0; y < N; y++) {
      const yOff = zOff + N*y;
      for (let x = 0; x < N; x++) {
        const inside = Math.abs(coord[x]) < w/2 && 
                      Math.abs(coord[y]) < w/2 && 
                      Math.abs(coord[z]) < w/2;
        V[yOff + x] = inside ? V0 : 0;
      }
    }
  }
}

/**
 * create a spherical potential well centered at origin.
 * formula: V(r) = V₀ if r < R, else 0  (where R = 0.25L, r = √(x²+y²+z²))
 * physics: spherically symmetric well, models atomic-like potentials.
 * @param {Float32Array} V - potential array to modify (dimensionless)
 * @param {Float32Array} coord - coordinate array
 * @param {number} N - grid size
 * @param {number} L - domain size
 * @param {number} V0 - well depth (default: -8, negative = attractive)
 */
export function presetSphere(V, coord, N, L, V0 = -8) {
  const R = L * 0.25; 
  const R2 = R * R;
  for (let z = 0; z < N; z++) {
    const zOff = N*N*z; 
    const z2 = coord[z] * coord[z];
    for (let y = 0; y < N; y++) {
      const yOff = zOff + N*y; 
      const y2 = coord[y] * coord[y];
      for (let x = 0; x < N; x++) {
        const r2 = coord[x] * coord[x] + y2 + z2;
        V[yOff + x] = (r2 < R2 ? V0 : 0);
      }
    }
  }
}

/**
 * create a 3D isotropic harmonic oscillator potential.
 * formula: V(r) = ½ω²r² = ½ω²(x² + y² + z²)
 * physics: exactly solvable; ground state energy E₀ = (3/2)ℏω = 1.5ω (with ℏ=1).
 * eigenstates are Hermite-Gaussian functions.
 * @param {Float32Array} V - potential array to modify (dimensionless)
 * @param {Float32Array} coord - coordinate array
 * @param {number} N - grid size
 * @param {number} omega - angular frequency (default: 1, dimensionless)
 */
export function presetHarmonic(V, coord, N, omega = 1) {
  for (let z = 0; z < N; z++) {
    const zOff = N*N*z; 
    const z2 = coord[z] * coord[z];
    for (let y = 0; y < N; y++) {
      const yOff = zOff + N*y; 
      const y2 = coord[y] * coord[y];
      for (let x = 0; x < N; x++) {
        const r2 = coord[x] * coord[x] + y2 + z2;
        V[yOff + x] = 0.5 * omega * omega * r2;
      }
    }
  }
}

/**
 * get all available potential presets
 * @returns {Object} map of preset names to functions
 */
export const POTENTIAL_PRESETS = {
  free: presetFree,
  planeBarrier: presetPlaneBarrier,
  boxWell: presetBoxWell,
  sphere: presetSphere,
  harmonic: presetHarmonic
};