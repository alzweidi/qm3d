// potential energy functions for quantum wave simulation
// various preset potentials: free space, barriers, wells, harmonic oscillator

/**
 * set potential to zero everywhere (free space)
 * @param {Float32Array} V - potential array to modify
 */
export function presetFree(V) {
  V.fill(0);
}

/**
 * create a plane barrier potential at x > 0
 * @param {Float32Array} V - potential array to modify
 * @param {Float32Array} coord - coordinate array
 * @param {number} N - grid size
 * @param {number} V0 - barrier height (default: 4)
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
 * create a cubic potential well in the center
 * @param {Float32Array} V - potential array to modify
 * @param {Float32Array} coord - coordinate array
 * @param {number} N - grid size
 * @param {number} L - domain size
 * @param {number} V0 - well depth (default: -6)
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
 * create a spherical potential well
 * @param {Float32Array} V - potential array to modify
 * @param {Float32Array} coord - coordinate array
 * @param {number} N - grid size
 * @param {number} L - domain size
 * @param {number} V0 - well depth (default: -8)
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
 * create a 3d harmonic oscillator potential
 * @param {Float32Array} V - potential array to modify
 * @param {Float32Array} coord - coordinate array
 * @param {number} N - grid size
 * @param {number} omega - oscillator frequency (default: 1)
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