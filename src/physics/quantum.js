// quantum wave function operations and simulation logic
// split-step method for time evolution, wave packet creation, normalization

import { fft3d } from './fft.js';

/**
 * create coordinate array for the simulation grid
 * @param {number} N - grid size
 * @param {number} L - domain size
 * @returns {Float32Array} coordinate array
 */
export function createCoordinateArray(N, L) {
  const arr = new Float32Array(N);
  const s = L / N; 
  const off = -L/2 + s * 0.5;
  for (let i = 0; i < N; i++) {
    arr[i] = off + i * s;
  }
  return arr;
}

function validateGridSize(N) {
  if (!(N >= 2) || ((N & (N - 1)) !== 0)) {
    throw new Error(`grid size N must be a power of two (>=2); got N=${N}`);
  }
}

/**
 * create k-space arrays for kinetic energy operator
 * @param {number} N - grid size
 * @param {number} L - domain size
 * @returns {Object} object containing kx2, ky2, kz2 arrays
 */
export function createKSpaceArrays(N, L) {
  validateGridSize(N);
  const kOf = (n) => {
    const m = (n <= N/2-1) ? n : n - N; 
    return (2*Math.PI/L) * m;
  };
  
  const kx2 = new Float32Array(N);
  const ky2 = new Float32Array(N);
  const kz2 = new Float32Array(N);
  
  for (let i = 0; i < N; i++) {
    const kx = kOf(i), ky = kOf(i), kz = kOf(i);
    kx2[i] = kx * kx;
    ky2[i] = ky * ky;
    kz2[i] = kz * kz;
  }
  
  return { kx2, ky2, kz2 };
}

/**
 * build exponential operators for kinetic energy evolution
 * @param {Float32Array} expK - output array for kinetic exponentials
 * @param {Object} kArrays - k-space arrays from createKSpaceArrays
 * @param {number} N - grid size
 * @param {number} dt - time step
 */
export function buildKineticExponentials(expK, kArrays, N, dt) {
  const { kx2, ky2, kz2 } = kArrays;
  
  for (let z = 0; z < N; z++) {
    const zOff = N*N*z;
    for (let y = 0; y < N; y++) {
      const yOff = zOff + N*y;
      for (let x = 0; x < N; x++) {
        const id = yOff + x;
        const k2 = kx2[x] + ky2[y] + kz2[z];
        const theta = -0.5 * k2 * dt;
        const cr = Math.cos(theta); 
        const ci = Math.sin(theta);
        const j = id << 1; 
        expK[j] = cr; 
        expK[j+1] = ci;
      }
    }
  }
}

/**
 * build exponential operators for potential energy evolution
 * @param {Float32Array} expVh - output array for potential exponentials
 * @param {Float32Array} V - potential array
 * @param {Float32Array} capS2 - absorbing boundary array
 * @param {number} dt - time step
 * @param {number} absorbStrength - absorption strength
 */
export function buildPotentialExponentials(expVh, V, capS2, dt, absorbStrength) {
  const half = -0.5 * dt;
  const len = V.length;
  
  const useCap = !!capS2 && absorbStrength > 0;
  for (let i = 0; i < len; i++) {
    const theta = half * V[i];
    const decay = useCap ? Math.exp(-0.5 * absorbStrength * capS2[i] * dt) : 1;
    const cr = Math.cos(theta) * decay; 
    const ci = Math.sin(theta) * decay;
    const j = i << 1; 
    expVh[j] = cr; 
    expVh[j+1] = ci;
  }
}

/**
 * create absorbing boundary (cap) array
 * @param {number} N - grid size
 * @param {number} absorbFrac - fraction of domain for absorption
 * @returns {Float32Array} cap array
 */
export function createAbsorbingBoundary(N, absorbFrac) {
  const cap = new Float32Array(N*N*N);
  if (!(absorbFrac > 0)) {
    return cap;
  }
  let nAbs = Math.floor(absorbFrac * N);
  nAbs = Math.max(4, nAbs);
  nAbs = Math.min(nAbs, Math.floor((N - 2) / 2));
  if (nAbs <= 0) return cap;
  
  for (let z = 0; z < N; z++) {
    const zOff = N*N*z;
    const dz = Math.max(0, (nAbs - z)/nAbs, (z - (N-1-nAbs))/nAbs);
    
    for (let y = 0; y < N; y++) {
      const yOff = zOff + N*y;
      const dy = Math.max(0, (nAbs - y)/nAbs, (y - (N-1-nAbs))/nAbs);
      const dy2 = dy * dy;
      
      for (let x = 0; x < N; x++) {
        const dxv = Math.max(0, (nAbs - x)/nAbs, (x - (N-1-nAbs))/nAbs);
        const s2 = dxv*dxv + dy2 + dz*dz;
        cap[yOff + x] = s2 * s2;
      }
    }
  }
  
  return cap;
}

/**
 * apply potential energy evolution for half time step
 * @param {Float32Array} psiRe - real part of wave function
 * @param {Float32Array} psiIm - imaginary part of wave function
 * @param {Float32Array} expVh - potential exponentials
 */
export function potentialHalfStep(psiRe, psiIm, expVh) {
  const len = psiRe.length;
  if (psiIm.length !== len) {
    throw new Error(`potentialHalfStep requires psiRe and psiIm to have the same length; got psiRe.length=${len}, psiIm.length=${psiIm.length}`);
  }
  if (expVh.length !== (len << 1)) {
    throw new Error(`potentialHalfStep requires expVh.length === 2 * psi length; got psi length=${len}, expVh.length=${expVh.length}`);
  }
  for (let i = 0; i < len; i++) {
    const pr = psiRe[i], pi = psiIm[i];
    const j = i << 1; 
    const vr = expVh[j], vi = expVh[j+1];
    const nr = pr*vr - pi*vi; 
    const ni = pr*vi + pi*vr;
    psiRe[i] = nr; 
    psiIm[i] = ni;
  }
}

/**
 * apply kinetic energy evolution for full time step
 * @param {Float32Array} psiRe - real part of wave function
 * @param {Float32Array} psiIm - imaginary part of wave function
 * @param {Float32Array} expK - kinetic exponentials
 * @param {number} N - grid size
 * @param {Float32Array} scratchRe - scratch array for fft
 * @param {Float32Array} scratchIm - Scratch array for FFT
 */
export function kineticFullStep(psiRe, psiIm, expK, N, scratchRe, scratchIm) {
  const len = psiRe.length;
  if (psiIm.length !== len) {
    throw new Error(`kineticFullStep requires psiRe and psiIm to have the same length; got psiRe.length=${len}, psiIm.length=${psiIm.length}`);
  }
  const expectedLen = N * N * N;
  if (len !== expectedLen) {
    throw new Error(`kineticFullStep requires psi arrays length N^3; got psi length=${len}, N=${N}`);
  }
  if (expK.length !== (len << 1)) {
    throw new Error(`kineticFullStep requires expK.length === 2 * psi length; got psi length=${len}, expK.length=${expK.length}`);
  }
  if (scratchRe.length < N || scratchIm.length < N) {
    throw new Error(`kineticFullStep requires scratchRe/scratchIm length d N; got scratchRe.length=${scratchRe.length}, scratchIm.length=${scratchIm.length}, N=${N}`);
  }

  // forward FFT to momentum space
  fft3d(psiRe, psiIm, N, false, scratchRe, scratchIm);
  
  // apply kinetic energy operator
  for (let i = 0; i < len; i++) {
    const j = i << 1; 
    const kr = expK[j], ki = expK[j+1];
    const pr = psiRe[i], pi = psiIm[i];
    const nr = pr*kr - pi*ki; 
    const ni = pr*ki + pi*kr;
    psiRe[i] = nr; 
    psiIm[i] = ni;
  }
  
  // inverse FFT back to position space
  fft3d(psiRe, psiIm, N, true, scratchRe, scratchIm);
}

/**
 * perform one time step using split-step method
 * @param {Float32Array} psiRe - real part of wave function
 * @param {Float32Array} psiIm - imaginary part of wave function
 * @param {Float32Array} expVh - potential exponentials
 * @param {Float32Array} expK - kinetic exponentials
 * @param {number} N - grid size
 * @param {Float32Array} scratchRe - scratch array for FFT
 * @param {Float32Array} scratchIm - scratch array for FFT
 */
export function timeStep(psiRe, psiIm, expVh, expK, N, scratchRe, scratchIm) {
  validateGridSize(N);
  potentialHalfStep(psiRe, psiIm, expVh);
  kineticFullStep(psiRe, psiIm, expK, N, scratchRe, scratchIm);
  potentialHalfStep(psiRe, psiIm, expVh);
}

/**
 * add a 3d gaussian wave packet to the wave function
 * @param {Float32Array} psiRe - real part of wave function
 * @param {Float32Array} psiIm - imaginary part of wave function
 * @param {Float32Array} coord - coordinate array
 * @param {number} N - grid size
 * @param {number} cx - center x position
 * @param {number} cy - center y position
 * @param {number} cz - center z position
 * @param {number} sx - width in x direction
 * @param {number} sy - width in y direction
 * @param {number} sz - width in z direction
 * @param {number} kx - wave vector x component
 * @param {number} ky - wave vector y component
 * @param {number} kz - wave vector z component
 * @param {number} scale - overall amplitude
 * @param {Float32Array} gX - scratch array for x gaussian
 * @param {Float32Array} gY - scratch array for y gaussian
 * @param {Float32Array} gZ - scratch array for z gaussian
 * @param {Float32Array} pX - scratch array for x phase
 * @param {Float32Array} pY - scratch array for y phase
 * @param {Float32Array} pZ - scratch array for z phase
 */
export function addPacket3D(...args) {
  let psiRe, psiIm, coord, N, cx, cy, cz, sx, sy, sz, kx, ky, kz, scale, gX, gY, gZ, pX, pY, pZ;
  if (args.length === 1 && args[0] && typeof args[0] === 'object' && ('psiRe' in args[0])) {
    ({ psiRe, psiIm, coord, N, cx, cy, cz, sx, sy, sz, kx, ky, kz, scale = 1, gX, gY, gZ, pX, pY, pZ } = args[0]);
  } else {
    [psiRe, psiIm, coord, N, cx, cy, cz, sx, sy, sz, kx, ky, kz, scale, gX, gY, gZ, pX, pY, pZ] = args;
  }
  if (!(sx > 0 && sy > 0 && sz > 0)) {
    throw new Error(`addPacket3D requires positive widths sx, sy, sz; got sx=${sx}, sy=${sy}, sz=${sz}`);
  }
  const s2x = sx*sx, s2y = sy*sy, s2z = sz*sz;
  
  // compute 1D gaussians and phases
  for (let x = 0; x < N; x++) { 
    const X = coord[x]; 
    gX[x] = Math.exp(-((X-cx)*(X-cx))/(2*s2x)); 
    pX[x] = kx*(X-cx); 
  }
  for (let y = 0; y < N; y++) { 
    const Y = coord[y]; 
    gY[y] = Math.exp(-((Y-cy)*(Y-cy))/(2*s2y)); 
    pY[y] = ky*(Y-cy); 
  }
  for (let z = 0; z < N; z++) { 
    const Z = coord[z]; 
    gZ[z] = Math.exp(-((Z-cz)*(Z-cz))/(2*s2z)); 
    pZ[z] = kz*(Z-cz); 
  }

  // add 3D gaussian packet
  for (let z = 0; z < N; z++) {
    const zOff = N*N*z; 
    const gz = gZ[z]; 
    const pz = pZ[z];
    for (let y = 0; y < N; y++) {
      const yOff = zOff + N*y; 
      const gy = gY[y]; 
      const py = pY[y];
      for (let x = 0; x < N; x++) {
        const id = yOff + x;
        const g = scale * gX[x] * gy * gz;
        const phase = pX[x] + py + pz;
        psiRe[id] += g * Math.cos(phase);
        psiIm[id] += g * Math.sin(phase);
      }
    }
  }
}

/**
 * normalise the wave function to unit probability
 * @param {Float32Array} psiRe - real part of wave function
 * @param {Float32Array} psiIm - imaginary part of wave function
 * @param {number} cellVol - volume of each grid cell
 */
export function renormalize(psiRe, psiIm, cellVol) {
  let sum = 0; 
  const len = psiRe.length;
  for (let i = 0; i < len; i++) { 
    const pr = psiRe[i], pi = psiIm[i]; 
    sum += pr*pr + pi*pi; 
  }
  const scale = 1 / Math.sqrt(sum * cellVol + 1e-30);
  for (let i = 0; i < len; i++) { 
    psiRe[i] *= scale; 
    psiIm[i] *= scale; 
  }
}

/**
 * calculate the total probability (norm) of the wave function
 * @param {Float32Array} psiRe - real part of wave function
 * @param {Float32Array} psiIm - imaginary part of wave function
 * @param {number} cellVol - volume of each grid cell
 * @returns {number} total probability
 */
export function calculateNorm(psiRe, psiIm, cellVol) {
  let norm = 0; 
  const len = psiRe.length;
  for (let i = 0; i < len; i++) { 
    const pr = psiRe[i], pi = psiIm[i]; 
    norm += pr*pr + pi*pi; 
  }
  return norm * cellVol;
}