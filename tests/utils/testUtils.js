// Shared test utilities to reduce duplication across unit tests
// Keep helpers here so SonarQube doesn't flag duplicate blocks in tests
import { fft1d, fft1d_p, makeFFTPlan, fft3d } from '../../src/physics/fft.js';
import { createCoordinateArray, createKSpaceArrays } from '../../src/physics/quantum.js';

/**
 * Linear index for 3D -> 1D mapping used in physics arrays
 */
export function lin(x, y, z, N) {
  return x + N * (y + N * z);
}

/**
 * RMS relative error between two complex arrays (re, im pairs)
 */
export function rmsRelativeErrorComplex(aRe, aIm, bRe, bIm) {
  let num = 0;
  let den = 0;
  const n = aRe.length;
  for (let i = 0; i < n; i++) {
    const dr = aRe[i] - bRe[i];
    const di = aIm[i] - bIm[i];
    num += dr * dr + di * di;
    const ar = bRe[i];
    const ai = bIm[i];
    den += ar * ar + ai * ai;
  }
  return Math.sqrt(num / Math.max(den, 1e-30));
}

/**
 * Fill real/imag arrays with random values in [-1, 1].
 */
export function fillRandomPair(re, im) {
  const n = re.length;
  for (let i = 0; i < n; i++) {
    re[i] = Math.random() * 2 - 1;
    im[i] = Math.random() * 2 - 1;
  }
  return [re, im];
}

/**
 * Expect planned FFT equals canonical FFT for forward or inverse case.
 * direction: 'forward' | 'inverse'
 */
export function expectFft1dPlanEquivalence(N, direction = 'forward') {
  const reA = new Float32Array(N);
  const imA = new Float32Array(N);
  fillRandomPair(reA, imA);

  const reB = reA.slice();
  const imB = imA.slice();

  const plan = makeFFTPlan(N);

  if (direction === 'forward') {
    fft1d(reA, imA, false);
    fft1d_p(reB, imB, false, plan);
  } else {
    // prepare with a forward first, then compare inverse paths
    fft1d(reA, imA, false);
    fft1d(reB, imB, false);
    fft1d(reA, imA, true);
    fft1d_p(reB, imB, true, plan);
  }

  for (let i = 0; i < N; i++) {
    if (Math.abs(reA[i] - reB[i]) > 1e-10 || Math.abs(imA[i] - imB[i]) > 1e-10) {
      throw new Error('Planned FFT differs from canonical FFT');
    }
  }
}

/**
 * Compute rms error of a 3D FFT round-trip (forward then inverse).
 */
export function compute3dRoundTripError(N) {
  const sz = N * N * N;
  const re = new Float32Array(sz);
  const im = new Float32Array(sz);
  fillRandomPair(re, im);
  const r0 = re.slice();
  const i0 = im.slice();
  const sRe = new Float32Array(N);
  const sIm = new Float32Array(N);
  fft3d(re, im, N, false, sRe, sIm);
  fft3d(re, im, N, true, sRe, sIm);
  return rmsRelativeErrorComplex(re, im, r0, i0);
}

/**
 * Allocate packet scratch arrays used by addPacket3D in multiple tests.
 */
export function allocPacketScratch(N) {
  return {
    gX: new Float32Array(N),
    gY: new Float32Array(N),
    gZ: new Float32Array(N),
    pX: new Float32Array(N),
    pY: new Float32Array(N),
    pZ: new Float32Array(N),
  };
}

/**
 * Create a minimal simulation context for testing.
 * Returns { N, L, dx, cellVol, size, coord, psiRe, psiIm, V, kArrays, sRe, sIm }
 */
export function createTestSimContext(N = 8, L = 10) {
  const dx = L / N;
  const cellVol = dx * dx * dx;
  const size = N * N * N;
  const coord = createCoordinateArray(N, L);
  const psiRe = new Float32Array(size);
  const psiIm = new Float32Array(size);
  const V = new Float32Array(size);
  const kArrays = createKSpaceArrays(N, L);
  const sRe = new Float32Array(N);
  const sIm = new Float32Array(N);
  return { N, L, dx, cellVol, size, coord, psiRe, psiIm, V, kArrays, sRe, sIm };
}

/**
 * Create exponential arrays for simulation testing.
 * Returns { expK, expVh }
 */
export function createTestExponentials(size) {
  return {
    expK: new Float32Array(2 * size),
    expVh: new Float32Array(2 * size),
  };
}
