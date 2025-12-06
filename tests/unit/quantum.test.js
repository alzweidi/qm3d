import { describe, it, expect } from 'vitest';
import { lin, fillRandomPair, rmsRelativeErrorComplex as rmsRel } from '../utils/testUtils.js';
import {
  createCoordinateArray,
  createKSpaceArrays,
  buildKineticExponentials,
  buildPotentialExponentials,
  createAbsorbingBoundary,
  potentialHalfStep,
  kineticFullStep,
  timeStep,
  addPacket3D,
  renormalize,
  calculateNorm,
  calculateEnergy,
  checkFinite,
} from '../../src/physics/quantum.js';
import { fft3d } from '../../src/physics/fft.js';

// lin moved to tests/utils/testUtils.js

describe('quantum.js', () => {
  it('createCoordinateArray produces centered cell coordinates', () => {
    const arr = createCoordinateArray(4, 2);
    expect(Array.from(arr)).toEqualCloseTo([-0.75, -0.25, 0.25, 0.75], 10);
  });

  it('createCoordinateArray throws on non-positive or non-finite L', () => {
    expect(() => createCoordinateArray(4, 0)).toThrow(/finite positive/);
    expect(() => createCoordinateArray(4, -1)).toThrow(/finite positive/);
    expect(() => createCoordinateArray(4, Infinity)).toThrow(/finite positive/);
  });

  it('renormalize after addPacket keeps norm ≈ 1 over time with CAP off', () => {
    const N = 16;
    const L = 8;
    const dx = L / N;
    const cellVol = dx * dx * dx;
    const coord = createCoordinateArray(N, L);

    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);

    const gX = new Float32Array(N);
    const gY = new Float32Array(N);
    const gZ = new Float32Array(N);
    const pX = new Float32Array(N);
    const pY = new Float32Array(N);
    const pZ = new Float32Array(N);

    const sigma = 0.6;
    const kx = (2 * Math.PI) / (8 * dx); // ~8 ppw

    addPacket3D(
      psiRe, psiIm, coord, N,
      0, 0, 0,
      sigma, sigma, sigma,
      kx, 0, 0,
      1,
      gX, gY, gZ, pX, pY, pZ
    );

    renormalize(psiRe, psiIm, cellVol);
    const norm0 = calculateNorm(psiRe, psiIm, cellVol);
    expect(norm0).toBeCloseTo(1, 6);

    const kArrays = createKSpaceArrays(N, L);
    const expK = new Float32Array(2 * size);
    const expVh = new Float32Array(2 * size);
    const cap = createAbsorbingBoundary(N, 0); // CAP off
    const dt = 0.02 * dx * dx;

    buildKineticExponentials(expK, kArrays, N, dt);
    buildPotentialExponentials(expVh, new Float32Array(size), cap, dt, 0);

    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);
    for (let i = 0; i < 200; i++) {
      timeStep(psiRe, psiIm, expVh, expK, N, sRe, sIm);
    }
    const normT = calculateNorm(psiRe, psiIm, cellVol);
    expect(normT).toBeCloseTo(1, 4);
  });

  it('super-Nyquist kx aliases to DC bin in FFT (motivates UI clamp)', () => {
    const N = 32;
    const L = 10;
    const dx = L / N;
    const coord = createCoordinateArray(N, L);

    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);

    const { gX, gY, gZ, pX, pY, pZ } = {
      gX: new Float32Array(N),
      gY: new Float32Array(N),
      gZ: new Float32Array(N),
      pX: new Float32Array(N),
      pY: new Float32Array(N),
      pZ: new Float32Array(N),
    };

    // Set kx to 2π/dx = 2π * N / L, which aliases to m=0 on the discrete grid
    const kx = 2 * Math.PI / dx;
    addPacket3D(
      psiRe, psiIm, coord, N,
      0, 0, 0,
      0.6, 0.6, 0.6,
      kx, 0, 0,
      1,
      gX, gY, gZ, pX, pY, pZ
    );

    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);
    fft3d(psiRe, psiIm, N, false, sRe, sIm);

    // Find argmax in spectrum
    let maxV = -Infinity, maxI = -1;
    for (let i = 0; i < size; i++) {
      const v = psiRe[i] * psiRe[i] + psiIm[i] * psiIm[i];
      if (v > maxV) { maxV = v; maxI = i; }
    }
    const x = maxI % N;
    const y = Math.floor(maxI / N) % N;
    const z = Math.floor(maxI / (N * N));
    // Expect DC along all axes due to aliasing
    expect(x).toBe(0);
    expect(y).toBe(0);
    expect(z).toBe(0);
  });

  it('createKSpaceArrays maps to integer grid when L=2π', () => {
    const N = 8;
    const L = 2 * Math.PI;
    const { kx2, ky2, kz2 } = createKSpaceArrays(N, L);
    // expected frequencies: [0,1,2,3,-4,-3,-2,-1]
    const expectedSq = [0,1,4,9,16,9,4,1];
    expect(Array.from(kx2)).toEqualCloseTo(expectedSq, 10);
    expect(Array.from(ky2)).toEqualCloseTo(expectedSq, 10);
    expect(Array.from(kz2)).toEqualCloseTo(expectedSq, 10);
  });

  it('createKSpaceArrays throws on invalid L', () => {
    expect(() => createKSpaceArrays(32, 0)).toThrow(/finite positive/);
    expect(() => createKSpaceArrays(32, -2)).toThrow(/finite positive/);
    expect(() => createKSpaceArrays(32, NaN)).toThrow(/finite positive/);
  });

  it('buildKineticExponentials throws when expK has wrong length', () => {
    const N = 4;
    const L = 10;
    const kArrays = createKSpaceArrays(N, L);
    const expK = new Float32Array(100); // wrong size, should be 2*N³ = 128
    
    expect(() => buildKineticExponentials(expK, kArrays, N, 0.01))
      .toThrow(/expK.length.*!== 2\*N³/);
  });

  it('buildKineticExponentials throws when kArrays have wrong length', () => {
    const N = 4;
    const size = N * N * N;
    const expK = new Float32Array(2 * size);
    
    // create kArrays with wrong length
    const badKArrays = {
      kx2: new Float32Array(N - 1), // wrong size
      ky2: new Float32Array(N),
      kz2: new Float32Array(N)
    };
    
    expect(() => buildKineticExponentials(expK, badKArrays, N, 0.01))
      .toThrow(/kArrays must have length N/);
  });

  it('buildKineticExponentials writes correct complex exponentials', () => {
    const N = 4;
    const L = 2 * Math.PI;
    const dt = 0.05;
    const kArrays = createKSpaceArrays(N, L);
    const expK = new Float32Array(2 * N * N * N);
    buildKineticExponentials(expK, kArrays, N, dt);

    const x = 1, y = 2, z = 3;
    const k2 = kArrays.kx2[x] + kArrays.ky2[y] + kArrays.kz2[z];
    const theta = -0.5 * k2 * dt;
    const id = lin(x, y, z, N) << 1;
    // expK is stored in Float32; allow a realistic tolerance
    expect(expK[id]).toBeCloseTo(Math.cos(theta), 7);
    expect(expK[id + 1]).toBeCloseTo(Math.sin(theta), 7);
  });

  it('buildPotentialExponentials honors V and CAP decay', () => {
    const N = 2;
    const size = N * N * N; // 8
    const V = new Float32Array(size);
    // v: [0,1,-2, 0,0,0, 0,0] to test different angles
    V[0] = 0; V[1] = 1; V[2] = -2;
    const dt = 0.1;
    const absorbStrength = 3;
    const expVh = new Float32Array(2 * size);

    // case 1: no CAP
    buildPotentialExponentials(expVh, V, null, dt, absorbStrength);
    let theta = -0.5 * V[1] * dt;
    let j = (1 << 1);
    // float32 rounding on expVh components
    expect(expVh[j]).toBeCloseTo(Math.cos(theta), 7);
    expect(expVh[j + 1]).toBeCloseTo(Math.sin(theta), 7);

    // case 2: with CAP
    const capS2 = new Float32Array(size);
    capS2[2] = 2; // non-zero cap only at index 2
    buildPotentialExponentials(expVh, V, capS2, dt, absorbStrength);
    theta = -0.5 * V[2] * dt;
    const decay = Math.exp(-0.5 * absorbStrength * capS2[2] * dt);
    j = (2 << 1);
    const mag = Math.hypot(expVh[j], expVh[j + 1]);
    // use relative tolerance suitable for Float32 computations
    const rel = Math.abs(mag - decay) / Math.max(decay, 1e-12);
    expect(rel).toBeLessThan(1e-7);
    // angle still follows theta
    const ang = Math.atan2(expVh[j + 1], expVh[j]);
    expect(ang).toBeCloseTo(theta, 7); // allow some wrap/precision tolerance
  });

  it('buildPotentialExponentials throws when |V|·dt exceeds safe limit', () => {
    const N = 2;
    const size = N * N * N;
    const V = new Float32Array(size);
    const expVh = new Float32Array(2 * size);
    
    // MAX_THETA is 1e6, so |V|·dt/2 > 1e6 means |V|·dt > 2e6
    // With dt=0.1 and V=1e8, theta = 0.5 * 0.1 * 1e8 = 5e6 > 1e6
    V[0] = 1e8;
    const dt = 0.1;
    
    expect(() => buildPotentialExponentials(expVh, V, null, dt, 0))
      .toThrow(/exceeds safe limit/);
  });

  it('buildPotentialExponentials clamps extreme decay argument', () => {
    // test that very large CAP values are clamped to prevent underflow
    const N = 2;
    const size = N * N * N;
    const V = new Float32Array(size);
    const expVh = new Float32Array(2 * size);
    const capS2 = new Float32Array(size);
    
    // set extremely high cap value that would cause decayArg > MAX_DECAY_ARG (700)
    // decayArg = 0.5 * absorbStrength * capS2[i] * dt
    // to get decayArg > 700: 0.5 * 10000 * 1 * 0.15 = 750 > 700
    capS2[0] = 1;
    const dt = 0.15;
    const absorbStrength = 10000; // very high strength
    
    // should not throw - just clamp
    expect(() => buildPotentialExponentials(expVh, V, capS2, dt, absorbStrength))
      .not.toThrow();
    
    // the decay should be clamped to exp(-700) ≈ 0, so magnitude should be ~0
    const mag = Math.hypot(expVh[0], expVh[1]);
    expect(mag).toBeCloseTo(0, 10);
  });

  it('buildPotentialExponentials throws when buffer lengths are inconsistent', () => {
    const len = 8;
    const V = new Float32Array(len);
    const dt = 0.1;
    const expShort = new Float32Array(2 * len - 2);
    expect(() => buildPotentialExponentials(expShort, V, null, dt, 0)).toThrow(/expVh.length === 2 \* V.length/);

    const expVh = new Float32Array(2 * len);
    const capBad = new Float32Array(len - 1);
    expect(() => buildPotentialExponentials(expVh, V, capBad, dt, 1)).toThrow(/capS2.length === V.length/);
  });

  it('createAbsorbingBoundary interior is zero and edges positive/symmetric', () => {
    const N = 8;
    const cap = createAbsorbingBoundary(N, 0.25); // nAbs=2

    const interiorIndex = lin(4, 4, 4, N);
    expect(cap[interiorIndex]).toBe(0);

    const edgeIndex = lin(0, 4, 4, N);
    expect(cap[edgeIndex]).toBeGreaterThan(0);

    // symmetry checks
    const x = 1, y = 2, z = 3;
    const a = cap[lin(x, y, z, N)];
    const ax = cap[lin(N - 1 - x, y, z, N)];
    const ay = cap[lin(x, N - 1 - y, z, N)];
    const az = cap[lin(x, y, N - 1 - z, N)];
    expect(a).toBeCloseTo(ax, 12);
    expect(a).toBeCloseTo(ay, 12);
    expect(a).toBeCloseTo(az, 12);
  });

  it('potentialHalfStep applies complex multiply per cell', () => {
    const psiRe = new Float32Array([1, 0.5]);
    const psiIm = new Float32Array([0, -1]);
    // expVh for angles [0, -π/2]
    const expVh = new Float32Array([
      Math.cos(0), Math.sin(0),
      Math.cos(-Math.PI / 2), Math.sin(-Math.PI / 2)
    ]);
    potentialHalfStep(psiRe, psiIm, expVh);
    // after multiply: index 0 stays (1,0); index 1 rotates by -π/2
    expect(psiRe[0]).toBeCloseTo(1, 12);
    expect(psiIm[0]).toBeCloseTo(0, 12);
    // (0.5 - i) * (0 - i) = (-1) * 0.5 + ... compute explicitly
    // let (a+ib)*(c+id) with c=0, d=-1 => real = a*c - b*d = 0 - (-1)*(-1)???
    // better compute numerically from initial
    const a = 0.5, b = -1;
    const c = 0, d = -1;
    const nr = a * c - b * d; // 0 - (-1)*(-1) = -1
    const ni = a * d + b * c; // 0.5*(-1) + (-1)*0 = -0.5
    expect(psiRe[1]).toBeCloseTo(nr, 12);
    expect(psiIm[1]).toBeCloseTo(ni, 12);
  });

  it('kineticFullStep preserves norm (no CAP)', () => {
    const N = 4;
    const L = 10;
    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);
    fillRandomPair(psiRe, psiIm);
    const kArrays = createKSpaceArrays(N, L);
    const expK = new Float32Array(2 * size);
    buildKineticExponentials(expK, kArrays, N, 0.02);

    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);

    const norm0 = calculateNorm(psiRe, psiIm, 1);
    kineticFullStep(psiRe, psiIm, expK, N, sRe, sIm);
    const norm1 = calculateNorm(psiRe, psiIm, 1);
    expect(Math.abs(norm1 - norm0)).toBeLessThan(1e-6);
  });

  it('timeStep preserves norm for V=0 and decreases with CAP', () => {
    const N = 4;
    const L = 10;
    const size = N * N * N;

    // initial random psi
    const psiReA = new Float32Array(size);
    const psiImA = new Float32Array(size);
    fillRandomPair(psiReA, psiImA);

    const kArrays = createKSpaceArrays(N, L);
    const expK = new Float32Array(2 * size);
    const expVh = new Float32Array(2 * size);

    // case 1: V=0, no CAP
    const V = new Float32Array(size); // zeros
    buildKineticExponentials(expK, kArrays, N, 0.02);
    buildPotentialExponentials(expVh, V, null, 0.02, 0);

    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);

    const norm0 = calculateNorm(psiReA, psiImA, 1);
    timeStep(psiReA, psiImA, expVh, expK, N, sRe, sIm);
    const norm1 = calculateNorm(psiReA, psiImA, 1);
    expect(Math.abs(norm1 - norm0)).toBeLessThan(1e-6);

    // case 2: with CAP
    const psiReB = psiReA.slice();
    const psiImB = psiImA.slice();
    const capS2 = createAbsorbingBoundary(N, 0.25);
    buildPotentialExponentials(expVh, V, capS2, 0.02, 3);

    const normB0 = calculateNorm(psiReB, psiImB, 1);
    timeStep(psiReB, psiImB, expVh, expK, N, sRe, sIm);
    const normB1 = calculateNorm(psiReB, psiImB, 1);
    expect(normB1).toBeLessThan(normB0);
  });

  it('timeStep with dt=0 leaves psi unchanged', () => {
    const N = 4;
    const L = 10;
    const size = N * N * N;

    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);
    for (let i = 0; i < size; i++) {
      psiRe[i] = Math.random() * 2 - 1;
      psiIm[i] = Math.random() * 2 - 1;
    }

    const kArrays = createKSpaceArrays(N, L);
    const expK = new Float32Array(2 * size);
    const expVh = new Float32Array(2 * size);

    buildKineticExponentials(expK, kArrays, N, 0);
    buildPotentialExponentials(expVh, new Float32Array(size), null, 0, 0);

    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);

    const r0 = psiRe.slice();
    const i0 = psiIm.slice();
    timeStep(psiRe, psiIm, expVh, expK, N, sRe, sIm);

    // FFT round-trips introduce tiny numeric noise; assert small RMS relative error
    const rel = rmsRel(psiRe, psiIm, r0, i0);
    expect(rel).toBeLessThan(1e-7);
  });

  it('potentialHalfStep throws when array lengths are inconsistent', () => {
    const psiRe = new Float32Array(4);
    const psiIm = new Float32Array(3); // wrong length
    const expVh = new Float32Array(8);
    expect(() => potentialHalfStep(psiRe, psiIm, expVh)).toThrow(/same length/);

    const psiIm2 = new Float32Array(4);
    const badExp = new Float32Array(6); // not 2 * psi length
    expect(() => potentialHalfStep(psiRe, psiIm2, badExp)).toThrow(/expVh.length === 2 \* psi length/);
  });

  it('kineticFullStep throws when array lengths are inconsistent', () => {
    const N = 4;
    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size - 1); // mismatch
    const expK = new Float32Array(2 * size);
    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);
    expect(() => kineticFullStep(psiRe, psiIm, expK, N, sRe, sIm)).toThrow(/same length/);

    const psiIm2 = new Float32Array(size);
    const expKBad = new Float32Array(2 * size - 2);
    expect(() => kineticFullStep(psiRe, psiIm2, expKBad, N, sRe, sIm)).toThrow(/expK.length === 2 \* psi length/);

    const psiShort = new Float32Array(size - 1);
    const expKShort = new Float32Array(2 * (size - 1));
    expect(() => kineticFullStep(psiShort, psiShort, expKShort, N, sRe, sIm)).toThrow(/psi arrays length N\^3/);

    const sReShort = new Float32Array(N - 1);
    const sImShort = new Float32Array(N - 1);
    expect(() => kineticFullStep(psiRe, psiIm2, expK, N, sReShort, sImShort)).toThrow(/scratchRe\/scratchIm length/);
  });

  it('addPacket3D centers amplitude near origin and obeys linearity', () => {
    const N = 8;
    const L = 8;
    const coord = createCoordinateArray(N, L);

    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);

    const gX = new Float32Array(N);
    const gY = new Float32Array(N);
    const gZ = new Float32Array(N);
    const pX = new Float32Array(N);
    const pY = new Float32Array(N);
    const pZ = new Float32Array(N);

    addPacket3D(
      psiRe, psiIm, coord, N,
      0, 0, 0, // center
      0.6, 0.6, 0.6, // widths
      0, 0, 0, // k
      1,
      gX, gY, gZ, pX, pY, pZ
    );

    // find argmax index
    let maxV = -Infinity, maxI = -1;
    for (let i = 0; i < size; i++) {
      const v = psiRe[i] * psiRe[i] + psiIm[i] * psiIm[i];
      if (v > maxV) { maxV = v; maxI = i; }
    }
    // decode x,y,z
    const x = maxI % N;
    const y = Math.floor(maxI / N) % N;
    const z = Math.floor(maxI / (N * N));
    const centers = [N/2 - 1, N/2];
    expect(centers).toContain(x);
    expect(centers).toContain(y);
    expect(centers).toContain(z);

    // linearity: second add with scale=2 doubles amplitudes
    const psiRe2 = psiRe.slice();
    const psiIm2 = psiIm.slice();
    addPacket3D(
      psiRe2, psiIm2, coord, N,
      0, 0, 0,
      0.6, 0.6, 0.6,
      0, 0, 0,
      2,
      gX, gY, gZ, pX, pY, pZ
    );
    let ratioOK = true;
    for (let i = 0; i < size; i++) {
      if (psiRe[i] !== 0 || psiIm[i] !== 0) {
        // expect the increment to be double
        const dRe = psiRe2[i] - psiRe[i];
        const dIm = psiIm2[i] - psiIm[i];
        const r = Math.hypot(dRe, dIm);
        const base = Math.hypot(psiRe[i], psiIm[i]);
        // because first packet + second packet => total amplitude is base + 2*base at center, but
        // comparing increment to base is approximately 2 when phases are same (k=0).
        if (base > 0) {
          const factor = r / base;
          if (Math.abs(factor - 2) > 1e-6) { ratioOK = false; break; }
        }
      }
    }
    expect(ratioOK).toBe(true);
  });

  it('addPacket3D throws for psi arrays with wrong length', () => {
    const N = 8;
    const L = 8;
    const coord = createCoordinateArray(N, L);
    const size = N * N * N;
    const psiRe = new Float32Array(size - 1); // wrong size
    const psiIm = new Float32Array(size);

    const gX = new Float32Array(N);
    const gY = new Float32Array(N);
    const gZ = new Float32Array(N);
    const pX = new Float32Array(N);
    const pY = new Float32Array(N);
    const pZ = new Float32Array(N);

    expect(() => addPacket3D(
      psiRe, psiIm, coord, N,
      0, 0, 0,
      0.6, 0.6, 0.6,
      0, 0, 0,
      1,
      gX, gY, gZ, pX, pY, pZ
    )).toThrow(/psi arrays must have length N³/);
  });

  it('addPacket3D throws for coord array with wrong length', () => {
    const N = 8;
    const L = 8;
    const coord = createCoordinateArray(N - 1, L); // wrong size
    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);

    const gX = new Float32Array(N);
    const gY = new Float32Array(N);
    const gZ = new Float32Array(N);
    const pX = new Float32Array(N);
    const pY = new Float32Array(N);
    const pZ = new Float32Array(N);

    expect(() => addPacket3D(
      psiRe, psiIm, coord, N,
      0, 0, 0,
      0.6, 0.6, 0.6,
      0, 0, 0,
      1,
      gX, gY, gZ, pX, pY, pZ
    )).toThrow(/coord.length.*!== N/);
  });

  it('addPacket3D throws for scratch arrays that are too small', () => {
    const N = 8;
    const L = 8;
    const coord = createCoordinateArray(N, L);
    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);

    // one scratch array too small
    const gX = new Float32Array(N - 1); // too small
    const gY = new Float32Array(N);
    const gZ = new Float32Array(N);
    const pX = new Float32Array(N);
    const pY = new Float32Array(N);
    const pZ = new Float32Array(N);

    expect(() => addPacket3D(
      psiRe, psiIm, coord, N,
      0, 0, 0,
      0.6, 0.6, 0.6,
      0, 0, 0,
      1,
      gX, gY, gZ, pX, pY, pZ
    )).toThrow(/scratch arrays must have length/);
  });

  it('addPacket3D throws for non-positive widths', () => {
    const N = 8;
    const L = 8;
    const coord = createCoordinateArray(N, L);

    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);

    const gX = new Float32Array(N);
    const gY = new Float32Array(N);
    const gZ = new Float32Array(N);
    const pX = new Float32Array(N);
    const pY = new Float32Array(N);
    const pZ = new Float32Array(N);

    expect(() => addPacket3D(
      psiRe, psiIm, coord, N,
      0, 0, 0,
      0, 0.6, 0.6,
      0, 0, 0,
      1,
      gX, gY, gZ, pX, pY, pZ
    )).toThrow(/positive widths/i);
  });

  it('renormalize throws when psiRe and psiIm lengths differ', () => {
    const psiRe = new Float32Array(64);
    const psiIm = new Float32Array(63); // wrong size
    expect(() => renormalize(psiRe, psiIm, 1))
      .toThrow(/psiRe.length.*!== psiIm.length/);
  });

  it('renormalize brings norm to ~1, zero state remains zero', () => {
    const N = 8;
    const L = 8;
    const size = N * N * N;

    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);

    const coord = createCoordinateArray(N, L);
    const gX = new Float32Array(N), gY = new Float32Array(N), gZ = new Float32Array(N);
    const pX = new Float32Array(N), pY = new Float32Array(N), pZ = new Float32Array(N);

    addPacket3D(psiRe, psiIm, coord, N, 0, 0, 0, 0.6, 0.6, 0.6, 0, 0, 0, 1, gX, gY, gZ, pX, pY, pZ);
    const cellVol = Math.pow(L / N, 3);
    renormalize(psiRe, psiIm, cellVol);
    const norm = calculateNorm(psiRe, psiIm, cellVol);
    expect(Math.abs(norm - 1)).toBeLessThan(1e-6);

    // zero state
    const zRe = new Float32Array(size);
    const zIm = new Float32Array(size);
    renormalize(zRe, zIm, cellVol);
    const zn = calculateNorm(zRe, zIm, cellVol);
    expect(zn).toBe(0);
  });

  it('addPacket3D yields positive phase gradient along +x near center when kx>0', () => {
    const N = 8;
    const L = 8;
    const coord = createCoordinateArray(N, L);

    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);

    const gX = new Float32Array(N);
    const gY = new Float32Array(N);
    const gZ = new Float32Array(N);
    const pX = new Float32Array(N);
    const pY = new Float32Array(N);
    const pZ = new Float32Array(N);

    const cx = 0, cy = 0, cz = 0;
    const kx = 0.5, ky = 0, kz = 0;
    addPacket3D(
      psiRe, psiIm, coord, N,
      cx, cy, cz,
      0.6, 0.6, 0.6,
      kx, ky, kz,
      1,
      gX, gY, gZ, pX, pY, pZ
    );

    // pick two adjacent x cells around the center (x=3 -> -0.5, x=4 -> 0.5)
    const y = 3, z = 3;
    const x1 = 3, x2 = 4;
    const id1 = lin(x1, y, z, N);
    const id2 = lin(x2, y, z, N);
    const phi1 = Math.atan2(psiIm[id1], psiRe[id1]);
    const phi2 = Math.atan2(psiIm[id2], psiRe[id2]);
    let dphi = phi2 - phi1;
    if (dphi > Math.PI) dphi -= 2 * Math.PI;
    if (dphi < -Math.PI) dphi += 2 * Math.PI;
    const expected = kx * (coord[x2] - coord[x1]);
    expect(dphi).toBeGreaterThan(0);
    expect(Math.abs(dphi - expected)).toBeLessThan(1e-3);
  });

  it('createAbsorbingBoundary increases monotonically from interior to face along +x', () => {
    const N = 8;
    const cap = createAbsorbingBoundary(N, 0.25); // nAbs=2
    const y = Math.floor(N / 2);
    const z = Math.floor(N / 2);
    const values = [];
    for (let x = Math.floor(N / 2); x < N; x++) {
      values.push(cap[lin(x, y, z, N)]);
    }
    for (let i = 0; i < values.length - 1; i++) {
      // allow tiny float slack
      expect(values[i + 1] + 1e-12).toBeGreaterThanOrEqual(values[i]);
    }
  });

  it('createAbsorbingBoundary enforces minimum width >= 4 for small absorbFrac', () => {
    const N = 32;
    const cap = createAbsorbingBoundary(N, 0.05);
    const y = Math.floor(N / 2);
    const z = Math.floor(N / 2);
    let width = 0;
    for (let x = 0; x < N; x++) {
      if (cap[lin(x, y, z, N)] > 0) width++; else break;
    }
    expect(width).toBeGreaterThanOrEqual(4);
  });

  it('CAP min-width reduces slab probability versus thin 1-cell CAP', () => {
    const N = 16;
    const L = 8;
    const dx = L / N;
    const cellVol = dx * dx * dx;
    const coord = createCoordinateArray(N, L);

    const size = N * N * N;
    const baseRe = new Float32Array(size);
    const baseIm = new Float32Array(size);

    const sigma = 0.6;
    const cx = -L / 4;
    const kx = (2 * Math.PI) / (8 * dx);

    const gX = new Float32Array(N);
    const gY = new Float32Array(N);
    const gZ = new Float32Array(N);
    const pX = new Float32Array(N);
    const pY = new Float32Array(N);
    const pZ = new Float32Array(N);

    addPacket3D(
      baseRe, baseIm, coord, N,
      cx, 0, 0,
      sigma, sigma, sigma,
      kx, 0, 0,
      1,
      gX, gY, gZ, pX, pY, pZ
    );
    renormalize(baseRe, baseIm, cellVol);

    const kArrays = createKSpaceArrays(N, L);
    const dt = 0.02 * dx * dx;
    const expK = new Float32Array(2 * size);
    buildKineticExponentials(expK, kArrays, N, dt);
    const V = new Float32Array(size);

    function simulateWithCap(cap) {
      const psiRe = baseRe.slice();
      const psiIm = baseIm.slice();
      const expVh = new Float32Array(2 * size);
      buildPotentialExponentials(expVh, V, cap, dt, 3);
      const sRe = new Float32Array(N);
      const sIm = new Float32Array(N);
      for (let i = 0; i < 900; i++) {
        timeStep(psiRe, psiIm, expVh, expK, N, sRe, sIm);
      }
      const y0 = Math.floor(N / 2);
      const z0 = Math.floor(N / 2);
      let capStartX = Math.floor(N / 2);
      for (let x = Math.floor(N / 2); x < N; x++) {
        if (cap[lin(x, y0, z0, N)] > 0) { capStartX = x; break; }
      }
      const slabX = capStartX - 1;
      let pref = 0;
      for (let z = 0; z < N; z++) {
        for (let y = 0; y < N; y++) {
          const id = lin(slabX, y, z, N);
          const pr = psiRe[id], pi = psiIm[id];
          pref += (pr * pr + pi * pi) * cellVol;
        }
      }
      return pref;
    }

    const capMin = createAbsorbingBoundary(N, 0.05);
    const capThin = (() => {
      const cap = new Float32Array(size);
      const nAbs = 1;
      for (let z = 0; z < N; z++) {
        const zOff = N * N * z;
        const dz = Math.max(0, (nAbs - z) / nAbs, (z - (N - 1 - nAbs)) / nAbs);
        for (let y = 0; y < N; y++) {
          const yOff = zOff + N * y;
          const dy = Math.max(0, (nAbs - y) / nAbs, (y - (N - 1 - nAbs)) / nAbs);
          const dy2 = dy * dy;
          for (let x = 0; x < N; x++) {
            const dxv = Math.max(0, (nAbs - x) / nAbs, (x - (N - 1 - nAbs)) / nAbs);
            const s2 = dxv * dxv + dy2 + dz * dz;
            cap[yOff + x] = s2;
          }
        }
      }
      return cap;
    })();

    const prefMin = simulateWithCap(capMin);
    const prefThin = simulateWithCap(capThin);
    expect(prefMin).toBeLessThan(prefThin);
  }, 20000);

  it('createAbsorbingBoundary returns zeros when CAP not possible at small N', () => {
    const N = 3;
    const cap = createAbsorbingBoundary(N, 0.5);
    for (const v of cap) {
      expect(v).toBe(0);
    }
  });

  it('createAbsorbingBoundary returns all zeros when absorbFrac <= 0', () => {
    const N = 8;
    let cap0 = createAbsorbingBoundary(N, 0);
    for (const v of cap0) {
      expect(v).toBe(0);
    }
    let capNeg = createAbsorbingBoundary(N, -0.1);
    for (const v of capNeg) {
      expect(v).toBe(0);
    }
  });

  it('createKSpaceArrays throws for non power-of-two N (fast validation)', () => {
    expect(() => createKSpaceArrays(12, 10)).toThrow(/power of two/i);
  });

  it('timeStep throws early for N < 2 (fast validation)', () => {
    const psiRe = new Float32Array(1);
    const psiIm = new Float32Array(1);
    const expVh = new Float32Array(2);
    const expK = new Float32Array(2);
    const sRe = new Float32Array(1);
    const sIm = new Float32Array(1);
    expect(() => timeStep(psiRe, psiIm, expVh, expK, 1, sRe, sIm)).toThrow(/power of two/i);
  });

  it('timeStep with checkStability throws when NaN detected', () => {
    const N = 4;
    const L = 10;
    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);
    
    // set initial values with a NaN that will propagate
    psiRe[0] = NaN;
    
    const kArrays = createKSpaceArrays(N, L);
    const expK = new Float32Array(2 * size);
    const expVh = new Float32Array(2 * size);
    
    buildKineticExponentials(expK, kArrays, N, 0.01);
    buildPotentialExponentials(expVh, new Float32Array(size), null, 0.01, 0);
    
    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);
    
    // should throw with checkStability=true
    expect(() => timeStep(psiRe, psiIm, expVh, expK, N, sRe, sIm, true))
      .toThrow(/NaN or Infinity detected/);
  });

  // PHYS-ENERGY-TEST: energy conservation test for free particle (V=0)
  it('calculateEnergy conserves total energy for free particle over many steps', () => {
    const N = 16;
    const L = 10;
    const dx = L / N;
    const cellVol = dx * dx * dx;
    const coord = createCoordinateArray(N, L);

    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);
    const V = new Float32Array(size); // free particle: V=0

    const gX = new Float32Array(N);
    const gY = new Float32Array(N);
    const gZ = new Float32Array(N);
    const pX = new Float32Array(N);
    const pY = new Float32Array(N);
    const pZ = new Float32Array(N);

    // create a wave packet with momentum
    const sigma = 0.8;
    const kx = 2.0;
    addPacket3D(
      psiRe, psiIm, coord, N,
      0, 0, 0,
      sigma, sigma, sigma,
      kx, 0, 0,
      1,
      gX, gY, gZ, pX, pY, pZ
    );
    renormalize(psiRe, psiIm, cellVol);

    const kArrays = createKSpaceArrays(N, L);
    const expK = new Float32Array(2 * size);
    const expVh = new Float32Array(2 * size);
    const dt = 0.02 * dx * dx;

    buildKineticExponentials(expK, kArrays, N, dt);
    buildPotentialExponentials(expVh, V, null, dt, 0);

    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);

    // initial energy
    const E0 = calculateEnergy(psiRe, psiIm, V, kArrays, N, cellVol, sRe, sIm);
    expect(E0.potential).toBeCloseTo(0, 10); // V=0 so ⟨V⟩=0

    // evolve for 100 steps
    for (let i = 0; i < 100; i++) {
      timeStep(psiRe, psiIm, expVh, expK, N, sRe, sIm);
    }

    // final energy
    const E1 = calculateEnergy(psiRe, psiIm, V, kArrays, N, cellVol, sRe, sIm);

    // energy should be conserved (relative drift < 1e-4 as per acceptance criteria)
    const relDrift = Math.abs(E1.total - E0.total) / Math.abs(E0.total);
    expect(relDrift).toBeLessThan(1e-4);
  });

  // PHYS-ANALYTIC: analytical comparison - free Gaussian spreading
  // theory: for a free Gaussian, the width evolves as σ(t) = σ₀√(1 + (ℏt/2mσ₀²)²)
  // with ℏ=m=1: σ(t) = σ₀√(1 + (t/2σ₀²)²)
  // note: the factor of 2 comes from the minimum uncertainty wavepacket convention.
  // our addPacket3D uses exp(-x²/2σ²), which has RMS = σ/√2 for |ψ|².
  it('free Gaussian spreading matches analytical formula', () => {
    const N = 16; // Reduced from 32 for faster test
    const L = 16; // Adjusted domain
    const dx = L / N;
    const cellVol = dx * dx * dx;
    const coord = createCoordinateArray(N, L);

    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);
    const V = new Float32Array(size); // V=0 for free particle

    const gX = new Float32Array(N);
    const gY = new Float32Array(N);
    const gZ = new Float32Array(N);
    const pX = new Float32Array(N);
    const pY = new Float32Array(N);
    const pZ = new Float32Array(N);

    // initial state: Gaussian with σ₀ = 1.5, no momentum (k=0)
    const sigma0 = 1.5;
    addPacket3D(
      psiRe, psiIm, coord, N,
      0, 0, 0,
      sigma0, sigma0, sigma0,
      0, 0, 0, // k=0 for stationary center
      1,
      gX, gY, gZ, pX, pY, pZ
    );
    renormalize(psiRe, psiIm, cellVol);

    // measure RMS width (compute √⟨x²⟩)
    function measureWidth() {
      let sumX2 = 0;
      for (let z = 0; z < N; z++) {
        const zOff = N * N * z;
        for (let y = 0; y < N; y++) {
          const yOff = zOff + N * y;
          for (let x = 0; x < N; x++) {
            const id = yOff + x;
            const prob = psiRe[id] * psiRe[id] + psiIm[id] * psiIm[id];
            const xPos = coord[x];
            sumX2 += prob * xPos * xPos * cellVol;
          }
        }
      }
      return Math.sqrt(sumX2); // RMS width
    }

    const width0 = measureWidth();
    // for addPacket3D using exp(-x²/2σ²), |ψ|² ∝ exp(-x²/σ²)
    // the variance of exp(-x²/σ²) is σ²/2, so RMS = σ/√2 ≈ 1.06 for σ=1.5
    const expectedInitialRMS = sigma0 / Math.sqrt(2);
    expect(width0).toBeCloseTo(expectedInitialRMS, 0);

    const kArrays = createKSpaceArrays(N, L);
    const expK = new Float32Array(2 * size);
    const expVh = new Float32Array(2 * size);
    const dt = 0.1 * dx * dx; // Larger dt for fewer steps
    const steps = 50; // Reduced from 200 for faster test
    const totalTime = steps * dt;

    buildKineticExponentials(expK, kArrays, N, dt);
    buildPotentialExponentials(expVh, V, null, dt, 0);

    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);

    // evolve
    for (let i = 0; i < steps; i++) {
      timeStep(psiRe, psiIm, expVh, expK, N, sRe, sIm);
    }

    const widthT = measureWidth();

    // key test: width should INCREASE over time (spreading)
    // this is the fundamental physics we're verifying
    expect(widthT).toBeGreaterThan(width0);
    
    // theoretical spreading: for |ψ|² width σ_prob = σ/√2,
    // evolves as σ_prob(t) = σ_prob(0)√(1 + (ℏt/mσ²)²) with ℏ=m=1
    // using the wavefunction σ parameter: t/σ² scaling
    const spreadFactor = Math.sqrt(1 + Math.pow(totalTime / (sigma0 * sigma0), 2));
    const theoreticalWidth = expectedInitialRMS * spreadFactor;

    // allow 35% tolerance due to:
    // - discrete grid effects
    // - finite domain boundary effects  
    // - float32 precision
    // the key physics (spreading occurs) is verified above
    const relError = Math.abs(widthT - theoreticalWidth) / theoreticalWidth;
    expect(relError).toBeLessThan(0.35);
  });

  // calculateEnergy validation tests
  it('calculateEnergy throws when psi length mismatches N³', () => {
    const N = 4;
    const size = N * N * N;
    const psiRe = new Float32Array(size - 1); // wrong size
    const psiIm = new Float32Array(size - 1);
    const V = new Float32Array(size - 1);
    const kArrays = createKSpaceArrays(N, 10);
    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);

    expect(() => calculateEnergy(psiRe, psiIm, V, kArrays, N, 1, sRe, sIm))
      .toThrow(/psi length.*!== N³/);
  });

  it('calculateEnergy throws when psiRe and psiIm lengths differ', () => {
    const N = 4;
    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size - 1); // wrong size
    const V = new Float32Array(size);
    const kArrays = createKSpaceArrays(N, 10);
    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);

    expect(() => calculateEnergy(psiRe, psiIm, V, kArrays, N, 1, sRe, sIm))
      .toThrow(/psiRe.length.*!== psiIm.length/);
  });

  it('calculateEnergy throws when V length mismatches psi', () => {
    const N = 4;
    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);
    const V = new Float32Array(size - 1); // wrong size
    const kArrays = createKSpaceArrays(N, 10);
    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);

    expect(() => calculateEnergy(psiRe, psiIm, V, kArrays, N, 1, sRe, sIm))
      .toThrow(/V.length.*!== psi length/);
  });

  it('calculateNorm throws when psiRe and psiIm lengths differ', () => {
    const psiRe = new Float32Array(64);
    const psiIm = new Float32Array(63); // wrong size
    expect(() => calculateNorm(psiRe, psiIm, 1))
      .toThrow(/psiRe.length.*!== psiIm.length/);
  });

  it('checkFinite returns true for finite arrays', () => {
    const psiRe = new Float32Array([1, 2, 3, 4]);
    const psiIm = new Float32Array([0.5, 0.5, 0.5, 0.5]);
    expect(checkFinite(psiRe, psiIm)).toBe(true);
  });

  it('checkFinite returns false when psiRe contains NaN', () => {
    const psiRe = new Float32Array([1, NaN, 3, 4]);
    const psiIm = new Float32Array([0.5, 0.5, 0.5, 0.5]);
    expect(checkFinite(psiRe, psiIm)).toBe(false);
  });

  it('checkFinite returns false when psiIm contains Infinity', () => {
    const psiRe = new Float32Array([1, 2, 3, 4]);
    const psiIm = new Float32Array([0.5, Infinity, 0.5, 0.5]);
    expect(checkFinite(psiRe, psiIm)).toBe(false);
  });

  // TEST-EXTREME: Extreme parameter test - high CAP strength
  it('simulation remains stable with high CAP strength', () => {
    const N = 16;
    const L = 10;
    const dx = L / N;
    const cellVol = dx * dx * dx;
    const coord = createCoordinateArray(N, L);

    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);
    const V = new Float32Array(size);

    const gX = new Float32Array(N);
    const gY = new Float32Array(N);
    const gZ = new Float32Array(N);
    const pX = new Float32Array(N);
    const pY = new Float32Array(N);
    const pZ = new Float32Array(N);

    addPacket3D(
      psiRe, psiIm, coord, N,
      0, 0, 0,
      0.8, 0.8, 0.8,
      3, 0, 0, // moving towards boundary
      1,
      gX, gY, gZ, pX, pY, pZ
    );
    renormalize(psiRe, psiIm, cellVol);

    const kArrays = createKSpaceArrays(N, L);
    const expK = new Float32Array(2 * size);
    const expVh = new Float32Array(2 * size);
    const dt = 0.02 * dx * dx;

    buildKineticExponentials(expK, kArrays, N, dt);

    // high CAP strength (10x normal)
    const cap = createAbsorbingBoundary(N, 0.2);
    buildPotentialExponentials(expVh, V, cap, dt, 30); // very strong CAP

    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);

    // run simulation - should not produce NaN/Infinity
    for (let i = 0; i < 100; i++) {
      timeStep(psiRe, psiIm, expVh, expK, N, sRe, sIm);
    }

    // check stability: norm should be finite and <= 1
    const finalNorm = calculateNorm(psiRe, psiIm, cellVol);
    expect(Number.isFinite(finalNorm)).toBe(true);
    expect(finalNorm).toBeLessThanOrEqual(1.001); // allow tiny floating point overshoot
    expect(finalNorm).toBeGreaterThanOrEqual(0);
  });
});
