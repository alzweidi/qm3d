import { describe, it, expect } from 'vitest';
import { lin, fillRandomPair, rmsRelativeErrorComplex as rmsRel, allocPacketScratch } from '../utils/testUtils.js';
import { presetHarmonic } from '../../src/physics/potentials.js';
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
  calculateEnergy
} from '../../src/physics/quantum.js';
import { fft3d } from '../../src/physics/fft.js';

// lin moved to tests/utils/testUtils.js

describe('quantum.js', () => {
  function meanX(psiRe, psiIm, coord, N, cellVol) {
    let sum = 0;
    for (let z = 0; z < N; z++) {
      for (let y = 0; y < N; y++) {
        for (let x = 0; x < N; x++) {
          const id = lin(x, y, z, N);
          const pr = psiRe[id], pi = psiIm[id];
          sum += coord[x] * (pr * pr + pi * pi);
        }
      }
    }
    return sum * cellVol;
  }
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

  it('calculateEnergy is zero for the all-zero state (V null)', () => {
    const N = 8;
    const L = 8;
    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);
    const E = calculateEnergy(psiRe, psiIm, null, N, L);
    expect(E).toBeCloseTo(0, 12);
  });

  it('calculateEnergy matches 0.5*k^2 for a free plane wave', () => {
    const N = 32;
    const L = 10;
    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);
    const V = new Float32Array(size);

    const coord = createCoordinateArray(N, L);
    const m = 3;
    const kx = (2 * Math.PI / L) * m;
    const amp = 1 / Math.sqrt(L * L * L);

    for (let z = 0; z < N; z++) {
      const zOff = N * N * z;
      for (let y = 0; y < N; y++) {
        const yOff = zOff + N * y;
        for (let x = 0; x < N; x++) {
          const id = yOff + x;
          const phase = kx * coord[x];
          psiRe[id] = amp * Math.cos(phase);
          psiIm[id] = amp * Math.sin(phase);
        }
      }
    }

    const dx = L / N;
    const cellVol = dx * dx * dx;
    expect(calculateNorm(psiRe, psiIm, cellVol)).toBeCloseTo(1, 4);

    const E = calculateEnergy(psiRe, psiIm, V, N, L);
    expect(E).toBeCloseTo(0.5 * kx * kx, 4);
  });

  it('energy drift is small over time for free propagation (CAP off)', () => {
    const N = 16;
    const L = 8;
    const dx = L / N;
    const cellVol = dx * dx * dx;
    const coord = createCoordinateArray(N, L);

    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);
    const { gX, gY, gZ, pX, pY, pZ } = allocPacketScratch(N);

    const sigma = 0.6;
    const kx = (2 * Math.PI) / (8 * dx);
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
    const cap = createAbsorbingBoundary(N, 0);
    const dt = 0.02 * dx * dx;
    const V = new Float32Array(size);
    buildKineticExponentials(expK, kArrays, N, dt);
    buildPotentialExponentials(expVh, V, cap, dt, 0);

    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);
    const E0 = calculateEnergy(psiRe, psiIm, V, N, L);
    for (let i = 0; i < 200; i++) {
      timeStep(psiRe, psiIm, expVh, expK, N, sRe, sIm);
    }
    const ET = calculateEnergy(psiRe, psiIm, V, N, L);
    const rel = Math.abs(ET - E0) / Math.max(1e-30, Math.abs(E0));
    expect(rel).toBeLessThan(5e-3);
  });

  it('timeStep is approximately time-reversible for free propagation (CAP off)', () => {
    const N = 16;
    const L = 8;
    const dx = L / N;
    const cellVol = dx * dx * dx;
    const coord = createCoordinateArray(N, L);

    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);
    const { gX, gY, gZ, pX, pY, pZ } = allocPacketScratch(N);

    const sigma = 0.6;
    const kx = (2 * Math.PI) / (8 * dx);
    addPacket3D(
      psiRe, psiIm, coord, N,
      0, 0, 0,
      sigma, sigma, sigma,
      kx, 0, 0,
      1,
      gX, gY, gZ, pX, pY, pZ
    );
    renormalize(psiRe, psiIm, cellVol);

    const psiRe0 = psiRe.slice();
    const psiIm0 = psiIm.slice();

    const kArrays = createKSpaceArrays(N, L);
    const cap = createAbsorbingBoundary(N, 0);
    const dt = 0.02 * dx * dx;
    const V = new Float32Array(size);

    const expK = new Float32Array(2 * size);
    const expVh = new Float32Array(2 * size);
    buildKineticExponentials(expK, kArrays, N, dt);
    buildPotentialExponentials(expVh, V, cap, dt, 0);

    const expKBack = new Float32Array(2 * size);
    const expVhBack = new Float32Array(2 * size);
    buildKineticExponentials(expKBack, kArrays, N, -dt);
    buildPotentialExponentials(expVhBack, V, cap, -dt, 0);

    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);

    for (let i = 0; i < 60; i++) {
      timeStep(psiRe, psiIm, expVh, expK, N, sRe, sIm);
    }
    for (let i = 0; i < 60; i++) {
      timeStep(psiRe, psiIm, expVhBack, expKBack, N, sRe, sIm);
    }

    const err = rmsRel(psiRe, psiIm, psiRe0, psiIm0);
    expect(err).toBeLessThan(2e-4);
  });

  it('calculateEnergy throws when array lengths are inconsistent', () => {
    const N = 2;
    const L = 2;
    const size = N * N * N;
    const psiRe = new Float32Array(size);
    const psiImShort = new Float32Array(size - 1);
    expect(() => calculateEnergy(psiRe, psiImShort, null, N, L)).toThrow(/same length/);

    const psiShort = new Float32Array(size - 1);
    expect(() => calculateEnergy(psiShort, psiShort, null, N, L)).toThrow(/N\^3/);

    const vBad = new Float32Array(size - 1);
    expect(() => calculateEnergy(psiRe, psiRe, vBad, N, L)).toThrow(/V\.length/);
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

  it('free packet center moves with group velocity ~ kx', () => {
    const N = 32;
    const L = 20;
    const dx = L / N;
    const cellVol = dx * dx * dx;
    const coord = createCoordinateArray(N, L);
    const size = N * N * N;

    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);
    const { gX, gY, gZ, pX, pY, pZ } = allocPacketScratch(N);

    const x0 = -L / 4;
    const kx = 1.0;
    addPacket3D(
      psiRe, psiIm, coord, N,
      x0, 0, 0,
      0.7, 0.7, 0.7,
      kx, 0, 0,
      1,
      gX, gY, gZ, pX, pY, pZ
    );
    renormalize(psiRe, psiIm, cellVol);

    const mean0 = meanX(psiRe, psiIm, coord, N, cellVol);

    const kArrays = createKSpaceArrays(N, L);
    const expK = new Float32Array(2 * size);
    const expVh = new Float32Array(2 * size);
    const dt = 0.02 * dx * dx;
    buildKineticExponentials(expK, kArrays, N, dt);
    buildPotentialExponentials(expVh, new Float32Array(size), null, dt, 0);

    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);
    const steps = 60;
    for (let i = 0; i < steps; i++) {
      timeStep(psiRe, psiIm, expVh, expK, N, sRe, sIm);
    }

    const mean1 = meanX(psiRe, psiIm, coord, N, cellVol);
    const shift = mean1 - mean0;
    const expected = kx * dt * steps;
    const rel = Math.abs(shift - expected) / Math.max(1e-6, Math.abs(expected));
    expect(shift).toBeGreaterThan(0);
    expect(rel).toBeLessThan(0.2);
  });

  it('harmonic potential pulls packet center toward origin', () => {
    const N = 32;
    const L = 20;
    const dx = L / N;
    const cellVol = dx * dx * dx;
    const coord = createCoordinateArray(N, L);
    const size = N * N * N;

    const psiRe = new Float32Array(size);
    const psiIm = new Float32Array(size);
    const { gX, gY, gZ, pX, pY, pZ } = allocPacketScratch(N);

    const x0 = -L / 5;
    addPacket3D(
      psiRe, psiIm, coord, N,
      x0, 0, 0,
      0.7, 0.7, 0.7,
      0, 0, 0,
      1,
      gX, gY, gZ, pX, pY, pZ
    );
    renormalize(psiRe, psiIm, cellVol);

    const mean0 = meanX(psiRe, psiIm, coord, N, cellVol);

    const V = new Float32Array(size);
    presetHarmonic(V, coord, N, 1);
    const kArrays = createKSpaceArrays(N, L);
    const expK = new Float32Array(2 * size);
    const expVh = new Float32Array(2 * size);
    const dt = 0.02 * dx * dx;
    buildKineticExponentials(expK, kArrays, N, dt);
    buildPotentialExponentials(expVh, V, null, dt, 0);

    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);
    const steps = 50;
    for (let i = 0; i < steps; i++) {
      timeStep(psiRe, psiIm, expVh, expK, N, sRe, sIm);
    }

    const mean1 = meanX(psiRe, psiIm, coord, N, cellVol);
    expect(mean1).toBeLessThan(0);
    expect(mean1).toBeGreaterThan(mean0);
    expect(Math.abs(mean1)).toBeLessThan(Math.abs(mean0));
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
});
