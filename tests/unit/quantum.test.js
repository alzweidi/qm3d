import { describe, it, expect } from 'vitest';
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
} from '../../src/physics/quantum.js';

// Helpers
function idx(x, y, z, N) {
  return y + N * (x + N * z); // CAREFUL: check implementation order
}
// In fft3d the index helper is idx(x, y, z) => x + N*(y + N*z)
// For physics arrays in quantum.js, linearization uses yOff + x within loops with z as outer:
// id = z*N*N + y*N + x  => id = x + N*(y + N*z) which matches fft.js helper.
function lin(x, y, z, N) {
  return x + N * (y + N * z);
}

function rmsRelativeError(a, b) {
  let num = 0;
  let den = 0;
  for (let i = 0; i < a.length; i++) {
    const d = a[i] - b[i];
    num += d * d;
    den += b[i] * b[i];
  }
  return Math.sqrt(num / Math.max(den, 1e-30));
}

describe('quantum.js', () => {
  it('createCoordinateArray produces centered cell coordinates', () => {
    const arr = createCoordinateArray(4, 2);
    expect(Array.from(arr)).toEqualCloseTo([-0.75, -0.25, 0.25, 0.75], 10);
  });

  it('createKSpaceArrays maps to integer grid when L=2π', () => {
    const N = 8;
    const L = 2 * Math.PI;
    const { kx2, ky2, kz2 } = createKSpaceArrays(N, L);
    // Expected frequencies: [0,1,2,3,-4,-3,-2,-1]
    const expectedSq = [0,1,4,9,16,9,4,1];
    expect(Array.from(kx2)).toEqualCloseTo(expectedSq, 10);
    expect(Array.from(ky2)).toEqualCloseTo(expectedSq, 10);
    expect(Array.from(kz2)).toEqualCloseTo(expectedSq, 10);
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
    // V: [0,1,-2, 0,0,0, 0,0] to test different angles
    V[0] = 0; V[1] = 1; V[2] = -2;
    const dt = 0.1;
    const absorbStrength = 3.0;
    const expVh = new Float32Array(2 * size);

    // Case 1: no CAP
    buildPotentialExponentials(expVh, V, null, dt, absorbStrength);
    let theta = -0.5 * V[1] * dt;
    let j = (1 << 1);
    // Float32 rounding on expVh components
    expect(expVh[j]).toBeCloseTo(Math.cos(theta), 7);
    expect(expVh[j + 1]).toBeCloseTo(Math.sin(theta), 7);

    // Case 2: with CAP
    const capS2 = new Float32Array(size);
    capS2[2] = 2.0; // non-zero cap only at index 2
    buildPotentialExponentials(expVh, V, capS2, dt, absorbStrength);
    theta = -0.5 * V[2] * dt;
    const decay = Math.exp(-0.5 * absorbStrength * capS2[2] * dt);
    j = (2 << 1);
    const mag = Math.hypot(expVh[j], expVh[j + 1]);
    // Use relative tolerance suitable for Float32 computations
    const rel = Math.abs(mag - decay) / Math.max(decay, 1e-12);
    expect(rel).toBeLessThan(1e-7);
    // Angle still follows theta
    const ang = Math.atan2(expVh[j + 1], expVh[j]);
    expect(ang).toBeCloseTo(theta, 7); // allow some wrap/precision tolerance
  });

  it('createAbsorbingBoundary interior is zero and edges positive/symmetric', () => {
    const N = 8;
    const cap = createAbsorbingBoundary(N, 0.25); // nAbs=2

    const interiorIndex = lin(4, 4, 4, N);
    expect(cap[interiorIndex]).toBe(0);

    const edgeIndex = lin(0, 4, 4, N);
    expect(cap[edgeIndex]).toBeGreaterThan(0);

    // Symmetry checks
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
    // After multiply: index 0 stays (1,0); index 1 rotates by -π/2
    expect(psiRe[0]).toBeCloseTo(1, 12);
    expect(psiIm[0]).toBeCloseTo(0, 12);
    // (0.5 - i) * (0 - i) = (-1) * 0.5 + ... compute explicitly
    // Let (a+ib)*(c+id) with c=0, d=-1 => real = a*c - b*d = 0 - (-1)*(-1)???
    // Better compute numerically from initial
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
    for (let i = 0; i < size; i++) {
      psiRe[i] = Math.random() * 2 - 1;
      psiIm[i] = Math.random() * 2 - 1;
    }
    const kArrays = createKSpaceArrays(N, L);
    const expK = new Float32Array(2 * size);
    buildKineticExponentials(expK, kArrays, N, 0.02);

    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);

    const norm0 = calculateNorm(psiRe, psiIm, 1.0);
    kineticFullStep(psiRe, psiIm, expK, N, sRe, sIm);
    const norm1 = calculateNorm(psiRe, psiIm, 1.0);
    expect(Math.abs(norm1 - norm0)).toBeLessThan(1e-6);
  });

  it('timeStep preserves norm for V=0 and decreases with CAP', () => {
    const N = 4;
    const L = 10;
    const size = N * N * N;

    // Initial random psi
    const psiReA = new Float32Array(size);
    const psiImA = new Float32Array(size);
    for (let i = 0; i < size; i++) {
      psiReA[i] = Math.random() * 2 - 1;
      psiImA[i] = Math.random() * 2 - 1;
    }

    const kArrays = createKSpaceArrays(N, L);
    const expK = new Float32Array(2 * size);
    const expVh = new Float32Array(2 * size);

    // Case 1: V=0, no CAP
    const V = new Float32Array(size); // zeros
    buildKineticExponentials(expK, kArrays, N, 0.02);
    buildPotentialExponentials(expVh, V, null, 0.02, 0);

    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);

    const norm0 = calculateNorm(psiReA, psiImA, 1.0);
    timeStep(psiReA, psiImA, expVh, expK, N, sRe, sIm);
    const norm1 = calculateNorm(psiReA, psiImA, 1.0);
    expect(Math.abs(norm1 - norm0)).toBeLessThan(1e-6);

    // Case 2: with CAP
    const psiReB = psiReA.slice();
    const psiImB = psiImA.slice();
    const capS2 = createAbsorbingBoundary(N, 0.25);
    buildPotentialExponentials(expVh, V, capS2, 0.02, 3.0);

    const normB0 = calculateNorm(psiReB, psiImB, 1.0);
    timeStep(psiReB, psiImB, expVh, expK, N, sRe, sIm);
    const normB1 = calculateNorm(psiReB, psiImB, 1.0);
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
    let num = 0, den = 0;
    for (let i = 0; i < size; i++) {
      const dr = psiRe[i] - r0[i];
      const di = psiIm[i] - i0[i];
      num += dr * dr + di * di;
      den += r0[i] * r0[i] + i0[i] * i0[i];
    }
    const rel = Math.sqrt(num / Math.max(den, 1e-30));
    expect(rel).toBeLessThan(1e-7);
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
      1.0,
      gX, gY, gZ, pX, pY, pZ
    );

    // Find argmax index
    let maxV = -Infinity, maxI = -1;
    for (let i = 0; i < size; i++) {
      const v = psiRe[i] * psiRe[i] + psiIm[i] * psiIm[i];
      if (v > maxV) { maxV = v; maxI = i; }
    }
    // Decode x,y,z
    const x = maxI % N;
    const y = Math.floor(maxI / N) % N;
    const z = Math.floor(maxI / (N * N));
    const centers = [N/2 - 1, N/2];
    expect(centers).toContain(x);
    expect(centers).toContain(y);
    expect(centers).toContain(z);

    // Linearity: second add with scale=2 doubles amplitudes
    const psiRe2 = psiRe.slice();
    const psiIm2 = psiIm.slice();
    addPacket3D(
      psiRe2, psiIm2, coord, N,
      0, 0, 0,
      0.6, 0.6, 0.6,
      0, 0, 0,
      2.0,
      gX, gY, gZ, pX, pY, pZ
    );
    let ratioOK = true;
    for (let i = 0; i < size; i++) {
      if (psiRe[i] !== 0 || psiIm[i] !== 0) {
        // Expect the increment to be double
        const dRe = psiRe2[i] - psiRe[i];
        const dIm = psiIm2[i] - psiIm[i];
        const r = Math.hypot(dRe, dIm);
        const base = Math.hypot(psiRe[i], psiIm[i]);
        // Because first packet + second packet => total amplitude is base + 2*base at center, but
        // comparing increment to base is approximately 2 when phases are same (k=0).
        if (base > 0) {
          const factor = r / base;
          if (Math.abs(factor - 2) > 1e-6) { ratioOK = false; break; }
        }
      }
    }
    expect(ratioOK).toBe(true);
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

    addPacket3D(psiRe, psiIm, coord, N, 0, 0, 0, 0.6, 0.6, 0.6, 0, 0, 0, 1.0, gX, gY, gZ, pX, pY, pZ);
    const cellVol = Math.pow(L / N, 3);
    renormalize(psiRe, psiIm, cellVol);
    const norm = calculateNorm(psiRe, psiIm, cellVol);
    expect(Math.abs(norm - 1)).toBeLessThan(1e-6);

    // Zero state
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
      1.0,
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
});
