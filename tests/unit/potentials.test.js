import { describe, it, expect } from 'vitest';
import {
  presetFree,
  presetPlaneBarrier,
  presetBoxWell,
  presetSphere,
  presetHarmonic,
} from '../../src/physics/potentials.js';
import { createCoordinateArray } from '../../src/physics/quantum.js';

function lin(x, y, z, N) {
  return x + N * (y + N * z);
}

describe('potentials.js', () => {
  it('presetFree sets all zeros', () => {
    const N = 4;
    const size = N * N * N;
    const V = new Float32Array(size);
    V.fill(7);
    presetFree(V);
    for (let i = 0; i < size; i++) expect(V[i]).toBe(0);
  });

  it('presetPlaneBarrier applies V0 for x>0 only', () => {
    const N = 4;
    const L = 2;
    const coord = createCoordinateArray(N, L);
    const size = N * N * N;
    const V = new Float32Array(size);
    const V0 = 5;
    presetPlaneBarrier(V, coord, N, V0);

    for (let z = 0; z < N; z++) {
      for (let y = 0; y < N; y++) {
        for (let x = 0; x < N; x++) {
          const id = lin(x, y, z, N);
          if (coord[x] > 0) expect(V[id]).toBe(V0);
          else expect(V[id]).toBe(0);
        }
      }
    }
  });

  it('presetBoxWell sets V0 inside central cube only', () => {
    const N = 6;
    const L = 12;
    const coord = createCoordinateArray(N, L);
    const V = new Float32Array(N * N * N);
    const V0 = -6;
    presetBoxWell(V, coord, N, L, V0);

    const w = L * 0.4;
    for (let z = 0; z < N; z++) {
      for (let y = 0; y < N; y++) {
        for (let x = 0; x < N; x++) {
          const id = lin(x, y, z, N);
          const inside = Math.abs(coord[x]) < w / 2 && Math.abs(coord[y]) < w / 2 && Math.abs(coord[z]) < w / 2;
          expect(V[id]).toBe(inside ? V0 : 0);
        }
      }
    }
  });

  it('presetSphere sets V0 for r<R only', () => {
    const N = 6;
    const L = 12;
    const coord = createCoordinateArray(N, L);
    const V = new Float32Array(N * N * N);
    const V0 = -8;
    presetSphere(V, coord, N, L, V0);

    const R = L * 0.25;
    const R2 = R * R;
    for (let z = 0; z < N; z++) {
      for (let y = 0; y < N; y++) {
        for (let x = 0; x < N; x++) {
          const id = lin(x, y, z, N);
          const r2 = coord[x] * coord[x] + coord[y] * coord[y] + coord[z] * coord[z];
          expect(V[id]).toBe(r2 < R2 ? V0 : 0);
        }
      }
    }
  });

  it('presetHarmonic yields 0.5*omega^2*r^2 and scales with omega^2', () => {
    const N = 4;
    const L = 8;
    const coord = createCoordinateArray(N, L);
    const V1 = new Float32Array(N * N * N);
    const V2 = new Float32Array(N * N * N);
    const w1 = 1;
    const w2 = 2;

    presetHarmonic(V1, coord, N, w1);
    presetHarmonic(V2, coord, N, w2);

    for (let i = 0; i < V1.length; i++) {
      expect(Math.abs(V2[i] - 4 * V1[i])).toBeLessThan(1e-12);
    }
  });
});
