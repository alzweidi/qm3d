import { describe, it, expect } from 'vitest';
import { createCoordinateArray, addPacket3D } from '../../src/physics/quantum.js';

describe('quantum.js addPacket3D object-arg overload', () => {
  it('matches positional-call output and covers object overload branch', () => {
    const N = 8;
    const L = 8;
    const coord = createCoordinateArray(N, L);

    const size = N * N * N;
    const psiRePos = new Float32Array(size);
    const psiImPos = new Float32Array(size);

    const psiReObj = new Float32Array(size);
    const psiImObj = new Float32Array(size);

    const gX = new Float32Array(N);
    const gY = new Float32Array(N);
    const gZ = new Float32Array(N);
    const pX = new Float32Array(N);
    const pY = new Float32Array(N);
    const pZ = new Float32Array(N);

    const cx = 0, cy = 0, cz = 0;
    const sx = 0.6, sy = 0.6, sz = 0.6;
    const kx = 0.5, ky = 0, kz = 0;

    // positional call (scale = 1.0)
    addPacket3D(
      psiRePos, psiImPos, coord, N,
      cx, cy, cz,
      sx, sy, sz,
      kx, ky, kz,
      1,
      gX, gY, gZ, pX, pY, pZ
    );

    // object overload call (omit scale to hit default value path)
    addPacket3D({
      psiRe: psiReObj,
      psiIm: psiImObj,
      coord,
      N,
      cx, cy, cz,
      sx, sy, sz,
      kx, ky, kz,
      gX, gY, gZ, pX, pY, pZ,
    });

    for (let i = 0; i < size; i++) {
      expect(Math.abs(psiReObj[i] - psiRePos[i])).toBeLessThan(1e-12);
      expect(Math.abs(psiImObj[i] - psiImPos[i])).toBeLessThan(1e-12);
    }
  });
});
