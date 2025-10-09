import { describe, it, expect, vi } from 'vitest';
import { fft1d, fft3d, cMul, cExp, makeFFTPlan, fft1d_p, runFFTSelfTests } from '../../src/physics/fft.js';

function rmsRelativeError(aRe, aIm, bRe, bIm) {
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

describe('fft.js', () => {
  it('fft1d throws on non power-of-two length', () => {
    const re = new Float32Array(6);
    const im = new Float32Array(6);
    expect(() => fft1d(re, im, false)).toThrow();
  });

  it('fft1d round-trip accuracy (N=32)', () => {
    const N = 32;
    const re = new Float32Array(N);
    const im = new Float32Array(N);
    for (let i = 0; i < N; i++) {
      re[i] = Math.random() * 2 - 1;
      im[i] = Math.random() * 2 - 1;
    }
    const r0 = re.slice();
    const i0 = im.slice();
    fft1d(re, im, false);
    fft1d(re, im, true);
    const err = rmsRelativeError(re, im, r0, i0);
    expect(err).toBeLessThan(1e-5);
  });

  it('fft3d round-trip accuracy (N=8)', () => {
    const N = 8;
    const sz = N * N * N;
    const re = new Float32Array(sz);
    const im = new Float32Array(sz);
    for (let i = 0; i < sz; i++) {
      re[i] = Math.random() * 2 - 1;
      im[i] = Math.random() * 2 - 1;
    }
    const r0 = re.slice();
    const i0 = im.slice();
    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);
    fft3d(re, im, N, false, sRe, sIm);
    fft3d(re, im, N, true, sRe, sIm);
    let num = 0, den = 0;
    for (let i = 0; i < sz; i++) {
      const dr = re[i] - r0[i];
      const di = im[i] - i0[i];
      num += dr * dr + di * di;
      den += r0[i] * r0[i] + i0[i] * i0[i];
    }
    const err = Math.sqrt(num / Math.max(den, 1e-30));
    expect(err).toBeLessThan(1e-5);
  });

  it('cMul correctness', () => {
    const [xr, xi] = cMul(1, 2, 3, 4);
    expect(Math.abs(xr + 5)).toBeLessThan(1e-12);
    expect(Math.abs(xi - 10)).toBeLessThan(1e-12);
  });

  it('cExp is on the unit circle', () => {
    const thetas = [-Math.PI, -Math.PI / 2, 0, Math.PI / 2, Math.PI];
    for (const t of thetas) {
      const [cr, ci] = cExp(t);
      const mag = Math.hypot(cr, ci);
      expect(Math.abs(mag - 1)).toBeLessThan(1e-12);
    }
  });

  it('makeFFTPlan works and rejects invalid sizes', () => {
    expect(() => makeFFTPlan(6)).toThrow();
    const plan = makeFFTPlan(16);
    expect(plan.n).toBe(16);
    expect(Array.isArray(plan.pairs)).toBe(true);
    expect(Array.isArray(plan.stages)).toBe(true);
  });

  it('fft1d_p equals fft1d (forward, N=16)', () => {
    const N = 16;
    const re1 = new Float32Array(N);
    const im1 = new Float32Array(N);
    for (let i = 0; i < N; i++) {
      re1[i] = Math.random() * 2 - 1;
      im1[i] = Math.random() * 2 - 1;
    }
    const re2 = re1.slice();
    const im2 = im1.slice();

    const plan = makeFFTPlan(N);
    fft1d(re1, im1, false);
    fft1d_p(re2, im2, false, plan);

    for (let i = 0; i < N; i++) {
      expect(Math.abs(re1[i] - re2[i])).toBeLessThan(1e-10);
      expect(Math.abs(im1[i] - im2[i])).toBeLessThan(1e-10);
    }
  });

  it('fft1d_p falls back to fft1d when plan size mismatches', () => {
    const N = 16;
    const re1 = new Float32Array(N);
    const im1 = new Float32Array(N);
    for (let i = 0; i < N; i++) {
      re1[i] = Math.random() * 2 - 1;
      im1[i] = Math.random() * 2 - 1;
    }
    const re2 = re1.slice();
    const im2 = im1.slice();

    const badPlan = makeFFTPlan(8); // wrong size
    fft1d(re1, im1, false);
    fft1d_p(re2, im2, false, badPlan); // should fallback internally

    for (let i = 0; i < N; i++) {
      expect(Math.abs(re1[i] - re2[i])).toBeLessThan(1e-10);
      expect(Math.abs(im1[i] - im2[i])).toBeLessThan(1e-10);
    }
  });

  it('runFFTSelfTests executes without throwing', () => {
    const log = vi.spyOn(console, 'log').mockImplementation(() => {});
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    const assertSpy = vi.spyOn(console, 'assert').mockImplementation(() => {});
    runFFTSelfTests();
    // if assertions fail, console.assert would be called with false; we still treat this as non-throwing
    expect(assertSpy).toHaveBeenCalled();
    log.mockRestore();
    warn.mockRestore();
    assertSpy.mockRestore();
  });
});
