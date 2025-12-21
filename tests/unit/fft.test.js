import { describe, it, expect, vi } from 'vitest';
import { rmsRelativeErrorComplex as rmsRelativeError, fillRandomPair, expectFft1dPlanEquivalence, compute3dRoundTripError } from '../utils/testUtils.js';
import { fft1d, fft3d, cMul, cExp, makeFFTPlan, fft1d_p, runFFTSelfTests } from '../../src/physics/fft.js';

// rmsRelativeError and fillRandomPair moved to tests/utils/testUtils.js

describe('fft.js', () => {
  it('fft1d throws on non power-of-two length', () => {
    const re = new Float32Array(6);
    const im = new Float32Array(6);
    expect(() => fft1d(re, im, false)).toThrow();
  });

  it('fft1d throws on mismatched re/im lengths', () => {
    const re = new Float32Array(8);
    const im = new Float32Array(4);
    expect(() => fft1d(re, im, false)).toThrow(/same length/i);
  });

  it('fft1d_p throws on mismatched re/im lengths', () => {
    const re = new Float32Array(8);
    const im = new Float32Array(4);
    const plan = makeFFTPlan(8);
    expect(() => fft1d_p(re, im, false, plan)).toThrow(/same length/i);
  });

  it('fft3d validates array lengths and scratch size', () => {
    const N = 4;
    const size = N * N * N;

    const reBad = new Float32Array(size);
    const imBad = new Float32Array(size - 1);
    const sRe = new Float32Array(N);
    const sIm = new Float32Array(N);
    expect(() => fft3d(reBad, imBad, N, false, sRe, sIm)).toThrow(/same length/i);

    const reShort = new Float32Array(size - 2);
    const imShort = new Float32Array(size - 2);
    expect(() => fft3d(reShort, imShort, N, false, sRe, sIm)).toThrow(/N\^3/i);

    const re = new Float32Array(size);
    const im = new Float32Array(size);
    const sReShort = new Float32Array(N - 1);
    expect(() => fft3d(re, im, N, false, sReShort, sIm)).toThrow(/lineRe\/lineIm length/i);
  });

  it('fft1d of a complex exponential peaks at the expected bin', () => {
    const N = 16;
    const m = 3;
    const re = new Float32Array(N);
    const im = new Float32Array(N);
    for (let n = 0; n < N; n++) {
      const ang = (2 * Math.PI * m * n) / N;
      re[n] = Math.cos(ang);
      im[n] = Math.sin(ang);
    }
    fft1d(re, im, false);
    let maxMag = -Infinity;
    let maxIdx = -1;
    let maxOther = 0;
    for (let k = 0; k < N; k++) {
      const mag = Math.hypot(re[k], im[k]);
      if (mag > maxMag) { maxMag = mag; maxIdx = k; }
      if (k !== m) maxOther = Math.max(maxOther, mag);
    }
    expect(maxIdx).toBe(m);
    expect(maxMag).toBeCloseTo(N, 6);
    expect(maxOther).toBeLessThan(1e-4);
  });

  it('fft1d_p equals fft1d (inverse, N=16)', () => {
    expectFft1dPlanEquivalence(16, 'inverse');
  });

  it('fft1d round-trip accuracy (N=32)', () => {
    const N = 32;
    const re = new Float32Array(N);
    const im = new Float32Array(N);
    fillRandomPair(re, im);
    const r0 = re.slice();
    const i0 = im.slice();
    fft1d(re, im, false);
    fft1d(re, im, true);
    const err = rmsRelativeError(re, im, r0, i0);
    expect(err).toBeLessThan(1e-5);
  });

  it('fft3d round-trip accuracy (N=8)', () => {
    const err = compute3dRoundTripError(8);
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
    expectFft1dPlanEquivalence(16, 'forward');
  });

  it('fft1d_p equals fft1d (inverse, N=16)', () => {
    const N = 16;
    const reF = new Float32Array(N);
    const imF = new Float32Array(N);
    for (let i = 0; i < N; i++) {
      reF[i] = Math.random() * 2 - 1;
      imF[i] = Math.random() * 2 - 1;
    }

    const re1 = reF.slice();
    const im1 = imF.slice();
    const re2 = reF.slice();
    const im2 = imF.slice();

    const plan = makeFFTPlan(N);

    // forward both using canonical fft
    fft1d(re1, im1, false);
    fft1d(re2, im2, false);

    // inverse: compare planned vs canonical (covers inverse normalization branch)
    fft1d(re1, im1, true);
    fft1d_p(re2, im2, true, plan);

    for (let i = 0; i < N; i++) {
      expect(Math.abs(re1[i] - re2[i])).toBeLessThan(1e-10);
      expect(Math.abs(im1[i] - im2[i])).toBeLessThan(1e-10);
    }
  });

  it('fft1d_p falls back to fft1d when plan size mismatches', () => {
    const N = 16;
    const re1 = new Float32Array(N);
    const im1 = new Float32Array(N);
    fillRandomPair(re1, im1);
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
    const err = vi.spyOn(console, 'error').mockImplementation(() => {});
    expect(() => runFFTSelfTests()).not.toThrow();
    expect(log).toHaveBeenCalled();
    log.mockRestore();
    err.mockRestore();
  });

  it('runFFTSelfTests propagates internal failure (fail hard)', () => {
    const log = vi.spyOn(console, 'log').mockImplementation(() => {});
    const err = vi.spyOn(console, 'error').mockImplementation(() => {});
    const rnd = vi.spyOn(Math, 'random').mockImplementation(() => { throw new Error('rand fail'); });
    expect(() => runFFTSelfTests()).toThrow();
    expect(err).toHaveBeenCalled();
    // cleanup
    rnd.mockRestore();
    log.mockRestore();
    err.mockRestore();
  });
});
