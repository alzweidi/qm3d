// @vitest-environment jsdom
import { describe, it, expect } from 'vitest';
import { estimateVsyncMs, computeDroppedFrames, buildHistogram } from '../../src/profiler/rafMonitor.js';

describe('profiler algorithms', () => {
  it('estimateVsyncMs returns a plausible vsync for 60Hz-like samples', () => {
    const dts = [16.2, 16.5, 16.6, 16.7, 16.8, 16.9, 17.0, 16.3, 16.4];
    const vs = estimateVsyncMs(dts);
    expect(vs).toBeGreaterThan(8);
    expect(vs).toBeLessThan(40);
    // allow tolerance for robust-estimation method
    expect(vs).toBeGreaterThan(15.5);
    expect(vs).toBeLessThan(17.5);
  });

  it('computeDroppedFrames yields 0/1/2 for ~1x/~2x/~3x frame durations', () => {
    const v = 16.67;
    expect(computeDroppedFrames(16.6, v)).toBe(0);
    expect(computeDroppedFrames(33.4, v)).toBe(1);
    expect(computeDroppedFrames(50.0, v)).toBe(2);
  });

  it('buildHistogram buckets values by edges (half-open intervals)', () => {
    const dts = [8, 10, 16.7, 17, 40, 100];
    const edges = [0, 12, 16.7, 33.3, 1000];
    const counts = buildHistogram(dts, edges);
    // [0-12), [12-16.7), [16.7-33.3), [33.3-1000)
    expect(Array.from(counts)).toEqual([2, 0, 2, 2]);
  });
});
