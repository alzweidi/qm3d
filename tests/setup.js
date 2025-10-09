// Vitest setup: deterministic RNG and jest-dom matchers
import '@testing-library/jest-dom/vitest';
import seedrandom from 'seedrandom';
import { expect } from 'vitest';
import { cleanup } from '@testing-library/react';
import { afterEach } from 'vitest';

// make Math.random deterministic across test runs
seedrandom('qm3d-fixed-seed', { global: true });

// ensure DOM is cleaned between tests (robustness across environments)
afterEach(() => {
  cleanup();
});

// custom matcher: compare arrays element-wise with precision (digits)
expect.extend({
  toEqualCloseTo(received, expected, precision = 6) {
    const rec = Array.from(received);
    const exp = Array.from(expected);
    if (rec.length !== exp.length) {
      return {
        pass: false,
        message: () => `length mismatch: got ${rec.length} vs ${exp.length}`,
      };
    }
    const tol = Math.pow(10, -precision);
    for (let i = 0; i < rec.length; i++) {
      if (!Number.isFinite(rec[i]) || !Number.isFinite(exp[i])) {
        if (rec[i] !== exp[i]) {
          return {
            pass: false,
            message: () => `index ${i}: expected ${exp[i]}, received ${rec[i]}`,
          };
        }
      } else if (Math.abs(rec[i] - exp[i]) > tol) {
        return {
          pass: false,
          message: () => `index ${i}: expected ${exp[i]} ± ${tol}, received ${rec[i]}`,
        };
      }
    }
    return { pass: true, message: () => 'arrays are close' };
  },
});
