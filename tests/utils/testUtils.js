// Shared test utilities to reduce duplication across unit tests
// Keep helpers here so SonarQube doesn't flag duplicate blocks in tests

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
