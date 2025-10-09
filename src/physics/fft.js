// fft implementations for quantum wave simulation
// complex number utilities and fast fourier transform functions

/**
 * complex number multiplication
 * @param {number} ar - real part of first complex number
 * @param {number} ai - imaginary part of first complex number  
 * @param {number} br - real part of second complex number
 * @param {number} bi - imaginary part of second complex number
 * @returns {number[]} [real, imaginary] result
 */
export const cMul = (ar, ai, br, bi) => [ar * br - ai * bi, ar * bi + ai * br];

/**
 * complex exponential: e^(i*theta) = cos(theta) + i*sin(theta)
 * @param {number} theta - angle in radians
 * @returns {number[]} [cos(theta), sin(theta)]
 */
export const cExp = (theta) => [Math.cos(theta), Math.sin(theta)];

/**
 * 1D Fast Fourier Transform using Cooley-Tukey algorithm
 * @param {Float32Array} re - real part array (modified in place)
 * @param {Float32Array} im - imaginary part array (modified in place)
 * @param {boolean} inverse - whether to perform inverse FFT
 */
export function fft1d(re, im, inverse = false) {
    const n = re.length;
    if ((n & (n - 1)) !== 0) throw new Error('fft1d requires power‑of‑two length');

    // bit-reversal permutation
    let j = 0;
    for (let i = 0; i < n; i++) {
        if (i < j) {
            const tr = re[i];
            const ti = im[i];
            re[i] = re[j];
            im[i] = im[j];
            re[j] = tr;
            im[j] = ti;
        }
        let m = n >> 1;
        while (m >= 1 && j >= m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    // cooley-tukey fft
    for (let len = 2; len <= n; len <<= 1) {
        const ang = (inverse ? +2 * Math.PI : -2 * Math.PI) / len;
        const wlen_r = Math.cos(ang);
        const wlen_i = Math.sin(ang);

        for (let i = 0; i < n; i += len) {
            let wr = 1, wi = 0;
            const half = len >> 1;

            for (let k = 0; k < half; k++) {
                const a = i + k, b = a + half;
                const u_r = re[a], u_i = im[a];
                const v_r = re[b] * wr - im[b] * wi;
                const v_i = re[b] * wi + im[b] * wr;

                re[a] = u_r + v_r;
                im[a] = u_i + v_i;
                re[b] = u_r - v_r;
                im[b] = u_i - v_i;

                const nwr = wr * wlen_r - wi * wlen_i;
                wi = wr * wlen_i + wi * wlen_r;
                wr = nwr;
            }
        }
    }

    // normalise for inverse transform
    if (inverse) {
        for (let i = 0; i < n; i++) {
            re[i] /= n;
            im[i] /= n;
        }
    }
}

/**
 * create an fft plan for optimized repeated transforms
 * @param {number} n - size of transform (must be power of 2)
 * @returns {Object} fft plan with precomputed values
 */
export function makeFFTPlan(n) {
    if ((n & (n - 1)) !== 0) throw new Error('fft plan requires power‑of‑two length');

    // precompute bit-reversal pairs
    const pairs = [];
    let j = 0;
    for (let i = 0; i < n; i++) {
        if (i < j) pairs.push([i, j]);
        let m = n >> 1;
        while (m >= 1 && j >= m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    // precompute twiddle factors for each stage
    const stages = [];
    for (let len = 2; len <= n; len <<= 1) {
        const angF = -2 * Math.PI / len;
        const angI = +2 * Math.PI / len;
        stages.push({
            len,
            half: len >> 1,
            wlen_r_f: Math.cos(angF),
            wlen_i_f: Math.sin(angF),
            wlen_r_i: Math.cos(angI),
            wlen_i_i: Math.sin(angI)
        });
    }

    return { n, pairs, stages };
}

/**
 * planned 1d fft using precomputed plan for better performance
 * @param {Float32Array} re - real part array
 * @param {Float32Array} im - imaginary part array  
 * @param {boolean} inverse - whether to perform inverse fft
 * @param {Object} plan - precomputed fft plan
 */
export function fft1d_p(re, im, inverse, plan) {
    const n = re.length;
    if (!plan || plan.n !== n) {
        fft1d(re, im, inverse);
        return;
    }

    // apply bit-reversal using precomputed pairs
    for (let t = 0; t < plan.pairs.length; t++) {
        const a = plan.pairs[t][0], b = plan.pairs[t][1];
        const tr = re[a], ti = im[a];
        re[a] = re[b];
        im[a] = im[b];
        re[b] = tr;
        im[b] = ti;
    }

    // apply stages using precomputed twiddle factors
    for (let s = 0; s < plan.stages.length; s++) {
        const st = plan.stages[s];
        const wlen_r = inverse ? st.wlen_r_i : st.wlen_r_f;
        const wlen_i = inverse ? st.wlen_i_i : st.wlen_i_f;
        const len = st.len, half = st.half;

        for (let i = 0; i < n; i += len) {
            let wr = 1, wi = 0;
            for (let k = 0; k < half; k++) {
                const a = i + k, b = a + half;
                const u_r = re[a], u_i = im[a];
                const v_r = re[b] * wr - im[b] * wi;
                const v_i = re[b] * wi + im[b] * wr;

                re[a] = u_r + v_r;
                im[a] = u_i + v_i;
                re[b] = u_r - v_r;
                im[b] = u_i - v_i;

                const nwr = wr * wlen_r - wi * wlen_i;
                wi = wr * wlen_i + wi * wlen_r;
                wr = nwr;
            }
        }
    }

    // normalise for inverse transform
    if (inverse) {
        for (let i = 0; i < n; i++) {
            re[i] /= n;
            im[i] /= n;
        }
    }
}

/**
 * 3d fast fourier transform
 * applies 1d fft along each dimension sequentially
 * @param {Float32Array} re - real part array (N³ elements)
 * @param {Float32Array} im - imaginary part array (N³ elements)
 * @param {number} N - grid size (N×N×N)
 * @param {boolean} inverse - whether to perform inverse fft
 * @param {Float32Array} lineRe - scratch array for 1d transforms
 * @param {Float32Array} lineIm - scratch array for 1d transforms
 * @param {Object} plan - optional fft plan for optimization
 */
export function fft3d(re, im, N, inverse, lineRe, lineIm, plan) {
    const idx = (x, y, z) => x + N * (y + N * z);

    // transform along x direction
    for (let z = 0; z < N; z++) {
        for (let y = 0; y < N; y++) {
            for (let x = 0; x < N; x++) {
                const id = idx(x, y, z);
                lineRe[x] = re[id];
                lineIm[x] = im[id];
            }
            fft1d(lineRe, lineIm, inverse);
            for (let x = 0; x < N; x++) {
                const id = idx(x, y, z);
                re[id] = lineRe[x];
                im[id] = lineIm[x];
            }
        }
    }

    // transform along y direction
    for (let z = 0; z < N; z++) {
        for (let x = 0; x < N; x++) {
            for (let y = 0; y < N; y++) {
                const id = idx(x, y, z);
                lineRe[y] = re[id];
                lineIm[y] = im[id];
            }
            fft1d(lineRe, lineIm, inverse);
            for (let y = 0; y < N; y++) {
                const id = idx(x, y, z);
                re[id] = lineRe[y];
                im[id] = lineIm[y];
            }
        }
    }

    // transform along z direction
    for (let y = 0; y < N; y++) {
        for (let x = 0; x < N; x++) {
            for (let z = 0; z < N; z++) {
                const id = idx(x, y, z);
                lineRe[z] = re[id];
                lineIm[z] = im[id];
            }
            fft1d(lineRe, lineIm, inverse);
            for (let z = 0; z < N; z++) {
                const id = idx(x, y, z);
                re[id] = lineRe[z];
                im[id] = lineIm[z];
            }
        }
    }
}

/**
 * run self-tests for fft implementation
 * verifies round-trip accuracy and complex arithmetic
 */
export function runFFTSelfTests() {
    try {
        // test 1d fft round-trip
        const n = 32;
        const re = new Float32Array(n);
        const im = new Float32Array(n);

        for (let i = 0; i < n; i++) {
            re[i] = Math.random() * 2 - 1;
            im[i] = Math.random() * 2 - 1;
        }

        const r1 = re.slice(), i1 = im.slice();
        fft1d(re, im, false);
        fft1d(re, im, true);

        let err = 0, norm = 0;
        for (let i = 0; i < n; i++) {
            const dr = re[i] - r1[i], di = im[i] - i1[i];
            err += dr * dr + di * di;
            norm += r1[i] * r1[i] + i1[i] * i1[i];
        }
        console.assert(Math.sqrt(err / (norm + 1e-12)) < 1e-5, "FFT(1D) round‑trip error too large");

        // test complex multiplication
        const cm = cMul(1, 2, 3, 4);
        console.assert(Math.abs(cm[0] + 5) < 1e-6 && Math.abs(cm[1] - 10) < 1e-6, "cMul incorrect");

        // test complex exponential
        const dtt = 0.01, Vc = 2.0;
        const eh = cExp(-0.5 * dtt * Vc);
        console.assert(Math.abs(eh[0] * eh[0] + eh[1] * eh[1] - 1) < 1e-6, "exp phase not unit magnitude");

        // test 3d fft round-trip
        const NN = 8;
        const sz = NN * NN * NN;
        const Re3 = new Float32Array(sz);
        const Im3 = new Float32Array(sz);

        for (let i = 0; i < sz; i++) {
            Re3[i] = Math.random() * 2 - 1;
            Im3[i] = Math.random() * 2 - 1;
        }

        const Re3b = Re3.slice(), Im3b = Im3.slice();
        const sRe = new Float32Array(NN), sIm = new Float32Array(NN);

        fft3d(Re3, Im3, NN, false, sRe, sIm);
        fft3d(Re3, Im3, NN, true, sRe, sIm);

        let e3 = 0, n3 = 0;
        for (let i = 0; i < sz; i++) {
            const dr = Re3[i] - Re3b[i], di = Im3[i] - Im3b[i];
            e3 += dr * dr + di * di;
            n3 += Re3b[i] * Re3b[i] + Im3b[i] * Im3b[i];
        }
        console.assert(Math.sqrt(e3 / (n3 + 1e-12)) < 1e-5, "FFT(3D) round‑trip error too large");

        console.log("✓ FFT self-tests passed");
    } catch (e) {
        console.warn("FFT self‑tests failed:", e);
    }
}