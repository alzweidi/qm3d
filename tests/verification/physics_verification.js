import {
  createCoordinateArray,
  createKSpaceArrays,
  buildKineticExponentials,
  buildPotentialExponentials,
  createAbsorbingBoundary,
  timeStep,
  addPacket3D,
  calculateNorm,
  renormalize
} from '../../src/physics/quantum.js';
import { presetFree, presetBoxWell } from '../../src/physics/potentials.js';

// helper to check for NaNs or Infinities
function checkNumericalStability(psiRe, psiIm, step) {
  for (let i = 0; i < psiRe.length; i++) {
    if (!Number.isFinite(psiRe[i]) || !Number.isFinite(psiIm[i])) {
      console.error(`[Stability Fail] NaN/Inf detected at step ${step}, index ${i}`);
      return false;
    }
  }
  return true;
}

// helper to run simulation
function runSimulation(config) {
  const { N, L, dtScale, steps, checkStability, checkNorm, checkBoundary } = config;
  const dx = L / N;
  const dt = dtScale * dx * dx;
  const size = N * N * N;
  const cellVol = dx * dx * dx;

  console.log(`\n--- Running Simulation: ${config.name} ---`);
  console.log(`N=${N}, L=${L}, dt=${dt.toExponential(3)}, steps=${steps}`);

  // arrays
  const psiRe = new Float32Array(size);
  const psiIm = new Float32Array(size);
  const V = new Float32Array(size);
  const expK = new Float32Array(2 * size);
  const expVh = new Float32Array(2 * size);
  const scratchRe = new Float32Array(N);
  const scratchIm = new Float32Array(N);
  
  // scratch arrays for addPacket3D
  const gX = new Float32Array(N);
  const gY = new Float32Array(N);
  const gZ = new Float32Array(N);
  const pX = new Float32Array(N);
  const pY = new Float32Array(N);
  const pZ = new Float32Array(N);

  // setup
  const coord = createCoordinateArray(N, L);
  const kArrays = createKSpaceArrays(N, L);
  buildKineticExponentials(expK, kArrays, N, dt);

  // potential
  if (config.potential === 'box') {
    presetBoxWell(V, coord, N, L);
  } else {
    presetFree(V);
  }
  
  // absorbing boundary
  let capS2 = null;
  if (config.absorb) {
    capS2 = createAbsorbingBoundary(N, 0.15);
  }
  buildPotentialExponentials(expVh, V, capS2, dt, config.absorbStrength || 0);

  // initial packet
  addPacket3D(
    psiRe, psiIm, coord, N,
    config.x0 || 0, config.y0 || 0, config.z0 || 0,
    config.sigma || 1.0, config.sigma || 1.0, config.sigma || 1.0,
    config.kx || 0, config.ky || 0, config.kz || 0,
    1.0,
    gX, gY, gZ, pX, pY, pZ
  );
  renormalize(psiRe, psiIm, cellVol);

  // loop
  let initialNorm = calculateNorm(psiRe, psiIm, cellVol);
  console.log(`Initial Norm: ${initialNorm.toFixed(6)}`);

  for (let s = 1; s <= steps; s++) {
    timeStep(psiRe, psiIm, expVh, expK, N, scratchRe, scratchIm);

    if (checkStability && !checkNumericalStability(psiRe, psiIm, s)) {
      return;
    }

    if (checkNorm && s % 10 === 0) {
      const norm = calculateNorm(psiRe, psiIm, cellVol);
      if (Math.abs(norm - 1.0) > 0.01) {
         console.warn(`[Norm Warning] Step ${s}: Norm = ${norm.toFixed(6)}`);
      }
    }
  }

  const finalNorm = calculateNorm(psiRe, psiIm, cellVol);
  console.log(`Final Norm: ${finalNorm.toFixed(6)}`);
  
  if (checkNorm) {
      if (Math.abs(finalNorm - 1.0) < 0.01) {
          console.log("[PASS] Norm Conservation");
      } else {
          console.error("[FAIL] Norm Conservation");
      }
  }
  
  if (checkBoundary) {
      // simple check: if packet moved towards boundary and norm decreased significantly (if absorbing)
      // or stayed same (if reflecting).
      // for this test, we'll just report the norm.
      console.log(`Boundary Check: Final Norm = ${finalNorm.toFixed(6)} (Expected behavior depends on config)`);
  }
}

// 1. numerical stability test
runSimulation({
  name: "Numerical Stability (Stress Test)",
  N: 32, L: 10, dtScale: 0.1, steps: 100,
  checkStability: true,
  potential: 'free',
  x0: 0, y0: 0, z0: 0,
  sigma: 0.5,
  kx: 5, ky: 5, kz: 5
});

// 2. norm conservation test (free particle)
runSimulation({
  name: "Norm Conservation (Free Particle)",
  N: 32, L: 20, dtScale: 0.05, steps: 50,
  checkNorm: true,
  potential: 'free',
  x0: 0, y0: 0, z0: 0,
  sigma: 2.0,
  kx: 0, ky: 0, kz: 0
});

// 3. boundary conditions (absorbing)
runSimulation({
  name: "Boundary Conditions (Absorbing)",
  N: 32, L: 10, dtScale: 0.05, steps: 100,
  checkBoundary: true,
  absorb: true,
  absorbStrength: 5.0,
  potential: 'free',
  x0: 0, y0: 0, z0: 0,
  sigma: 1.0,
  kx: 10, ky: 0, kz: 0 // moving fast towards +x boundary
});

// 4. boundary conditions (reflecting - box well)
runSimulation({
  name: "Boundary Conditions (Reflecting - Box Well)",
  N: 32, L: 10, dtScale: 0.05, steps: 100,
  checkBoundary: true,
  absorb: false,
  potential: 'box', // infinite potential at boundaries
  x0: 0, y0: 0, z0: 0,
  sigma: 1.0,
  kx: 10, ky: 0, kz: 0 // moving fast towards +x boundary
});