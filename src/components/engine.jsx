// main quantum wave engine component
// orchestrates physics simulation, rendering, and user interface

import React, { useCallback, useEffect, useMemo, useRef, useState } from "react";
import Controls from './controls.jsx';
import { runFFTSelfTests } from '../physics/fft.js';
import { 
  presetFree, 
  presetPlaneBarrier, 
  presetBoxWell, 
  presetSphere, 
  presetHarmonic 
} from '../physics/potentials.js';
import {
  createCoordinateArray,
  createKSpaceArrays,
  buildKineticExponentials,
  buildPotentialExponentials,
  createAbsorbingBoundary,
  timeStep,
  addPacket3D,
  renormalize,
} from '../physics/quantum.js';
import {
  initialiseThreeJS,
  createQuantumPointCloud,
  createBoundingBox,
  updatePointCloud,
  handleResize,
  renderScene,
  disposeThreeJS
} from '../rendering/visualisation.js';
import { captureScreenshot, startRecording } from '../rendering/capture.js';

export default function QuantumWaveEngine() {
  // simulation parameters
  const [N, setN] = useState(32);
  const [L, setL] = useState(10);
  const dx = useMemo(() => L / N, [L, N]);
  const [dtScale, setDtScale] = useState(0.08);
  const dt = dtScale * dx * dx;
  const kMax = useMemo(() => Math.PI / dx * 0.9, [dx]);
  const [stepsPerFrame, setStepsPerFrame] = useState(1);

  // absorbing boundaries
  const [absorbFrac, setAbsorbFrac] = useState(0.15);
  const [absorbStrength, setAbsorbStrength] = useState(3);

  // wave packet parameters
  const [sigma, setSigma] = useState(0.6);
  const targetPPW = 8;
  const kDefault = useMemo(() => {
    const kMaxSafe = Math.PI / dx * 0.6;
    const kPpw = (2 * Math.PI) / (targetPPW * dx);
    return Math.min(kMaxSafe, kPpw);
  }, [dx]);
  const [k0x, setK0x] = useState(kDefault);
  const [k0y, setK0y] = useState(0);
  const [k0z, setK0z] = useState(0);
  const [amp, setAmp] = useState(1);
  const [x0, setX0] = useState(-L / 4);
  const [y0, setY0] = useState(0);
  const [z0, setZ0] = useState(0);

  // visualisation
  const [densityScale, setDensityScale] = useState(1.2);
  const [showPhase, setShowPhase] = useState(true);
  const [running, setRunning] = useState(true);
  const [isRecording, setIsRecording] = useState(false);
  const [isCapturingShot, setIsCapturingShot] = useState(false);
  const [recordElapsedSec, setRecordElapsedSec] = useState(0);

  // refs for dom and three.js
  const mountRef = useRef(null);
  const threeRefs = useRef({});
  const maxDRef = useRef(1e-9);
  const recorderRef = useRef(null);
  const recordTimerRef = useRef(null);

  const recTimeLabel = useMemo(() => {
    const s = recordElapsedSec;
    const mm = String(Math.floor(s / 60)).padStart(2, '0');
    const ss = String(s % 60).padStart(2, '0');
    return `${mm}:${ss}`;
  }, [recordElapsedSec]);

  const getCaptureSize = useCallback((baseWidth) => {
    const el = mountRef.current;
    if (!el) return { width: baseWidth, height: Math.round(baseWidth * 9 / 16) };
    const w = el.clientWidth || 1;
    const h = el.clientHeight || 480;
    const aspect = h / Math.max(1, w);
    return { width: baseWidth, height: Math.round(baseWidth * aspect) };
  }, []);

  // physics arrays
  const size = useMemo(() => N * N * N, [N]);
  const psiRe = useRef(new Float32Array(size));
  const psiIm = useRef(new Float32Array(size));
  const V = useRef(new Float32Array(size));
  const expK = useRef(new Float32Array(2 * size));
  const expVh = useRef(new Float32Array(2 * size));
  const capS2Ref = useRef(new Float32Array(size));

  // scratch arrays for computations
  const gXRef = useRef(new Float32Array(N));
  const gYRef = useRef(new Float32Array(N));
  const gZRef = useRef(new Float32Array(N));
  const pXRef = useRef(new Float32Array(N));
  const pYRef = useRef(new Float32Array(N));
  const pZRef = useRef(new Float32Array(N));
  const scratchRe = useRef(new Float32Array(N));
  const scratchIm = useRef(new Float32Array(N));

  // coordinate and k-space arrays
  const coord = useMemo(() => createCoordinateArray(N, L), [N, L]);
  const kArrays = useMemo(() => createKSpaceArrays(N, L), [N, L]);

  const updateVisualisation = useCallback(() => {
    const { ampAttr, phaseAttr } = threeRefs.current;
    if (ampAttr && phaseAttr) {
      updatePointCloud(
        psiRe.current,
        psiIm.current,
        ampAttr,
        phaseAttr,
        maxDRef,
        densityScale,
        showPhase
      );
    }
  }, [densityScale, showPhase]);

  const rebuildPotentialExponentials = useCallback(() => {
    buildPotentialExponentials(
      expVh.current,
      V.current,
      capS2Ref.current,
      dt,
      absorbStrength
    );
  }, [dt, absorbStrength]);

  const initialiseVisualisation = useCallback(() => {
    disposeThreeJS(threeRefs.current);

    const threeObjects = initialiseThreeJS(mountRef.current, L);
    const pointCloudObjects = createQuantumPointCloud(
      coord, N, threeObjects.maxPointSize, true
    );
    const boxHelper = createBoundingBox(L);

    threeObjects.scene.add(pointCloudObjects.points);
    threeObjects.scene.add(boxHelper);

    threeRefs.current = {
      ...threeObjects,
      ...pointCloudObjects,
      boxHelper
    };

    handleResize(
      threeObjects.renderer,
      threeObjects.camera,
      mountRef.current,
      pointCloudObjects.points
    );
  }, [L, coord, N]);

  // update potential exponentials
  useEffect(() => {
    rebuildPotentialExponentials();
  }, [dt, absorbStrength, rebuildPotentialExponentials]);

  // initialise arrays when n or l changes
  useEffect(() => {
    const newSize = N * N * N;
    psiRe.current = new Float32Array(newSize);
    psiIm.current = new Float32Array(newSize);
    V.current = new Float32Array(newSize);
    expK.current = new Float32Array(2 * newSize);
    expVh.current = new Float32Array(2 * newSize);

    gXRef.current = new Float32Array(N);
    gYRef.current = new Float32Array(N);
    gZRef.current = new Float32Array(N);
    pXRef.current = new Float32Array(N);
    pYRef.current = new Float32Array(N);
    pZRef.current = new Float32Array(N);
    scratchRe.current = new Float32Array(N);
    scratchIm.current = new Float32Array(N);

    V.current.fill(0);
    psiRe.current.fill(0);
    psiIm.current.fill(0);
    maxDRef.current = 1e-9;

    initialiseVisualisation();
  }, [N, L, initialiseVisualisation]);

  // clamp packet center positions to domain
  useEffect(() => {
    const clamp = (v) => Math.max(-L/2, Math.min(L/2, v));
    setX0(x => clamp(x));
    setY0(y => clamp(y));
    setZ0(z => clamp(z));
  }, [L]);

  useEffect(() => {
    const clamp = (v) => Math.max(-kMax, Math.min(kMax, v));
    setK0x(k => clamp(k));
    setK0y(k => clamp(k));
    setK0z(k => clamp(k));
  }, [kMax]);

  const prevKDefaultRef = useRef(kDefault);
  useEffect(() => {
    const eps = 1e-6;
    setK0x(v => (Math.abs(v - prevKDefaultRef.current) < eps ? kDefault : v));
    prevKDefaultRef.current = kDefault;
  }, [kDefault]);

  // update absorbing boundaries
  useEffect(() => {
    capS2Ref.current = createAbsorbingBoundary(N, absorbFrac);
    rebuildPotentialExponentials();
  }, [N, absorbFrac, rebuildPotentialExponentials]);

  // update kinetic exponentials
  useEffect(() => {
    buildKineticExponentials(expK.current, kArrays, N, dt);
  }, [N, L, dt, kArrays]);

  const renderOnce = useCallback(() => {
    const { renderer, scene, camera, controls } = threeRefs.current;
    renderScene(renderer, scene, camera, controls);
  }, []);

  const renderDirect = useCallback(() => {
    const { renderer, scene, camera } = threeRefs.current;
    if (renderer && scene && camera) {
      renderer.render(scene, camera);
    }
  }, []);

  // update visualisation when showPhase changes
  useEffect(() => {
    const uShow = threeRefs.current?.points?.material?.uniforms?.uShowPhase;
    if (uShow) {
      uShow.value = showPhase ? 1 : 0;
    }
    updateVisualisation();
    renderOnce();
  }, [showPhase, updateVisualisation, renderOnce]);

  // handle window resize
  useEffect(() => {
    const onWindowResize = () => {
      const { renderer, camera, points } = threeRefs.current;
      handleResize(renderer, camera, mountRef.current, points);
      renderOnce();
    };

    globalThis.window?.addEventListener('resize', onWindowResize);

    return () => {
      globalThis.window?.removeEventListener('resize', onWindowResize);
      disposeThreeJS(threeRefs.current);
    };
  }, [renderOnce]);


  useEffect(() => {
    if (isRecording) {
      setRecordElapsedSec(0);
      const id = setInterval(() => setRecordElapsedSec((v) => v + 1), 1000);
      recordTimerRef.current = id;
      return () => {
        clearInterval(id);
        recordTimerRef.current = null;
      };
    }
    return undefined;
  }, [isRecording]);


  // animation loop
  useEffect(() => {
    let raf = 0;

    const tick = () => {
      if (typeof document !== 'undefined' && document.hidden) {
        raf = 0;
        return;
      }

      const { controls } = threeRefs.current;
      controls?.update();

      if (running) {
        for (let s = 0; s < stepsPerFrame; s++) {
          timeStep(
            psiRe.current,
            psiIm.current,
            expVh.current,
            expK.current,
            N,
            scratchRe.current,
            scratchIm.current
          );
        }
        updateVisualisation();
        renderOnce();
        {
          const { scene, camera } = threeRefs.current;
          if (recorderRef.current && scene && camera) {
            recorderRef.current.renderFrame(scene, camera);
          }
        }
        raf = requestAnimationFrame(tick);
      }
    };

    const onControlsChange = () => {
      if (!running) renderDirect();
      {
        const { scene, camera } = threeRefs.current;
        if (recorderRef.current && scene && camera) {
          recorderRef.current.renderFrame(scene, camera);
        }
      }
    };

    const { controls } = threeRefs.current;
    controls?.addEventListener?.('change', onControlsChange);

    const onVis = () => {
      if (typeof document !== 'undefined' && !document.hidden && running && !raf) {
        raf = requestAnimationFrame(tick);
      }
    };

    if (typeof document !== 'undefined') {
      document.addEventListener('visibilitychange', onVis);
    }

    if (running) {
      raf = requestAnimationFrame(tick);
    } else {
      renderOnce();
    }

    return () => {
      cancelAnimationFrame(raf);
      controls?.removeEventListener?.('change', onControlsChange);
      if (typeof document !== 'undefined') {
        document.removeEventListener('visibilitychange', onVis);
      }
    };
  }, [running, stepsPerFrame, dt, densityScale, showPhase, N, updateVisualisation, renderOnce, renderDirect]);

  // potential preset functions
  function handlePresetFree() {
    presetFree(V.current);
    rebuildPotentialExponentials();
  }

  function handlePresetPlaneBarrier() {
    presetPlaneBarrier(V.current, coord, N);
    rebuildPotentialExponentials();
  }

  function handlePresetBoxWell() {
    presetBoxWell(V.current, coord, N, L);
    rebuildPotentialExponentials();
  }

  function handlePresetSphere() {
    presetSphere(V.current, coord, N, L);
    rebuildPotentialExponentials();
  }

  function handlePresetHarmonic() {
    presetHarmonic(V.current, coord, N);
    rebuildPotentialExponentials();
  }

  function handleAddPacket() {
    const clampK = (v) => Math.max(-kMax, Math.min(kMax, v));
    const kkx = clampK(k0x);
    const kky = clampK(k0y);
    const kkz = clampK(k0z);
    addPacket3D(
      psiRe.current,
      psiIm.current,
      coord,
      N,
      x0, y0, z0,
      sigma, sigma, sigma,
      kkx, kky, kkz,
      amp,
      gXRef.current,
      gYRef.current,
      gZRef.current,
      pXRef.current,
      pYRef.current,
      pZRef.current
    );
    const cellVol = dx * dx * dx;
    renormalize(psiRe.current, psiIm.current, cellVol);
    updateVisualisation();
    renderOnce();
  }

  function handleResetPsi() {
    psiRe.current.fill(0);
    psiIm.current.fill(0);
    maxDRef.current = 1e-9;
    updateVisualisation();
    renderOnce();
  }

  async function handleScreenshot() {
    const { scene, camera, points, renderer } = threeRefs.current;
    if (!scene || !camera) return;
    const { width, height } = getCaptureSize(3840);
    const ts = new Date();
    const pad = (n) => String(n).padStart(2, '0');
    const filename = `quantum_screenshot_${ts.getFullYear()}${pad(ts.getMonth() + 1)}${pad(ts.getDate())}_${pad(ts.getHours())}${pad(ts.getMinutes())}${pad(ts.getSeconds())}.png`;
    try {
      setIsCapturingShot(true);
      await captureScreenshot({ scene, camera, points, renderer, width, height, dpr: 1, ssaa: 2, filename });
    } catch (e) {
      console.error(e);
    } finally {
      setIsCapturingShot(false);
    }
  }

  function handleStartRecording() {
    if (recorderRef.current) return;
    const { scene, camera, points } = threeRefs.current;
    if (!scene || !camera) return;
    const { width, height } = getCaptureSize(1920);
    try {
      const session = startRecording({ scene, camera, points, width, height, dpr: 1, fps: 60, ssaa: 2, videoBitsPerSecond: 35000000 });
      recorderRef.current = session;
      setIsRecording(true);
    } catch (e) {
      console.error(e);
      recorderRef.current = null;
      setIsRecording(false);
    }
  }

  async function handleStopRecording() {
    const session = recorderRef.current;
    if (!session) return;
    try {
      const blob = await session.stop();
      recorderRef.current = null;
      setIsRecording(false);
      if (typeof document !== 'undefined') {
        const ts = new Date();
        const pad = (n) => String(n).padStart(2, '0');
        const filename = `quantum_recording_${ts.getFullYear()}${pad(ts.getMonth() + 1)}${pad(ts.getDate())}_${pad(ts.getHours())}${pad(ts.getMinutes())}${pad(ts.getSeconds())}.webm`;
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        setTimeout(() => {
          URL.revokeObjectURL(url);
          a.remove();
        }, 0);
      }
    } catch (e) {
      console.error(e);
      recorderRef.current = null;
      setIsRecording(false);
    }
  }

  // run self-tests on mount
  useEffect(() => {
    // avoid running FFT self-tests during unit/component test runs
    if (import.meta?.env?.MODE !== 'test') {
      runFFTSelfTests();
    }
  }, []);

  return (
    <div className="w-full max-w-6xl mx-auto p-4 md:p-6">
      <h1 className="text-2xl md:text-3xl font-semibold text-white tracking-tight mb-2">
        interactive quantum wave simulator: 3d (split-step fft)
      </h1>
      <p className="text-slate-300 mb-4">
        add a 3d gaussian packet, choose a potential, and orbit the volume. points encode |ψ|² by size and arg(ψ) by hue.
      </p>

      <div className="grid lg:grid-cols-3 gap-4 md:gap-6">
        <div className="lg:col-span-2 card p-4">
          <div 
            ref={mountRef} 
            className="w-full" 
            style={{ height: "min(60vh, 600px)" }} 
          />
          <div className="flex gap-2 mt-3 flex-wrap">
            <button 
              className={`btn ${running ? 'btn--primary' : 'btn--secondary'}`} 
              onClick={() => setRunning(r => !r)}
            >
              {running ? "pause" : "run"}
            </button>
            <button 
              className="btn btn--secondary" 
              onClick={handleResetPsi}
            >
              reset ψ
            </button>
            <button 
              className="btn btn--primary" 
              onClick={handleAddPacket}
            >
              add packet
            </button>
            <button 
              className="btn btn--ghost" 
              onClick={handlePresetFree}
            >
              free
            </button>
            <button 
              className="btn btn--ghost" 
              onClick={handlePresetPlaneBarrier}
            >
              plane barrier
            </button>
            <button 
              className="btn btn--ghost" 
              onClick={handlePresetBoxWell}
            >
              box well
            </button>
            <button 
              className="btn btn--ghost" 
              onClick={handlePresetSphere}
            >
              spherical well
            </button>
            <button 
              className="btn btn--ghost" 
              onClick={handlePresetHarmonic}
            >
              harmonic
            </button>
            <button 
              className="btn btn--secondary" 
              onClick={handleScreenshot}
              disabled={isCapturingShot}
            >
              {isCapturingShot ? "saving..." : "screenshot"}
            </button>
            <button 
              className="btn btn--primary" 
              onClick={handleStartRecording}
              disabled={isRecording}
            >
              {isRecording ? "recording..." : "start recording"}
            </button>
            <button 
              className="btn btn--secondary" 
              onClick={handleStopRecording}
              disabled={!isRecording}
            >
              stop recording
            </button>
            {isRecording && (
              <div className="flex items-center gap-2 ml-auto">
                <span className="w-2 h-2 bg-red-500 rounded-full inline-block" />
                <span className="text-red-500 font-medium">REC {recTimeLabel}</span>
              </div>
            )}
          </div>
        </div>

        <Controls
          N={N} setN={setN}
          L={L} setL={setL}
          dtScale={dtScale} setDtScale={setDtScale}
          dt={dt}
          kMax={kMax}
          stepsPerFrame={stepsPerFrame} setStepsPerFrame={setStepsPerFrame}
          absorbFrac={absorbFrac} setAbsorbFrac={setAbsorbFrac}
          absorbStrength={absorbStrength} setAbsorbStrength={setAbsorbStrength}
          sigma={sigma} setSigma={setSigma}
          k0x={k0x} setK0x={setK0x}
          k0y={k0y} setK0y={setK0y}
          k0z={k0z} setK0z={setK0z}
          amp={amp} setAmp={setAmp}
          x0={x0} setX0={setX0}
          y0={y0} setY0={setY0}
          z0={z0} setZ0={setZ0}
          densityScale={densityScale} setDensityScale={setDensityScale}
          showPhase={showPhase} setShowPhase={setShowPhase}
          running={running} setRunning={setRunning}
          onAddPacket={handleAddPacket}
          onResetPsi={handleResetPsi}
          onPresetFree={handlePresetFree}
          onPresetPlaneBarrier={handlePresetPlaneBarrier}
          onPresetBoxWell={handlePresetBoxWell}
          onPresetSphere={handlePresetSphere}
          onPresetHarmonic={handlePresetHarmonic}
        />
      </div>

      <div className="card p-3 md:p-4 text-slate-400 text-xs mt-4 leading-relaxed space-y-2">
        <p>
          tips: keep n=32 for smooth interactivity. use smaller <code>dt</code> for accuracy. the fft implies periodic boundaries. the cap mitigates wrap-around by damping outgoing flux. harmonic and spherical wells are good sanity checks (bound states, breathing modes).
        </p>
        <details className="mt-2">
          <summary className="cursor-pointer text-slate-300 hover:text-white font-medium">help: controls explained</summary>
          <div className="mt-2 text-slate-400 prose prose-sm prose-invert max-w-none">
            <p>mouse/touch navigation: left-drag orbits, wheel/trackpad zooms, right-drag pans. this only moves the camera, not the physics.</p>
            <p>points encode probability density and phase: size ∝ |ψ|, colour hue maps phase arg(ψ) from -π to π around the colour wheel.</p>
            <p><strong>grid n (32/64):</strong> samples along x, y, z (n³ total). larger n resolves finer structure but needs more memory/compute.</p>
            <p><strong>domain l:</strong> physical side length of the cube. spatial step is dx = l/n. increasing l makes the same packet look smaller.</p>
            <p><strong>dt scale:</strong> sets time step via dt = dtScale × dx². smaller dt reduces splitting error.</p>
            <p><strong>steps per frame:</strong> how many physics steps between renders.</p>
            <p><strong>periodic boundaries:</strong> the fft method is periodic. waves exiting one face re-enter at the opposite face unless absorbed (see cap).</p>
            <p><strong>σ (sigma):</strong> spatial width of gaussian packet. smaller σ spreads faster (fourier uncertainty trade-off).</p>
            <p><strong>centre x/y/z:</strong> where you drop the packet inside the box (clamped to ±l/2).</p>
            <p><strong>k0x/y/z:</strong> wave-vector components. set the direction and speed the blob travels.</p>
            <p><strong>amplitude:</strong> linear scale of the packet you add.</p>
            <p><strong>cap width:</strong> thickness of the absorbing sponge near the box faces. reduces wrap-around reflections.</p>
            <p><strong>cap strength:</strong> how strongly the cap damps outgoing flux each step. too high can nibble the interior.</p>
            <p><strong>potentials:</strong> free (v=0 everywhere), plane barrier (step on +x side), box well (cubic trap), spherical well (round trap), harmonic (spring to centre).</p>
            <p><strong>density scale:</strong> multiplies point size after auto-normalisation.</p>
            <p><strong>phase hue:</strong> toggles colour-by-phase.</p>
          </div>
        </details>
      </div>
    </div>
  );
}