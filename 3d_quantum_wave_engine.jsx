import React, { useEffect, useMemo, useRef, useState } from "react";
import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";

const cMul = (ar, ai, br, bi) => [ar*br - ai*bi, ar*bi + ai*br];
const cExp = (theta) => [Math.cos(theta), Math.sin(theta)];

function fft1d(re, im, inverse=false) {
  const n = re.length;
  if ((n & (n - 1)) !== 0) throw new Error('fft1d requires power‑of‑two length');
  let j = 0;
  for (let i = 0; i < n; i++) {
    if (i < j) { const tr = re[i]; const ti = im[i]; re[i] = re[j]; im[i] = im[j]; re[j] = tr; im[j] = ti; }
    let m = n >> 1;
    while (m >= 1 && j >= m) { j -= m; m >>= 1; }
    j += m;
  }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = (inverse ? +2*Math.PI : -2*Math.PI) / len;
    const wlen_r = Math.cos(ang), wlen_i = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let wr = 1, wi = 0;
      const half = len >> 1;
      for (let k = 0; k < half; k++) {
        const a = i + k, b = a + half;
        const u_r = re[a], u_i = im[a];
        const v_r = re[b]*wr - im[b]*wi;
        const v_i = re[b]*wi + im[b]*wr;
        re[a] = u_r + v_r; im[a] = u_i + v_i;
        re[b] = u_r - v_r; im[b] = u_i - v_i;
        const nwr = wr*wlen_r - wi*wlen_i; wi = wr*wlen_i + wi*wlen_r; wr = nwr;
      }
    }
  }
  if (inverse) { for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; } }
}

function makeFFTPlan(n){
  if ((n & (n-1)) !== 0) throw new Error('fft plan requires power‑of‑two length');
  const pairs = [];
  let j = 0;
  for (let i = 0; i < n; i++){
    if (i < j) pairs.push([i, j]);
    let m = n >> 1;
    while (m >= 1 && j >= m) { j -= m; m >>= 1; }
    j += m;
  }
  const stages = [];
  for (let len = 2; len <= n; len <<= 1){
    const angF = -2*Math.PI/len, angI = +2*Math.PI/len;
    stages.push({ len, half: len>>1, wlen_r_f: Math.cos(angF), wlen_i_f: Math.sin(angF), wlen_r_i: Math.cos(angI), wlen_i_i: Math.sin(angI) });
  }
  return { n, pairs, stages };
}

function fft1d_p(re, im, inverse, plan){
  const n = re.length;
  if (!plan || plan.n !== n) { fft1d(re, im, inverse); return; }
  for (let t=0; t<plan.pairs.length; t++){
    const a = plan.pairs[t][0], b = plan.pairs[t][1];
    const tr = re[a], ti = im[a]; re[a] = re[b]; im[a] = im[b]; re[b] = tr; im[b] = ti;
  }
  for (let s=0; s<plan.stages.length; s++){
    const st = plan.stages[s];
    const wlen_r = inverse ? st.wlen_r_i : st.wlen_r_f;
    const wlen_i = inverse ? st.wlen_i_i : st.wlen_i_f;
    const len = st.len, half = st.half;
    for (let i = 0; i < n; i += len){
      let wr = 1, wi = 0;
      for (let k = 0; k < half; k++){
        const a = i + k, b = a + half;
        const u_r = re[a], u_i = im[a];
        const v_r = re[b]*wr - im[b]*wi;
        const v_i = re[b]*wi + im[b]*wr;
        re[a] = u_r + v_r; im[a] = u_i + v_i;
        re[b] = u_r - v_r; im[b] = u_i - v_i;
        const nwr = wr*wlen_r - wi*wlen_i; wi = wr*wlen_i + wi*wlen_r; wr = nwr;
      }
    }
  }
  if (inverse) { for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; } }
}

function fft3d(re, im, N, inverse, lineRe, lineIm, plan){
  const idx = (x,y,z) => x + N*(y + N*z);
  for (let z=0; z<N; z++) for (let y=0; y<N; y++) {
    for (let x=0; x<N; x++){ const id=idx(x,y,z); lineRe[x]=re[id]; lineIm[x]=im[id]; }
    fft1d(lineRe, lineIm, inverse);
    for (let x=0; x<N; x++){ const id=idx(x,y,z); re[id]=lineRe[x]; im[id]=lineIm[x]; }
  }
  for (let z=0; z<N; z++) for (let x=0; x<N; x++) {
    for (let y=0; y<N; y++){ const id=idx(x,y,z); lineRe[y]=re[id]; lineIm[y]=im[id]; }
    fft1d(lineRe, lineIm, inverse);
    for (let y=0; y<N; y++){ const id=idx(x,y,z); re[id]=lineRe[y]; im[id]=lineIm[y]; }
  }
  for (let y=0; y<N; y++) for (let x=0; x<N; x++) {
    for (let z=0; z<N; z++){ const id=idx(x,y,z); lineRe[z]=re[id]; lineIm[z]=im[id]; }
    fft1d(lineRe, lineIm, inverse);
    for (let z=0; z<N; z++){ const id=idx(x,y,z); re[id]=lineRe[z]; im[id]=lineIm[z]; }
  }
}

export default function QuantumWaveSimulator3D(){
  const [N, setN] = useState(32);
  const [L, setL] = useState(10);
  const dx = useMemo(()=> L/(N), [L,N]);
  const [dtScale, setDtScale] = useState(0.08);
  const dt = dtScale*dx*dx;
  const cellVol = useMemo(()=> Math.pow(L/N, 3), [L, N]);
  const [stepsPerFrame, setStepsPerFrame] = useState(1);

  const [absorbFrac, setAbsorbFrac] = useState(0.15);
  const [absorbStrength, setAbsorbStrength] = useState(3.0);

  const [sigma, setSigma] = useState(0.6);
  const [k0x, setK0x] = useState(8.0);
  const [k0y, setK0y] = useState(0.0);
  const [k0z, setK0z] = useState(0.0);
  const [amp, setAmp] = useState(1.0);
  const [x0, setX0] = useState( -L/4 );
  const [y0, setY0] = useState( 0 );
  const [z0, setZ0] = useState( 0 );

  const [densityScale, setDensityScale] = useState(1.2);
  const [showPhase, setShowPhase] = useState(true);

  const mountRef = useRef(null);
  const rendererRef = useRef(null);
  const sceneRef = useRef(null);
  const cameraRef = useRef(null);
  const controlsRef = useRef(null);
  const pointsRef = useRef(null);
  const boxHelperRef = useRef(null);

  const [running, setRunning] = useState(true);

  const size = useMemo(()=>N*N*N, [N]);
  const psiRe = useRef(new Float32Array(size));
  const psiIm = useRef(new Float32Array(size));
  const V    = useRef(new Float32Array(size));
  const expK = useRef(new Float32Array(2*size));
  const expVh = useRef(new Float32Array(2*size));
  const capS2Ref = useRef(new Float32Array(size));

  const gXRef = useRef(new Float32Array(N)), gYRef = useRef(new Float32Array(N)), gZRef = useRef(new Float32Array(N));
  const pXRef = useRef(new Float32Array(N)), pYRef = useRef(new Float32Array(N)), pZRef = useRef(new Float32Array(N));

  const coord = useMemo(()=>{
    const arr = new Float32Array(N);
    const s = L/N; const off = -L/2 + s*0.5;
    for (let i=0;i<N;i++) arr[i] = off + i*s;
    return arr;
  }, [N, L]);

  const geomRef = useRef(null);
  const ampAttrRef = useRef(null);
  const phaseAttrRef = useRef(null);

  const maxDRef = useRef(1e-9);

  const scratchRe = useRef(new Float32Array(N));
  const scratchIm = useRef(new Float32Array(N));

  const idx = (x,y,z)=> x + N*(y + N*z);
  const kx2Ref = useRef(new Float32Array(N));
  const ky2Ref = useRef(new Float32Array(N));
  const kz2Ref = useRef(new Float32Array(N));

  useEffect(()=>{
    psiRe.current = new Float32Array(N*N*N);
    psiIm.current = new Float32Array(N*N*N);
    V.current    = new Float32Array(N*N*N);
    expK.current = new Float32Array(2*N*N*N);
    expVh.current= new Float32Array(2*N*N*N);

    gXRef.current = new Float32Array(N); gYRef.current = new Float32Array(N); gZRef.current = new Float32Array(N);
    pXRef.current = new Float32Array(N); pYRef.current = new Float32Array(N); pZRef.current = new Float32Array(N);

    scratchRe.current = new Float32Array(N);
    scratchIm.current = new Float32Array(N);

    V.current.fill(0);
    expK.current.fill(0); expVh.current.fill(0);
    psiRe.current.fill(0); psiIm.current.fill(0);
    maxDRef.current = 1e-9;

    initThree();
    rebuildExpVHalf();
  },[N, L]);

  useEffect(()=>{
    const clamp = (v)=> Math.max(-L/2, Math.min(L/2, v));
    setX0(x=>clamp(x));
    setY0(y=>clamp(y));
    setZ0(z=>clamp(z));
  }, [L]);

  useEffect(()=>{
    const cap = new Float32Array(N*N*N);
    const nAbs = Math.floor(absorbFrac * N);
    for (let z=0; z<N; z++){
      const zOff = N*N*z;
      const dz = Math.max(0, (nAbs - z)/nAbs, (z - (N-1-nAbs))/nAbs);
      for (let y=0; y<N; y++){
        const yOff = zOff + N*y;
        const dy = Math.max(0, (nAbs - y)/nAbs, (y - (N-1-nAbs))/nAbs);
        const dy2 = dy*dy;
        for (let x=0; x<N; x++){
          const dxv = Math.max(0, (nAbs - x)/nAbs, (x - (N-1-nAbs))/nAbs);
          const s2 = dxv*dxv + dy2 + dz*dz;
          cap[yOff + x] = s2;
        }
      }
    }
    capS2Ref.current = cap;
    rebuildExpVHalf();
  }, [N, absorbFrac]);

  useEffect(()=>{
    const kOf = (n)=>{ const m = (n <= N/2-1) ? n : n - N; return (2*Math.PI/L)*m; };
    for (let i=0;i<N;i++){
      const kx = kOf(i), ky = kOf(i), kz = kOf(i);
      kx2Ref.current[i] = kx*kx;
      ky2Ref.current[i] = ky*ky;
      kz2Ref.current[i] = kz*kz;
    }
  }, [N, L]);

  useEffect(()=>{
    const twoSize = 2*N*N*N;
    if (!expK.current || expK.current.length !== twoSize) expK.current = new Float32Array(twoSize);
    for (let z=0; z<N; z++){
      const zOff = N*N*z;
      for (let y=0; y<N; y++){
        const yOff = zOff + N*y;
        for (let x=0; x<N; x++){
          const id = yOff + x;
          const k2 = kx2Ref.current[x] + ky2Ref.current[y] + kz2Ref.current[z];
          const theta = -0.5 * k2 * dt;
          const cr = Math.cos(theta); const ci = Math.sin(theta);
          const j = id<<1; expK.current[j]=cr; expK.current[j+1]=ci;
        }
      }
    }
  },[N, L, dt]);

  useEffect(()=>{ rebuildExpVHalf(); }, [dt, N, absorbStrength, absorbFrac]);

  function rebuildExpVHalf(){
    const half = -0.5*dt;
    const len = N*N*N;
    if (!expVh.current || expVh.current.length !== 2*len) expVh.current = new Float32Array(2*len);
    const cap = capS2Ref.current;
    for (let i=0;i<len;i++){
      const theta = half * V.current[i];
      const decay = Math.exp(-0.5 * absorbStrength * (cap ? cap[i] : 0) * dt);
      const cr = Math.cos(theta) * decay; const ci = Math.sin(theta) * decay;
      const j = i<<1; expVh.current[j] = cr; expVh.current[j+1] = ci;
    }
  }

  function disposeThree(){
    if (controlsRef.current) controlsRef.current.dispose();
    if (pointsRef.current){
      const mat = pointsRef.current.material; if (mat && mat.dispose) mat.dispose();
      if (pointsRef.current.geometry) pointsRef.current.geometry.dispose();
      if (sceneRef.current) sceneRef.current.remove(pointsRef.current);
    }
    if (boxHelperRef.current){
      try { boxHelperRef.current.geometry && boxHelperRef.current.geometry.dispose && boxHelperRef.current.geometry.dispose(); } catch(_){ }
      try { boxHelperRef.current.material && boxHelperRef.current.material.dispose && boxHelperRef.current.material.dispose(); } catch(_){ }
      if (sceneRef.current) sceneRef.current.remove(boxHelperRef.current);
      boxHelperRef.current = null;
    }
    if (rendererRef.current){ rendererRef.current.dispose(); }
    rendererRef.current = null; pointsRef.current = null; geomRef.current = null; ampAttrRef.current = null; phaseAttrRef.current = null;
  }

  function initThree(){
    const mount = mountRef.current;
    disposeThree();
    if (mount) mount.innerHTML = "";

    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x0b1220);
    const camera = new THREE.PerspectiveCamera(45, 1, 0.01, 1000);
    camera.position.set(0, 0, L*1.8);

    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: false });
    const dpr = Math.min((typeof window !== 'undefined' && window.devicePixelRatio) ? window.devicePixelRatio : 1, 2);
    renderer.setPixelRatio(dpr);
    if (mount) {
      renderer.setSize(mount.clientWidth, mount.clientHeight || 480);
      mount.appendChild(renderer.domElement);
    }
    rendererRef.current = renderer;

    const gl = renderer.getContext();
    const range = gl && gl.getParameter ? gl.getParameter(gl.ALIASED_POINT_SIZE_RANGE) : [1, 64];
    const maxPoint = Array.isArray(range) || (range && range.length) ? range[1] : 64;

    const controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controlsRef.current = controls;

    const positions = new Float32Array(size*3);
    const amps      = new Float32Array(size);
    const phases    = new Float32Array(size);
    let p = 0;
    for (let z=0; z<N; z++) for (let y=0; y<N; y++) for (let x=0; x<N; x++){
      positions[p++] = coord[x];
      positions[p++] = coord[y];
      positions[p++] = coord[z];
    }

    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    const ampAttr   = new THREE.BufferAttribute(amps, 1);
    const phaseAttr = new THREE.BufferAttribute(phases, 1);
    ampAttr.setUsage(THREE.DynamicDrawUsage);
    phaseAttr.setUsage(THREE.DynamicDrawUsage);
    geometry.setAttribute('aAmp', ampAttr);
    geometry.setAttribute('aPhase', phaseAttr);
    geometry.computeBoundingSphere();

    const material = new THREE.ShaderMaterial({
      transparent: true,
      depthWrite: false,
      depthTest: false,
      blending: THREE.AdditiveBlending,
      vertexColors: false,
      uniforms: {
        uSizeScale:   { value: 120.0 },
        uShowPhase:   { value: showPhase ? 1.0 : 0.0 },
        uMaxPointSize:{ value: maxPoint || 64 }
      },
      vertexShader: `
        attribute float aAmp;
        attribute float aPhase;
        varying float vAmp;
        varying float vPhase;
        uniform float uSizeScale; // ≈ DPR / tan(fov/2) scaled
        uniform float uMaxPointSize;
        void main(){
          vAmp = aAmp;
          vPhase = aPhase;
          vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
          float dist = max(1.0, -mvPosition.z);
          float sizePx = vAmp * uSizeScale / dist; // size ∝ amplitude
          gl_PointSize = min(sizePx, uMaxPointSize);
          gl_Position = projectionMatrix * mvPosition;
        }
      `,
      fragmentShader: `
        #ifdef GL_FRAGMENT_PRECISION_HIGH
        precision highp float;
        #else
        precision mediump float;
        #endif
        varying float vAmp;
        varying float vPhase;
        uniform float uShowPhase;

        // Helper for HSL→RGB (pulled out to avoid nested function errors in WebGL1)
        float hsl_f(float n, float h, float a, float l){
          float k = mod(n + h*12.0, 12.0);
          return l - a * max(-1.0, min(min(k - 3.0, 9.0 - k), 1.0));
        }
        vec3 hsl2rgb(float h, float s, float l){
          float a = s * min(l, 1.0 - l);
          return vec3(hsl_f(0.0, h, a, l), hsl_f(8.0, h, a, l), hsl_f(4.0, h, a, l));
        }

        void main(){
          if (vAmp <= 0.0) discard; // early discard to reduce overdraw
          vec2 d = gl_PointCoord - vec2(0.5);
          float r = length(d);
          if (r > 0.5) discard; // round sprite
          float a = smoothstep(0.5, 0.0, r);
          // Map phase [-π,π] → hue [0,1] (or fixed hue when disabled)
          float hue = (uShowPhase > 0.5) ? ((vPhase / 3.14159265358979323846) * 0.5 + 0.5) : 0.6;
          vec3 color = hsl2rgb(hue, 1.0, 0.5);
          gl_FragColor = vec4(color, a);
        }
      `
    });

    const pts = new THREE.Points(geometry, material);
    pts.frustumCulled = false;
    scene.add(pts);

    const box = new THREE.Box3(new THREE.Vector3(-L/2,-L/2,-L/2), new THREE.Vector3(L/2,L/2,L/2));
    const helper = new THREE.Box3Helper(box, new THREE.Color(0x334155));
    scene.add(helper);
    boxHelperRef.current = helper;

    sceneRef.current = scene;
    cameraRef.current = camera;
    pointsRef.current = pts;
    geomRef.current = geometry;
    ampAttrRef.current = ampAttr;
    phaseAttrRef.current = phaseAttr;

    onResize();
  }

  function renderOnce(){
    controlsRef.current?.update();
    if (rendererRef.current && sceneRef.current && cameraRef.current)
      rendererRef.current.render(sceneRef.current, cameraRef.current);
  }
  useEffect(()=>{
    const pts = pointsRef.current; if (!pts) return;
    const mat = pts.material;
    if (mat && mat.uniforms && mat.uniforms.uShowPhase){ mat.uniforms.uShowPhase.value = showPhase ? 1.0 : 0.0; }
    updatePointCloud();
    renderOnce();
  }, [showPhase]);
  useEffect(()=>{
    const onWindowResize = () => onResize();
    if (typeof window !== 'undefined') window.addEventListener('resize', onWindowResize);
    return ()=> { if (typeof window !== 'undefined') window.removeEventListener('resize', onWindowResize); disposeThree(); };
  },[]);

  function onResize(){
    const mount = mountRef.current; if (!mount) return;
    const renderer = rendererRef.current, camera = cameraRef.current;
    if (!renderer || !camera) return;
    const w = mount.clientWidth; const h = mount.clientHeight || 480;
    const ndpr = Math.min((typeof window !== 'undefined' && window.devicePixelRatio) ? window.devicePixelRatio : 1, 2);
    if (renderer.getPixelRatio && renderer.getPixelRatio() !== ndpr) renderer.setPixelRatio(ndpr);
    renderer.setSize(w, h);
    camera.aspect = w / h; camera.updateProjectionMatrix();
    const pts = pointsRef.current; if (pts){
      const mat = pts.material;
      const dpr = (renderer.getPixelRatio && renderer.getPixelRatio()) || ((typeof window !== 'undefined' && window.devicePixelRatio) ? window.devicePixelRatio : 1);
      const fovTan = Math.tan((camera.fov * Math.PI/180) * 0.5);
      if (mat.uniforms && mat.uniforms.uSizeScale) mat.uniforms.uSizeScale.value = (dpr / Math.max(1e-4, fovTan)) * 100.0; // tuned factor
    }
    renderOnce();
  }

  function potentialHalfStep(){
    const ev = expVh.current; const len = psiRe.current.length;
    for (let i=0;i<len;i++){
      const pr = psiRe.current[i], pi = psiIm.current[i];
      const j = i<<1; const vr = ev[j], vi = ev[j+1];
      const nr = pr*vr - pi*vi; const ni = pr*vi + pi*vr;
      psiRe.current[i] = nr; psiIm.current[i] = ni;
    }
  }

  function kineticFullStep(){
    fft3d(psiRe.current, psiIm.current, N, false, scratchRe.current, scratchIm.current);
    const ek = expK.current; const len = psiRe.current.length;
    for (let i=0;i<len;i++){
      const j = i<<1; const kr = ek[j], ki = ek[j+1];
      const pr = psiRe.current[i], pi = psiIm.current[i];
      const nr = pr*kr - pi*ki; const ni = pr*ki + pi*kr;
      psiRe.current[i]=nr; psiIm.current[i]=ni;
    }
    fft3d(psiRe.current, psiIm.current, N, true, scratchRe.current, scratchIm.current);
  }

  function step(){
    potentialHalfStep();
    kineticFullStep();
    potentialHalfStep();
  }

  function addPacket3D(cx, cy, cz, sx, sy, sz, kx, ky, kz, scale=1){
    const gX = gXRef.current, gY = gYRef.current, gZ = gZRef.current;
    const pX = pXRef.current, pY = pYRef.current, pZ = pZRef.current;
    if (gX.length !== N) { gXRef.current = new Float32Array(N); gYRef.current = new Float32Array(N); gZRef.current = new Float32Array(N); }
    if (pX.length !== N) { pXRef.current = new Float32Array(N); pYRef.current = new Float32Array(N); pZRef.current = new Float32Array(N); }
    const s2x=sx*sx, s2y=sy*sy, s2z=sz*sz;
    for (let x=0;x<N;x++){ const X=coord[x]; gXRef.current[x]=Math.exp(-((X-cx)*(X-cx))/(2*s2x)); pXRef.current[x]=kx*(X-cx); }
    for (let y=0;y<N;y++){ const Y=coord[y]; gYRef.current[y]=Math.exp(-((Y-cy)*(Y-cy))/(2*s2y)); pYRef.current[y]=ky*(Y-cy); }
    for (let z=0;z<N;z++){ const Z=coord[z]; gZRef.current[z]=Math.exp(-((Z-cz)*(Z-cz))/(2*s2z)); pZRef.current[z]=kz*(Z-cz); }

    for (let z=0; z<N; z++){
      const zOff = N*N*z; const gz = gZRef.current[z]; const pz = pZRef.current[z];
      for (let y=0; y<N; y++){
        const yOff = zOff + N*y; const gy = gYRef.current[y]; const py = pYRef.current[y];
        for (let x=0; x<N; x++){
          const id = yOff + x;
          const g = scale * gXRef.current[x]*gy*gz;
          const phase = pXRef.current[x] + py + pz;
          psiRe.current[id] += g * Math.cos(phase);
          psiIm.current[id] += g * Math.sin(phase);
        }
      }
    }

    let norm = 0; const len = psiRe.current.length;
    for (let i=0;i<len;i++){ const pr=psiRe.current[i], pi=psiIm.current[i]; norm += pr*pr + pi*pi; }
    norm *= cellVol;
    if (norm < 1e-12) renormalize();
  }

  function renormalize(){
    let sum = 0; const len = psiRe.current.length;
    for (let i=0;i<len;i++){ const pr=psiRe.current[i], pi=psiIm.current[i]; sum += pr*pr + pi*pi; }
    const scale = 1/Math.sqrt(sum * cellVol + 1e-30);
    for (let i=0;i<len;i++){ psiRe.current[i]*=scale; psiIm.current[i]*=scale; }
  }

  function updatePointCloud(){
    const geom = geomRef.current, ampAttr = ampAttrRef.current, phaseAttr = phaseAttrRef.current;
    if (!geom || !ampAttr || !phaseAttr) return;
    const amps = ampAttr.array;
    const phasesArr = showPhase ? phaseAttr.array : null;

    const prevMax = Math.max(1e-20, maxDRef.current);
    const invNorm = densityScale / Math.sqrt(prevMax);
    const eps = 1e-4 * prevMax;

    let running = maxDRef.current * 0.98;
    for (let i=0;i<psiRe.current.length;i++){
      const pr = psiRe.current[i], pi = psiIm.current[i];
      const dens = pr*pr + pi*pi;
      if (dens > running) running = dens;
      const amp = dens > eps ? Math.sqrt(dens) * invNorm : 0.0;
      amps[i] = amp;
      if (showPhase && amp > 0.0) phasesArr[i] = Math.atan2(pi, pr);
    }
    maxDRef.current = running;

    ampAttr.needsUpdate = true;
    if (showPhase) phaseAttr.needsUpdate = true;
  }

  useEffect(()=>{
    let raf = 0;

    const tick = () => {
      if (typeof document !== 'undefined' && document.hidden) { raf = 0; return; }
      controlsRef.current?.update(); 
      if (running){
        for (let s=0;s<stepsPerFrame;s++) step();
        updatePointCloud();
        renderOnce();
        raf = requestAnimationFrame(tick);
      }
    };

    const controls = controlsRef.current;
    const onControlsChange = () => { if (!running) renderOnce(); };
    if (controls && controls.addEventListener) controls.addEventListener('change', onControlsChange);

    const onVis = () => { if (typeof document !== 'undefined' && !document.hidden && running && !raf) { raf = requestAnimationFrame(tick); } };
    if (typeof document !== 'undefined') document.addEventListener('visibilitychange', onVis);

    if (running) { raf = requestAnimationFrame(tick); } else { renderOnce(); }

    return ()=> {
      cancelAnimationFrame(raf);
      if (controls && controls.removeEventListener) controls.removeEventListener('change', onControlsChange);
      if (typeof document !== 'undefined') document.removeEventListener('visibilitychange', onVis);
    };
  },[running, stepsPerFrame, dt, densityScale, showPhase]);

  function presetFree(){ V.current.fill(0); rebuildExpVHalf(); }
  function presetPlaneBarrier(V0=4){
    for (let z=0; z<N; z++){
      const zOff = N*N*z;
      for (let y=0; y<N; y++){
        const yOff = zOff + N*y;
        for (let x=0; x<N; x++){
          V.current[yOff + x] = ( (coord[x]) > 0 ? V0 : 0);
        }
      }
    }
    rebuildExpVHalf();
  }
  function presetBoxWell(V0=-6){
    const w = L*0.4;
    for (let z=0; z<N; z++){
      const zOff = N*N*z;
      for (let y=0; y<N; y++){
        const yOff = zOff + N*y;
        for (let x=0; x<N; x++){
          const inside = Math.abs(coord[x])<w/2 && Math.abs(coord[y])<w/2 && Math.abs(coord[z])<w/2;
          V.current[yOff + x] = inside ? V0 : 0;
        }
      }
    }
    rebuildExpVHalf();
  }
  function presetSphere(V0=-8){
    const R = L*0.25; const R2 = R*R;
    for (let z=0; z<N; z++){
      const zOff = N*N*z; const z2 = coord[z]*coord[z];
      for (let y=0; y<N; y++){
        const yOff = zOff + N*y; const y2 = coord[y]*coord[y];
        for (let x=0; x<N; x++){
          const r2 = coord[x]*coord[x] + y2 + z2;
          V.current[yOff + x] = (r2 < R2 ? V0 : 0);
        }
      }
    }
    rebuildExpVHalf();
  }
  function presetHarmonic(omega=1){
    for (let z=0; z<N; z++){
      const zOff = N*N*z; const z2 = coord[z]*coord[z];
      for (let y=0; y<N; y++){
        const yOff = zOff + N*y; const y2 = coord[y]*coord[y];
        for (let x=0; x<N; x++){
          const r2 = coord[x]*coord[x] + y2 + z2;
          V.current[yOff + x] = 0.5*omega*omega*r2;
        }
      }
    }
    rebuildExpVHalf();
  }

  function handleAddPacket(){
    addPacket3D(x0, y0, z0, sigma, sigma, sigma, k0x, k0y, k0z, amp);
    updatePointCloud();
    renderOnce();
  }

  const Control = ({ label, children }) => (
    <div className="flex items-center justify-between gap-3 py-1">
      <div className="text-sm text-slate-300 w-48">{label}</div>
      <div className="flex-1">{children}</div>
    </div>
  );

  function runSelfTests(){
    try{
      const n = 32;
      const re = new Float32Array(n); const im = new Float32Array(n);
      for (let i=0;i<n;i++){ re[i] = Math.random()*2-1; im[i] = Math.random()*2-1; }
      const r1 = re.slice(), i1 = im.slice();
      fft1d(re, im, false); fft1d(re, im, true);
      let err = 0, norm = 0; for (let i=0;i<n;i++){ const dr = re[i]-r1[i], di = im[i]-i1[i]; err += dr*dr+di*di; norm += r1[i]*r1[i]+i1[i]*i1[i]; }
      console.assert(Math.sqrt(err/(norm+1e-12)) < 1e-5, "FFT(1D) round‑trip error too large");

      const cm = cMul(1,2,3,4); 
      console.assert(Math.abs(cm[0] + 5) < 1e-6 && Math.abs(cm[1] - 10) < 1e-6, "cMul incorrect");

      const dtt = 0.01, Vc = 2.0; const eh = cExp(-0.5*dtt*Vc);
      console.assert(Math.abs(eh[0]*eh[0] + eh[1]*eh[1] - 1) < 1e-6, "exp phase not unit magnitude");

      const NN = 8; const sz = NN*NN*NN;
      const Re3 = new Float32Array(sz); const Im3 = new Float32Array(sz);
      for (let i=0;i<sz;i++){ Re3[i] = Math.random()*2-1; Im3[i] = Math.random()*2-1; }
      const Re3b = Re3.slice(), Im3b = Im3.slice();
      const sRe = new Float32Array(NN), sIm = new Float32Array(NN);
      fft3d(Re3, Im3, NN, false, sRe, sIm);
      fft3d(Re3, Im3, NN, true,  sRe, sIm);
      let e3 = 0, n3 = 0; for (let i=0;i<sz;i++){ const dr=Re3[i]-Re3b[i], di=Im3[i]-Im3b[i]; e3 += dr*dr+di*di; n3 += Re3b[i]*Re3b[i]+Im3b[i]*Im3b[i]; }
      console.assert(Math.sqrt(e3/(n3+1e-12)) < 1e-5, "FFT(3D) round‑trip error too large");
    } catch(e){ console.warn("Self‑tests failed:", e); }
  }

  useEffect(()=>{ runSelfTests(); }, []);

  return (
    <div className="w-full max-w-6xl mx-auto p-4">
      <h1 className="text-2xl font-semibold text-white mb-2">Interactive Quantum Wave Simulator — 3D (Split–Step FFT)</h1>
      <p className="text-slate-300 mb-3">Add a 3D Gaussian packet, choose a potential, and orbit the volume. Points encode |ψ|² by size and arg(ψ) by hue.</p>

      <div className="grid lg:grid-cols-3 gap-4">
        <div className="lg:col-span-2 bg-slate-900/60 rounded-2xl p-3 shadow">
          <div ref={mountRef} className="w-full" style={{height: "min(60vh, 600px)"}} />
          <div className="flex gap-2 mt-3 flex-wrap">
            <button className="px-3 py-1.5 rounded-xl bg-emerald-600 hover:bg-emerald-500 text-white" onClick={()=>setRunning(r=>!r)}>
              {running?"Pause":"Run"}
            </button>
            <button className="px-3 py-1.5 rounded-xl bg-slate-700 hover:bg-slate-600 text-white" onClick={()=>{ psiRe.current.fill(0); psiIm.current.fill(0); maxDRef.current=1e-9; updatePointCloud(); renderOnce(); }}>
              Reset ψ
            </button>
            <button className="px-3 py-1.5 rounded-xl bg-indigo-600 hover:bg-indigo-500 text-white" onClick={handleAddPacket}>
              Add Packet
            </button>
            <button className="px-3 py-1.5 rounded-xl bg-slate-700 hover:bg-slate-600 text-white" onClick={presetFree}>Free</button>
            <button className="px-3 py-1.5 rounded-xl bg-slate-700 hover:bg-slate-600 text-white" onClick={()=>presetPlaneBarrier()}>Plane barrier</button>
            <button className="px-3 py-1.5 rounded-xl bg-slate-700 hover:bg-slate-600 text-white" onClick={()=>presetBoxWell()}>Box well</button>
            <button className="px-3 py-1.5 rounded-xl bg-slate-700 hover:bg-slate-600 text-white" onClick={()=>presetSphere()}>Spherical well</button>
            <button className="px-3 py-1.5 rounded-xl bg-slate-700 hover:bg-slate-600 text-white" onClick={()=>presetHarmonic()}>Harmonic</button>
          </div>
        </div>

        <div className="bg-slate-900/60 rounded-2xl p-4 shadow">
          <h2 className="text-lg font-medium text-white mb-2">Controls</h2>

          <Control label={`Grid N = ${N}`}>
            <input type="range" min={32} max={64} step={32} value={N} onChange={e=>setN(parseInt(e.target.value))} className="w-full"/>
          </Control>
          <Control label={`Domain L = ${L.toFixed(1)}`}>
            <input type="range" min={6} max={16} step={0.5} value={L} onChange={e=>setL(parseFloat(e.target.value))} className="w-full"/>
          </Control>
          <Control label={`dt scale = ${dtScale.toFixed(3)} (dt≈${dt.toExponential(1)})`}>
            <input type="range" min={0.02} max={0.15} step={0.001} value={dtScale} onChange={e=>setDtScale(parseFloat(e.target.value))} className="w-full"/>
          </Control>
          <Control label={`Steps/frame = ${stepsPerFrame}`}>
            <input type="range" min={1} max={4} step={1} value={stepsPerFrame} onChange={e=>setStepsPerFrame(parseInt(e.target.value))} className="w-full"/>
          </Control>

          <div className="h-px bg-slate-700 my-3" />

          <h3 className="text-white text-sm mb-1">Packet</h3>
          <Control label={`σ = ${sigma.toFixed(2)}`}>
            <input type="range" min={0.3} max={1.5} step={0.05} value={sigma} onChange={e=>setSigma(parseFloat(e.target.value))} className="w-full"/>
          </Control>
          <Control label={`Center x = ${x0.toFixed(2)}`}>
            <input type="range" min={-L/2} max={L/2} step={0.1} value={x0} onChange={e=>setX0(parseFloat(e.target.value))} className="w-full"/>
          </Control>
          <Control label={`Center y = ${y0.toFixed(2)}`}>
            <input type="range" min={-L/2} max={L/2} step={0.1} value={y0} onChange={e=>setY0(parseFloat(e.target.value))} className="w-full"/>
          </Control>
          <Control label={`Center z = ${z0.toFixed(2)}`}>
            <input type="range" min={-L/2} max={L/2} step={0.1} value={z0} onChange={e=>setZ0(parseFloat(e.target.value))} className="w-full"/>
          </Control>
          <Control label={`k0x = ${k0x.toFixed(1)}`}>
            <input type="range" min={-20} max={20} step={0.5} value={k0x} onChange={e=>setK0x(parseFloat(e.target.value))} className="w-full"/>
          </Control>
          <Control label={`k0y = ${k0y.toFixed(1)}`}>
            <input type="range" min={-20} max={20} step={0.5} value={k0y} onChange={e=>setK0y(parseFloat(e.target.value))} className="w-full"/>
          </Control>
          <Control label={`k0z = ${k0z.toFixed(1)}`}>
            <input type="range" min={-20} max={20} step={0.5} value={k0z} onChange={e=>setK0z(parseFloat(e.target.value))} className="w-full"/>
          </Control>
          <Control label={`Amplitude = ${amp.toFixed(2)}`}>
            <input type="range" min={0.5} max={2.0} step={0.1} value={amp} onChange={e=>setAmp(parseFloat(e.target.value))} className="w-full"/>
          </Control>

          <div className="h-px bg-slate-700 my-3" />

          <h3 className="text-white text-sm mb-1">Absorbing boundaries (CAP)</h3>
          <Control label={`Width = ${(absorbFrac*100).toFixed(0)}%`}>
            <input type="range" min={0.05} max={0.30} step={0.01} value={absorbFrac} onChange={e=>setAbsorbFrac(parseFloat(e.target.value))} className="w-full"/>
          </Control>
          <Control label={`Strength = ${absorbStrength.toFixed(2)}`}>
            <input type="range" min={1} max={8} step={0.1} value={absorbStrength} onChange={e=>setAbsorbStrength(parseFloat(e.target.value))} className="w-full"/>
          </Control>

          <div className="h-px bg-slate-700 my-3" />

          <h3 className="text-white text-sm mb-1">Rendering</h3>
          <Control label={`Density scale = ${densityScale.toFixed(2)}`}>
            <input type="range" min={0.4} max={2.0} step={0.05} value={densityScale} onChange={e=>setDensityScale(parseFloat(e.target.value))} className="w-full"/>
          </Control>
          <Control label={`Phase hue = ${showPhase?"on":"off"}`}>
            <input type="checkbox" checked={showPhase} onChange={e=>setShowPhase(e.target.checked)} />
          </Control>
        </div>
      </div>

      <div className="text-slate-400 text-xs mt-4 leading-relaxed">
        <p>
          tips: keep n=32 for smooth interactivity; use smaller <code>dt</code> for accuracy. The FFT implies periodic boundaries; the CAP
          mitigates wrap‑around by damping outgoing flux. harmonic and spherical wells are good sanity checks (bound states / breathing modes).
        </p>
      </div>
    </div>
  );
}
