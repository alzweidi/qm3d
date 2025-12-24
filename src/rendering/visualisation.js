// three.js visualisation setup and management
// scene creation, point cloud rendering, and camera controls

import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { createQuantumShaderMaterial } from './shaders.js';

/**
 * initialise three.js scene, camera, renderer, and controls
 * @param {HTMLElement} mountElement - DOM element to mount renderer
 * @param {number} L - domain size for camera positioning
 * @returns {Object} three.js objects and references
 */
export function initialiseThreeJS(mountElement, L) {
  // clean up any existing content
  if (mountElement) mountElement.replaceChildren();

  // create scene
  const scene = new THREE.Scene();
  // monochrome ui theme: medium grey background
  scene.background = new THREE.Color(0x2a2a2a);

  // create camera
  const camera = new THREE.PerspectiveCamera(45, 1, 0.01, 1000);
  camera.position.set(0, 0, L * 1.8);

  // create renderer
  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: false });
  const dpr = Math.min(
    globalThis.window?.devicePixelRatio ?? 1,
    2
  );
  renderer.setPixelRatio(dpr);

  if (mountElement) {
    renderer.setSize(mountElement.clientWidth, mountElement.clientHeight || 480);
    mountElement.appendChild(renderer.domElement);
  }

  // get GPU capabilities
  const gl = renderer.getContext();
  const range = gl?.getParameter ?
    gl.getParameter(gl.ALIASED_POINT_SIZE_RANGE) : [1, 64];
  const maxPointSize = Array.isArray(range) && range.length > 1 ?
    range[1] : 64;

  // create orbit controls
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;

  return {
    scene,
    camera,
    renderer,
    controls,
    maxPointSize
  };
}

/**
 * create point cloud geometry and material for quantum wave visualisation
 * @param {Float32Array} coord - coordinate array
 * @param {number} N - grid size
 * @param {number} maxPointSize - maximum point size
 * @param {boolean} showPhase - whether to show phase colouring
 * @returns {Object} geometry, material, and attribute references
 */
export function createQuantumPointCloud(coord, N, maxPointSize, showPhase = true) {
  const size = N * N * N;
  
  // create position array
  const positions = new Float32Array(size * 3);
  const amps = new Float32Array(size);
  const phases = new Float32Array(size);
  
  let p = 0;
  for (let z = 0; z < N; z++) {
    for (let y = 0; y < N; y++) {
      for (let x = 0; x < N; x++) {
        positions[p++] = coord[x];
        positions[p++] = coord[y];
        positions[p++] = coord[z];
      }
    }
  }

  // create geometry
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
  
  const ampAttr = new THREE.BufferAttribute(amps, 1);
  const phaseAttr = new THREE.BufferAttribute(phases, 1);
  ampAttr.setUsage(THREE.DynamicDrawUsage);
  phaseAttr.setUsage(THREE.DynamicDrawUsage);
  
  geometry.setAttribute('aAmp', ampAttr);
  geometry.setAttribute('aPhase', phaseAttr);
  geometry.computeBoundingSphere();

  // create material
  const material = createQuantumShaderMaterial(maxPointSize, showPhase);

  // create points object
  const points = new THREE.Points(geometry, material);
  points.frustumCulled = false;

  return {
    points,
    geometry,
    material,
    ampAttr,
    phaseAttr
  };
}

/**
 * create bounding box helper for the simulation domain
 * @param {number} L - domain size
 * @returns {THREE.Box3Helper} box helper object
 */
export function createBoundingBox(L) {
  const box = new THREE.Box3(
    new THREE.Vector3(-L/2, -L/2, -L/2), 
    new THREE.Vector3(L/2, L/2, L/2)
  );
  // neutral grey box lines for monochrome aesthetic
  return new THREE.Box3Helper(box, new THREE.Color(0x5a5a5a));
}

/**
 * update point cloud visualisation with current wave function data
 * @param {Float32Array} psiRe - real part of wave function
 * @param {Float32Array} psiIm - imaginary part of wave function
 * @param {THREE.BufferAttribute} ampAttr - amplitude attribute
 * @param {THREE.BufferAttribute} phaseAttr - phase attribute
 * @param {Object} maxDRef - reference to maximum density tracker
 * @param {number} densityScale - density scaling factor
 * @param {boolean} showPhase - whether to update phase data
 */
export function updatePointCloud(psiRe, psiIm, ampAttr, phaseAttr, maxDRef, densityScale, showPhase) {
  const amps = ampAttr.array;
  const phasesArr = showPhase ? phaseAttr.array : null;

  const prevMax = Math.max(1e-20, maxDRef.current);
  const invNorm = densityScale / Math.sqrt(prevMax);
  const eps = 1e-4 * prevMax;

  let runningMax = maxDRef.current * 0.98;
  
  for (let i = 0; i < psiRe.length; i++) {
    const pr = psiRe[i], pi = psiIm[i];
    const dens = pr*pr + pi*pi;
    
    if (dens > runningMax) runningMax = dens;
    
    const amp = dens > eps ? Math.sqrt(dens) * invNorm : 0;
    amps[i] = amp;
    
    if (showPhase && amp > 0) {
      phasesArr[i] = Math.atan2(pi, pr);
    }
  }
  
  maxDRef.current = runningMax;

  ampAttr.needsUpdate = true;
  if (showPhase) phaseAttr.needsUpdate = true;
}

/**
 * handle window resize events
 * @param {THREE.WebGLRenderer} renderer - three.js renderer
 * @param {THREE.PerspectiveCamera} camera - three.js camera
 * @param {HTMLElement} mountElement - mount element
 * @param {THREE.Points} points - points object for shader uniform updates
 */
export function handleResize(renderer, camera, mountElement, points) {
  if (!mountElement || !renderer || !camera) return;
  
  const w = mountElement.clientWidth;
  const h = mountElement.clientHeight || 480;
  const ndpr = Math.min(
    globalThis.window?.devicePixelRatio ?? 1,
    2
  );
  
  const currentDpr = renderer.getPixelRatio?.();
  if (currentDpr !== undefined && currentDpr !== ndpr) {
    renderer.setPixelRatio(ndpr);
  }
  
  renderer.setSize(w, h);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
  
  // update shader uniforms for point sizing
  if (points) {
    const material = points.material;
    const dpr = renderer.getPixelRatio?.() ?? globalThis.window?.devicePixelRatio ?? 1;
    const fovTan = Math.tan((camera.fov * Math.PI/180) * 0.5);

    if (material?.uniforms?.uSizeScale) {
      material.uniforms.uSizeScale.value = (dpr / Math.max(1e-4, fovTan)) * 100;
    }
  }
}

/**
 * render the scene once
 * @param {THREE.WebGLRenderer} renderer - three.js renderer
 * @param {THREE.Scene} scene - three.js scene
 * @param {THREE.PerspectiveCamera} camera - three.js camera
 * @param {OrbitControls} controls - orbit controls
 */
export function renderScene(renderer, scene, camera, controls) {
  controls?.update();
  if (renderer && scene && camera) {
    renderer.render(scene, camera);
  }
}

function safeDispose(obj) {
  if (obj && typeof obj.dispose === 'function') {
    obj.dispose();
  }
}

function disposeBufferAttributes(geometry) {
  if (!geometry) return;
  const attrs = geometry.attributes;
  if (attrs) {
    for (const attr of Object.values(attrs)) {
      if (attr && typeof attr.dispose === 'function') attr.dispose();
    }
  }
  if (typeof geometry.dispose === 'function') geometry.dispose();
}

function disposePointsResources(points, geometry, material) {
  if (!points) return;
  safeDispose(material);
  disposeBufferAttributes(geometry);
}

function disposeBoxHelper(boxHelper) {
  if (!boxHelper) return;
  boxHelper.geometry?.dispose?.();
  boxHelper.material?.dispose?.();
}

/**
 * clean up three.js resources
 * @param {Object} threeRefs - object containing three.js references
 */
export function disposeThreeJS(threeRefs) {
  const { controls, points, boxHelper, renderer, geometry, material } = threeRefs;

  safeDispose(controls);
  disposePointsResources(points, geometry, material);
  disposeBoxHelper(boxHelper);
  safeDispose(renderer);
}