// Three.js visualization setup and management
// Scene creation, point cloud rendering, and camera controls

import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { createQuantumShaderMaterial } from './shaders.js';

/**
 * Initialize Three.js scene, camera, renderer, and controls
 * @param {HTMLElement} mountElement - DOM element to mount renderer
 * @param {number} L - Domain size for camera positioning
 * @returns {Object} Three.js objects and references
 */
export function initializeThreeJS(mountElement, L) {
  // Clean up any existing content
  if (mountElement) mountElement.innerHTML = "";

  // Create scene
  const scene = new THREE.Scene();
  // monochrome ui theme: medium grey background
  scene.background = new THREE.Color(0x2a2a2a);

  // Create camera
  const camera = new THREE.PerspectiveCamera(45, 1, 0.01, 1000);
  camera.position.set(0, 0, L * 1.8);

  // Create renderer
  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: false });
  const dpr = Math.min(
    (typeof window !== 'undefined' && window.devicePixelRatio) ? 
    window.devicePixelRatio : 1, 
    2
  );
  renderer.setPixelRatio(dpr);

  if (mountElement) {
    renderer.setSize(mountElement.clientWidth, mountElement.clientHeight || 480);
    mountElement.appendChild(renderer.domElement);
  }

  // Get GPU capabilities
  const gl = renderer.getContext();
  const range = gl && gl.getParameter ? 
    gl.getParameter(gl.ALIASED_POINT_SIZE_RANGE) : [1, 64];
  const maxPointSize = Array.isArray(range) || (range && range.length) ? 
    range[1] : 64;

  // Create orbit controls
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
 * Create point cloud geometry and material for quantum wave visualization
 * @param {Float32Array} coord - Coordinate array
 * @param {number} N - Grid size
 * @param {number} maxPointSize - Maximum point size
 * @param {boolean} showPhase - Whether to show phase coloring
 * @returns {Object} Geometry, material, and attribute references
 */
export function createQuantumPointCloud(coord, N, maxPointSize, showPhase = true) {
  const size = N * N * N;
  
  // Create position array
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

  // Create geometry
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
  
  const ampAttr = new THREE.BufferAttribute(amps, 1);
  const phaseAttr = new THREE.BufferAttribute(phases, 1);
  ampAttr.setUsage(THREE.DynamicDrawUsage);
  phaseAttr.setUsage(THREE.DynamicDrawUsage);
  
  geometry.setAttribute('aAmp', ampAttr);
  geometry.setAttribute('aPhase', phaseAttr);
  geometry.computeBoundingSphere();

  // Create material
  const material = createQuantumShaderMaterial(maxPointSize, showPhase);

  // Create points object
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
 * Create bounding box helper for the simulation domain
 * @param {number} L - Domain size
 * @returns {THREE.Box3Helper} Box helper object
 */
export function createBoundingBox(L) {
  const box = new THREE.Box3(
    new THREE.Vector3(-L/2, -L/2, -L/2), 
    new THREE.Vector3(L/2, L/2, L/2)
  );
  // Neutral gray box lines for monochrome aesthetic
  return new THREE.Box3Helper(box, new THREE.Color(0x5a5a5a));
}

/**
 * Update point cloud visualization with current wave function data
 * @param {Float32Array} psiRe - Real part of wave function
 * @param {Float32Array} psiIm - Imaginary part of wave function
 * @param {THREE.BufferAttribute} ampAttr - Amplitude attribute
 * @param {THREE.BufferAttribute} phaseAttr - Phase attribute
 * @param {Object} maxDRef - Reference to maximum density tracker
 * @param {number} densityScale - Density scaling factor
 * @param {boolean} showPhase - Whether to update phase data
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
    
    const amp = dens > eps ? Math.sqrt(dens) * invNorm : 0.0;
    amps[i] = amp;
    
    if (showPhase && amp > 0.0) {
      phasesArr[i] = Math.atan2(pi, pr);
    }
  }
  
  maxDRef.current = runningMax;

  ampAttr.needsUpdate = true;
  if (showPhase) phaseAttr.needsUpdate = true;
}

/**
 * Handle window resize events
 * @param {THREE.WebGLRenderer} renderer - Three.js renderer
 * @param {THREE.PerspectiveCamera} camera - Three.js camera
 * @param {HTMLElement} mountElement - Mount element
 * @param {THREE.Points} points - Points object for shader uniform updates
 */
export function handleResize(renderer, camera, mountElement, points) {
  if (!mountElement || !renderer || !camera) return;
  
  const w = mountElement.clientWidth;
  const h = mountElement.clientHeight || 480;
  const ndpr = Math.min(
    (typeof window !== 'undefined' && window.devicePixelRatio) ? 
    window.devicePixelRatio : 1, 
    2
  );
  
  if (renderer.getPixelRatio && renderer.getPixelRatio() !== ndpr) {
    renderer.setPixelRatio(ndpr);
  }
  
  renderer.setSize(w, h);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
  
  // Update shader uniforms for point sizing
  if (points) {
    const material = points.material;
    const dpr = (renderer.getPixelRatio && renderer.getPixelRatio()) || 
      ((typeof window !== 'undefined' && window.devicePixelRatio) ? 
       window.devicePixelRatio : 1);
    const fovTan = Math.tan((camera.fov * Math.PI/180) * 0.5);
    
    if (material.uniforms && material.uniforms.uSizeScale) {
      material.uniforms.uSizeScale.value = (dpr / Math.max(1e-4, fovTan)) * 100.0;
    }
  }
}

/**
 * Render the scene once
 * @param {THREE.WebGLRenderer} renderer - Three.js renderer
 * @param {THREE.Scene} scene - Three.js scene
 * @param {THREE.PerspectiveCamera} camera - Three.js camera
 * @param {OrbitControls} controls - Orbit controls
 */
export function renderScene(renderer, scene, camera, controls) {
  controls?.update();
  if (renderer && scene && camera) {
    renderer.render(scene, camera);
  }
}

/**
 * Clean up Three.js resources
 * @param {Object} threeRefs - Object containing Three.js references
 */
export function disposeThreeJS(threeRefs) {
  const { controls, points, boxHelper, renderer, geometry, material } = threeRefs;
  
  if (controls) controls.dispose();
  
  if (points) {
    if (material && material.dispose) material.dispose();
    if (geometry) {
      // Dispose buffer attributes
      const attrs = geometry.attributes;
      Object.values(attrs).forEach(attr => {
        if (attr && attr.dispose) attr.dispose();
      });
      geometry.dispose();
    }
  }
  
  if (boxHelper) {
    try { 
      if (boxHelper.geometry && boxHelper.geometry.dispose) {
        boxHelper.geometry.dispose();
      }
    } catch(_) {}
    try { 
      if (boxHelper.material && boxHelper.material.dispose) {
        boxHelper.material.dispose();
      }
    } catch(_) {}
  }
  
  if (renderer) {
    renderer.dispose();
  }
}