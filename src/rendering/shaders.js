// shader definitions for quantum wave visualisation
// vertex and fragment shaders for point cloud rendering with phase colouring

import * as THREE from "three";

/**
 * vertex shader for quantum wave point cloud
 * handles amplitude-based sizing and position transformation
 */
export const vertexShader = `
  attribute float aAmp;
  attribute float aPhase;
  varying float vAmp;
  varying float vPhase;
  uniform float uSizeScale; // ≈ DPR / tan(fov/2) scaled
  uniform float uMaxPointSize;
  
  void main() {
    vAmp = aAmp;
    vPhase = aPhase;
    vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
    float dist = max(1.0, -mvPosition.z);
    float sizePx = vAmp * uSizeScale / dist; // size ∝ amplitude
    gl_PointSize = min(sizePx, uMaxPointSize);
    gl_Position = projectionMatrix * mvPosition;
  }
`;

/**
 * fragment shader for quantum wave point cloud
 * renders circular sprites with HSL phase colouring
 */
export const fragmentShader = `
  #ifdef GL_FRAGMENT_PRECISION_HIGH
  precision highp float;
  #else
  precision mediump float;
  #endif
  
  varying float vAmp;
  varying float vPhase;
  uniform float uShowPhase;

  // helper for HSL→RGB (pulled out to avoid nested function errors in WebGL1)
  float hsl_f(float n, float h, float a, float l) {
    float k = mod(n + h*12.0, 12.0);
    return l - a * max(-1.0, min(min(k - 3.0, 9.0 - k), 1.0));
  }
  
  vec3 hsl2rgb(float h, float s, float l) {
    float a = s * min(l, 1.0 - l);
    return vec3(hsl_f(0.0, h, a, l), hsl_f(8.0, h, a, l), hsl_f(4.0, h, a, l));
  }

  void main() {
    if (vAmp <= 0.0) discard; // early discard to reduce overdraw
    
    vec2 d = gl_PointCoord - vec2(0.5);
    float r = length(d);
    if (r > 0.5) discard; // round sprite
    
    float a = smoothstep(0.5, 0.0, r);
    
    // map phase [-π,π] → hue [0,1] (or fixed hue when disabled)
    float hue = (uShowPhase > 0.5) ? 
      ((vPhase / 3.14159265358979323846) * 0.5 + 0.5) : 
      0.6;
    
    vec3 color = hsl2rgb(hue, 1.0, 0.5);
    gl_FragColor = vec4(color, a);
  }
`;

/**
 * shader uniforms configuration
 */
export const shaderUniforms = {
  uSizeScale: { value: 120 },
  uShowPhase: { value: 1 },
  uMaxPointSize: { value: 64 }
};

/**
 * create shader material for quantum wave visualisation
 * @param {number} maxPointSize - maximum point size supported by GPU
 * @param {boolean} showPhase - whether to show phase colouring
 * @returns {THREE.ShaderMaterial} configured shader material
 */
export function createQuantumShaderMaterial(maxPointSize = 64, showPhase = true) {
  return new THREE.ShaderMaterial({
    transparent: true,
    depthWrite: false,
    depthTest: false,
    blending: THREE.AdditiveBlending,
    vertexColors: false,
    uniforms: {
      uSizeScale: { value: 120 },
      uShowPhase: { value: showPhase ? 1 : 0 },
      uMaxPointSize: { value: maxPointSize }
    },
    vertexShader,
    fragmentShader
  });
}