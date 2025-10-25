import * as THREE from "three";

function createOffscreenRenderer(width, height, dpr = 1, preserve = false) {
  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: false, preserveDrawingBuffer: preserve });
  renderer.setPixelRatio(dpr);
  renderer.setSize(width, height, false);
  return renderer;
}

function getMaxPointSize(renderer) {
  const gl = renderer?.getContext?.();
  try {
    const range = gl?.getParameter?.(gl.ALIASED_POINT_SIZE_RANGE);
    if (Array.isArray(range) && range.length > 1) return range[1];
  } catch (e) { void e; }
  return 64;
}

function computeSizeScale(dpr, fovDeg) {
  const fovTan = Math.tan((fovDeg * Math.PI / 180) * 0.5);
  return (dpr / Math.max(1e-4, fovTan)) * 100;
}

function syncCameraToSize(srcCam, dstCam, width, height) {
  dstCam.fov = srcCam.fov;
  dstCam.near = srcCam.near;
  dstCam.far = srcCam.far;
  dstCam.up.copy(srcCam.up);
  dstCam.position.copy(srcCam.position);
  dstCam.quaternion.copy(srcCam.quaternion);
  dstCam.aspect = width / height;
  dstCam.updateProjectionMatrix();
}

function setUniformsForOffscreen(points, renderer, camera, maxPointSize) {
  const dpr = renderer.getPixelRatio?.() ?? 1;
  const mat = points?.material;
  const uniforms = mat?.uniforms;
  const prev = {
    sizeScale: uniforms?.uSizeScale?.value,
    maxSize: uniforms?.uMaxPointSize?.value,
  };
  if (uniforms?.uSizeScale) uniforms.uSizeScale.value = computeSizeScale(dpr, camera.fov);
  if (uniforms?.uMaxPointSize && typeof maxPointSize === 'number') uniforms.uMaxPointSize.value = maxPointSize;
  return prev;
}

function restoreUniforms(points, prev) {
  const uniforms = points?.material?.uniforms;
  if (!uniforms || !prev) return;
  if (uniforms.uSizeScale && typeof prev.sizeScale === 'number') uniforms.uSizeScale.value = prev.sizeScale;
  if (uniforms.uMaxPointSize && typeof prev.maxSize === 'number') uniforms.uMaxPointSize.value = prev.maxSize;
}

export async function captureScreenshot({ scene, camera, points, width, height, dpr = 1, filename }) {
  if (!scene || !camera) throw new Error('Scene and camera are required');
  const offscreen = createOffscreenRenderer(width, height, dpr, true);
  const camClone = new THREE.PerspectiveCamera();
  syncCameraToSize(camera, camClone, width, height);
  const maxPs = getMaxPointSize(offscreen);
  const prev = setUniformsForOffscreen(points, offscreen, camClone, maxPs);
  try {
    offscreen.render(scene, camClone);
  } finally {
    restoreUniforms(points, prev);
  }
  const canvas = offscreen.domElement;
  const blob = await new Promise((resolve, reject) => {
    try {
      canvas.toBlob((b) => {
        if (b) resolve(b); else reject(new Error('toBlob failed'));
      }, 'image/png');
    } catch (e) {
      reject(e);
    }
  });
  offscreen.dispose();
  if (filename && typeof document !== 'undefined') {
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
  return blob;
}

function selectMimeType(candidates) {
  if (typeof window === 'undefined' || typeof window.MediaRecorder === 'undefined') return null;
  const isSupported = window.MediaRecorder.isTypeSupported?.bind(window.MediaRecorder);
  for (const t of candidates) {
    if (!t) continue;
    if (!isSupported || isSupported(t)) return t;
  }
  return '';
}

export function startRecording({ scene, camera, points, width, height, dpr = 1, fps = 60, mimeCandidates = ['video/webm;codecs=vp9', 'video/webm;codecs=vp8', 'video/webm'], videoBitsPerSecond = 25000000 }) {
  if (!scene || !camera) throw new Error('Scene and camera are required');
  const mimeType = selectMimeType(mimeCandidates);
  if (mimeType === null) throw new Error('MediaRecorder is not available');
  if (mimeType === '') throw new Error('No supported MediaRecorder mime type');

  const offscreen = createOffscreenRenderer(width, height, dpr, false);
  const stream = offscreen.domElement.captureStream?.(fps);
  if (!stream) {
    offscreen.dispose();
    throw new Error('Canvas captureStream is not available');
  }

  const camClone = new THREE.PerspectiveCamera();
  syncCameraToSize(camera, camClone, width, height);
  const maxPs = getMaxPointSize(offscreen);

  const chunks = [];
  const recorder = new MediaRecorder(stream, { mimeType, videoBitsPerSecond });

  recorder.ondataavailable = (e) => {
    if (e?.data && e.data.size > 0) chunks.push(e.data);
  };

  let stopped = false;
  const stop = () => new Promise((resolve, reject) => {
    if (stopped) return resolve(new Blob([], { type: mimeType }));
    stopped = true;
    recorder.onstop = () => {
      try {
        const blob = new Blob(chunks, { type: mimeType });
        try { stream.getTracks?.().forEach(t => { try { t.stop(); } catch (e2) { void e2; } }); } catch (e2) { void e2; }
        offscreen.dispose();
        resolve(blob);
      } catch (e) {
        try { stream.getTracks?.().forEach(t => { try { t.stop(); } catch (e2) { void e2; } }); } catch (e2) { void e2; }
        offscreen.dispose();
        reject(e);
      }
    };
    try { recorder.stop(); } catch (e) { reject(e); }
  });

  const renderFrame = (srcScene, srcCamera) => {
    if (!srcScene || !srcCamera) return;
    syncCameraToSize(srcCamera, camClone, width, height);
    const prev = setUniformsForOffscreen(points, offscreen, camClone, maxPs);
    try { offscreen.render(srcScene, camClone); } finally { restoreUniforms(points, prev); }
  };

  try { recorder.start(); } catch (e) {
    try { stream.getTracks?.().forEach(t => { try { t.stop(); } catch (e2) { void e2; } }); } catch (e2) { void e2; }
    offscreen.dispose();
    throw e;
  }

  return { renderFrame, stop, mimeType };
}
