import * as THREE from "three";

function createOffscreenRenderer(width, height, dpr = 1, preserve = false) {
  const renderer = new THREE.WebGLRenderer({
    antialias: true,
    alpha: false,
    preserveDrawingBuffer: preserve,
    powerPreference: 'high-performance'
  });

  if ("outputColorSpace" in renderer && THREE.SRGBColorSpace) {
    renderer.outputColorSpace = THREE.SRGBColorSpace;
  } else if ("outputEncoding" in renderer && THREE.sRGBEncoding) {
    renderer.outputEncoding = THREE.sRGBEncoding;
  }
  if ("NoToneMapping" in THREE) {
    renderer.toneMapping = THREE.NoToneMapping;
  }

  renderer.setPixelRatio(dpr);
  renderer.setSize(width, height, false);
  return renderer;
}

function getMaxPointSize(renderer) {
  const gl = renderer?.getContext?.();
  const range = gl?.getParameter?.(gl.ALIASED_POINT_SIZE_RANGE);
  if (Array.isArray(range) && range.length > 1) return range[1];
  return 64;
}

function computeSizeScale(dpr, fovDeg) {
  const fovTan = Math.tan((fovDeg * Math.PI / 180) * 0.5);
  return (dpr / Math.max(1e-4, fovTan)) * 100;
}

function syncCameraToSize(srcCam, dstCam, width, height) {
  dstCam.copy(srcCam, false);
  dstCam.aspect = width / height;
  if (typeof dstCam.clearViewOffset === 'function') dstCam.clearViewOffset();
  dstCam.updateProjectionMatrix();
  dstCam.updateMatrixWorld(true);
}

function setUniformsForOffscreen(points, renderer, camera, maxPointSize, dprOverride) {
  const dpr = (typeof dprOverride === 'number' && dprOverride > 0) ? dprOverride : (renderer.getPixelRatio?.() ?? 1);
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

function downscaleCanvas(srcCanvas, outWidth, outHeight) {
  const dst = document.createElement('canvas');
  dst.width = outWidth;
  dst.height = outHeight;
  const ctx = dst.getContext('2d', { willReadFrequently: true });
  if (ctx) {
    if ('imageSmoothingEnabled' in ctx) ctx.imageSmoothingEnabled = true;
    if ('imageSmoothingQuality' in ctx) ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(srcCanvas, 0, 0, srcCanvas.width, srcCanvas.height, 0, 0, outWidth, outHeight);
  }
  return dst;
}

function canvasToBlob(canvas) {
  return new Promise((resolve, reject) => {
    try {
      if (canvas.toBlob) {
        canvas.toBlob((b) => {
          if (b) {
            resolve(b);
            return;
          }
          try {
            const url = canvas.toDataURL('image/png');
            const arr = url.split(',');
            const mime = arr[0].match(/:(.*?);/)[1];
            const bstr = atob(arr[1]);
            let n = bstr.length;
            const u8arr = new Uint8Array(n);
            while (n--) u8arr[n] = (bstr.codePointAt(n) ?? 0) & 255;
            resolve(new Blob([u8arr], { type: mime }));
          } catch (error_) {
            reject(error_);
          }
        }, 'image/png');
      } else {
        const url = canvas.toDataURL('image/png');
        const arr = url.split(',');
        const mime = arr[0].match(/:(.*?);/)[1];
        const bstr = atob(arr[1]);
        let n = bstr.length;
        const u8arr = new Uint8Array(n);
        while (n--) u8arr[n] = (bstr.codePointAt(n) ?? 0) & 255;
        resolve(new Blob([u8arr], { type: mime }));
      }
    } catch (e) {
      reject(e);
    }
  });
}

export async function captureScreenshot({ scene, camera, points, renderer, width, height, dpr = 1, filename, ssaa = 1 }) {
  if (!scene || !camera) throw new Error('Scene and camera are required');
  if (!renderer) throw new Error('Renderer is required for screenshot');
  const ssaaFactor = Math.max(1, Math.floor(ssaa));
  const targetScale = (dpr || 1) * ssaaFactor;
  const camClone = new THREE.PerspectiveCamera();
  syncCameraToSize(camera, camClone, width, height);

  const gl = renderer.getContext?.();
  const maxRB = gl?.getParameter?.(gl.MAX_RENDERBUFFER_SIZE) ?? Infinity;
  const maxTex = gl?.getParameter?.(gl.MAX_TEXTURE_SIZE) ?? Infinity;
  const maxVD = gl?.getParameter?.(gl.MAX_VIEWPORT_DIMS);
  const maxView = Array.isArray(maxVD) ? Math.min(maxVD[0] || Infinity, maxVD[1] || Infinity) : Infinity;
  const limit = Math.max(1, Math.min(maxRB, maxTex, maxView));
  const baseSafe = Math.max(0.25, Math.min(targetScale, limit / Math.max(1, width), limit / Math.max(1, height)));
  const maxPixels = 16000000;
  const areaScale = Math.sqrt(Math.max(1, maxPixels) / Math.max(1, width * height));
  const safeScale = Math.min(baseSafe, areaScale);

  const rtW = Math.max(1, Math.ceil(width * safeScale));
  const rtH = Math.max(1, Math.ceil(height * safeScale));
  const rt = new THREE.WebGLRenderTarget(rtW, rtH, { depthBuffer: true, stencilBuffer: false });
  // ensure render target stores sRGB so readback matches on-screen appearance
  try {
    if (rt && rt.texture) {
      if (THREE?.SRGBColorSpace && 'colorSpace' in rt.texture) {
        rt.texture.colorSpace = THREE.SRGBColorSpace;
      } else if (THREE?.sRGBEncoding && 'encoding' in rt.texture) {
        // fallback for older three versions
        rt.texture.encoding = THREE.sRGBEncoding;
      }
    }
  } catch (_) { /* no-op */ }

  // Option A diagnostic: force pixel ratio = 1 for the main renderer during RT render
  const prevDpr = renderer.getPixelRatio?.();
  if (typeof renderer.setPixelRatio === 'function') renderer.setPixelRatio(1);

  const maxPs = getMaxPointSize(renderer);
  // match on-screen sprite sizing: use ORIGINAL screen DPR scaled by internal safeScale
  const screenDpr = (prevDpr !== undefined ? prevDpr : (globalThis.window?.devicePixelRatio ?? 1));
  const effectiveDpr = screenDpr * safeScale;
  const prevUniforms = setUniformsForOffscreen(points, renderer, camClone, maxPs, effectiveDpr);
  const prevTarget = renderer.getRenderTarget ? renderer.getRenderTarget() : null;
  const prevScissor = renderer.getScissor ? renderer.getScissor(new THREE.Vector4()) : null;
  const prevScissorTest = renderer.getScissorTest ? renderer.getScissorTest() : false;
  const prevViewport = renderer.getViewport ? renderer.getViewport(new THREE.Vector4()) : null;
  let outCanvas;
  try {
    renderer.setRenderTarget(rt);
    renderer.setScissorTest(false);
    renderer.setViewport(0, 0, rtW, rtH);
    renderer.clear(true, true, true);
    renderer.render(scene, camClone);
    const pixels = new Uint8Array(rtW * rtH * 4);
    renderer.readRenderTargetPixels(rt, 0, 0, rtW, rtH, pixels);
    const src = document.createElement('canvas');
    src.width = rtW; src.height = rtH;
    const sctx = src.getContext('2d', { willReadFrequently: true });
    const img = sctx.createImageData(rtW, rtH);
    const rowBytes = rtW * 4;
    for (let y = 0; y < rtH; y++) {
      const srcRow = (rtH - 1 - y) * rowBytes;
      const dstRow = y * rowBytes;
      img.data.set(pixels.subarray(srcRow, srcRow + rowBytes), dstRow);
    }
    sctx.putImageData(img, 0, 0);
    outCanvas = downscaleCanvas(src, width, height);
  } finally {
    restoreUniforms(points, prevUniforms);
    renderer.setRenderTarget(prevTarget || null);
    if (prevViewport) renderer.setViewport(prevViewport);
    if (prevScissor) renderer.setScissor(prevScissor);
    renderer.setScissorTest(prevScissorTest);
    // restore original pixel ratio
    if (prevDpr !== undefined && typeof renderer.setPixelRatio === 'function') {
      renderer.setPixelRatio(prevDpr);
    }
    if (typeof rt.dispose === 'function') rt.dispose();
  }
  const outBlob = await canvasToBlob(outCanvas);
  if (filename && typeof document !== 'undefined') {
    const url = URL.createObjectURL(outBlob);
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
  return outBlob;
}

function selectMimeType(candidates) {
  if (globalThis.MediaRecorder === undefined) return null;
  const isSupported = globalThis.MediaRecorder.isTypeSupported?.bind(globalThis.MediaRecorder);
  for (const t of candidates) {
    if (!t) continue;
    if (!isSupported || isSupported(t)) return t;
  }
  return '';
}

function stopTracks(stream) {
  const tracks = stream?.getTracks?.();
  if (Array.isArray(tracks)) {
    for (const t of tracks) {
      if (typeof t.stop === 'function') t.stop();
    }
  }
}

function createPresentCanvas(width, height) {
  const presentCanvas = document.createElement('canvas');
  presentCanvas.width = width;
  presentCanvas.height = height;
  const presentCtx = presentCanvas.getContext('2d');
  if (presentCtx) {
    if ('imageSmoothingEnabled' in presentCtx) presentCtx.imageSmoothingEnabled = true;
    if ('imageSmoothingQuality' in presentCtx) presentCtx.imageSmoothingQuality = 'high';
  }
  return { presentCanvas, presentCtx };
}

export function startRecording({ scene, camera, points, width, height, dpr = 1, fps = 60, mimeCandidates = ['video/webm;codecs=vp9', 'video/webm;codecs=vp8', 'video/webm'], videoBitsPerSecond = 25000000, ssaa = 1 }) {
  if (!scene || !camera) throw new Error('Scene and camera are required');
  const mimeType = selectMimeType(mimeCandidates);
  if (mimeType === null) throw new Error('MediaRecorder is not available');
  if (mimeType === '') throw new Error('No supported MediaRecorder mime type');

  const ssaaFactor = Math.max(1, Math.floor(ssaa));
  const composedDpr = (dpr || 1) * ssaaFactor;
  const offscreen = createOffscreenRenderer(width, height, composedDpr, false);

  let presentCanvas = null;
  let presentCtx = null;
  let stream = null;
  if (ssaaFactor > 1) {
    ({ presentCanvas, presentCtx } = createPresentCanvas(width, height));
    stream = presentCanvas.captureStream?.(fps);
  } else {
    stream = offscreen.domElement.captureStream?.(fps);
  }
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
      let blob;
      try { blob = new Blob(chunks, { type: mimeType }); }
      catch (e) { stopTracks(stream); offscreen.dispose(); return reject(e); }
      stopTracks(stream);
      offscreen.dispose();
      resolve(blob);
    };
    try { recorder.stop(); } catch (e) { stopTracks(stream); offscreen.dispose(); reject(e); }
  });

  const renderFrame = (srcScene, srcCamera) => {
    if (!srcScene || !srcCamera) return;
    syncCameraToSize(srcCamera, camClone, width, height);
    const prev = setUniformsForOffscreen(points, offscreen, camClone, maxPs);
    try {
      offscreen.render(srcScene, camClone);
      if (presentCtx && presentCanvas) {
        presentCtx.clearRect(0, 0, width, height);
        presentCtx.drawImage(offscreen.domElement, 0, 0, width, height);
      }
    } finally {
      restoreUniforms(points, prev);
    }
  };

  try { recorder.start(); } catch (e) { stopTracks(stream); offscreen.dispose(); throw e; }

  return { renderFrame, stop, mimeType };
}
