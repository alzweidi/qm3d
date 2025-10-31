const NOW = () => (globalThis.performance?.now?.() ?? Date.now());

function clamp(v, min, max) { return Math.max(min, Math.min(max, v)); }

function makeRing(capacity) {
  return {
    capacity,
    times: new Float64Array(capacity),
    dts: new Float32Array(capacity),
    head: 0,
    count: 0,
    push(t, dt) {
      this.times[this.head] = t;
      this.dts[this.head] = dt;
      this.head = (this.head + 1) % this.capacity;
      if (this.count < this.capacity) this.count++;
    },
    tailIndex() {
      return (this.head - this.count + this.capacity) % this.capacity;
    },
    evictOlderThan(minT) {
      while (this.count > 0) {
        const tail = this.tailIndex();
        if (this.times[tail] >= minT) break;
        this.count--;
      }
    },
    forEach(cb) {
      let n = this.count;
      let idx = this.tailIndex();
      for (let i = 0; i < n; i++) {
        cb(this.times[idx], this.dts[idx]);
        idx = (idx + 1) % this.capacity;
      }
    },
    toArray() {
      const arr = [];
      this.forEach((t, dt) => arr.push({ t, dt }));
      return arr;
    },
  };
}

export function estimateVsyncMs(samples) {
  if (!samples.length) return 16.667;
  const sorted = samples.slice().sort((a,b)=>a-b);
  const q = Math.max(1, Math.floor(sorted.length * 0.25));
  const subset = sorted.slice(0, q);
  const med = subset[Math.floor(subset.length/2)] ?? sorted[0];
  return clamp(med, 8, 40);
}

export function computeDroppedFrames(dt, vsyncMs) {
  const v = vsyncMs || 16.667;
  const missed = Math.round((dt + v * 0.25) / v) - 1;
  return missed > 0 ? missed : 0;
}

export function buildHistogram(dtsMs, edges) {
  const counts = new Uint16Array(edges.length - 1);
  for (let i = 0; i < dtsMs.length; i++) {
    const x = dtsMs[i];
    let bin = edges.length - 2;
    for (let b = 0; b < edges.length - 1; b++) {
      if (x < edges[b+1]) { bin = b; break; }
    }
    counts[bin]++;
  }
  return counts;
}

export function createMonitor({ windowMs = 5000, notifyMs = 250, histogramEdgesMs } = {}) {
  const edges = histogramEdgesMs || [0,8,12,16.7,20,25,33.3,41,50,66,84,100,1e9];
  const capacity = Math.max(120, Math.ceil(windowMs / 4));
  const ring = makeRing(capacity);
  let rafId = 0;
  let running = false;
  let lastT = 0;
  let lastNotify = 0;
  let lifetimeDropped = 0;
  let lifetimeFrames = 0;
  let onSnapshot = null;

  const handleVis = () => {
    if (typeof document !== 'undefined' && document.hidden) {
      if (rafId) { cancelAnimationFrame(rafId); rafId = 0; }
      lastT = 0;
      // no notification needed here; HUD can infer idle from hidden
    } else if (running && !rafId) {
      lastT = 0;
      rafId = requestAnimationFrame(tick);
    }
  };

  const tick = (now) => {
    if (!running) { rafId = 0; return; }
    const t = Number.isFinite(now) ? now : NOW();
    if (lastT) {
      const dt = Math.max(0.1, t - lastT);
      ring.push(t, dt);
      ring.evictOlderThan(t - windowMs);
      lifetimeFrames++;
      const vsync = estimateVsyncMs(Array.from({length: ring.count}, (_, i) => {
        const idx = (ring.tailIndex() + i) % ring.capacity;
        return ring.dts[idx];
      }));
      lifetimeDropped += computeDroppedFrames(dt, vsync);
      if (!lastNotify || (t - lastNotify) >= notifyMs) {
        lastNotify = t;
        if (onSnapshot) {
          const dts = Array.from({length: ring.count}, (_, i) => {
            const idx = (ring.tailIndex() + i) % ring.capacity;
            return ring.dts[idx];
          });
          const fpsValues = dts.map(x => 1000 / x);
          const sum = dts.reduce((a,b)=>a+b,0);
          const avgDt = dts.length ? sum / dts.length : 0;
          const avgFps = avgDt ? 1000 / avgDt : 0;
          const minFps = fpsValues.length ? Math.min(...fpsValues) : 0;
          const maxFps = fpsValues.length ? Math.max(...fpsValues) : 0;
          const counts = buildHistogram(dts, edges);
          let windowDropped = 0;
          for (let i = 0; i < dts.length; i++) windowDropped += computeDroppedFrames(dts[i], vsync);
          onSnapshot({
            now: t,
            windowMs,
            samples: ring.count,
            lastDtMs: dts.length ? dts[dts.length - 1] : 0,
            avgFps,
            minFps,
            maxFps,
            histogram: { edges, counts: Array.from(counts) },
            vsyncMs: vsync,
            windowDropped,
            lifetimeDropped,
            lifetimeFrames,
            hidden: typeof document !== 'undefined' ? !!document.hidden : false,
          });
        }
      }
    }
    lastT = t;
    rafId = requestAnimationFrame(tick);
  };

  return {
    start() {
      if (running) return;
      running = true;
      lastT = 0;
      lastNotify = 0;
      if (typeof document !== 'undefined') {
        document.addEventListener('visibilitychange', handleVis);
      }
      rafId = requestAnimationFrame(tick);
    },
    stop() {
      running = false;
      if (rafId) { cancelAnimationFrame(rafId); rafId = 0; }
      if (typeof document !== 'undefined') {
        document.removeEventListener('visibilitychange', handleVis);
      }
      lastT = 0;
    },
    setOnSnapshot(cb) { onSnapshot = cb; },
  };
}
