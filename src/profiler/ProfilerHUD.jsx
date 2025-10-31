import React from 'react';
import { useProfilerEnabled, useProfilerSnapshot } from './store.js';

function formatF(f) { return Number.isFinite(f) ? f.toFixed(1) : '0.0'; }

export default function ProfilerHUD() {
  const enabled = useProfilerEnabled();
  const snap = useProfilerSnapshot();
  if (!enabled || !snap) return null;

  const { histogram, minFps, avgFps, maxFps, windowDropped, lifetimeDropped } = snap;
  const edges = histogram.edges;
  const counts = histogram.counts;
  const bins = counts.length;
  const width = 140;
  const height = 36;
  const pad = 2;
  const barW = Math.max(1, Math.floor((width - pad * 2) / bins));
  const maxCount = counts.reduce((a,b)=>Math.max(a,b), 0) || 1;

  return (
    <div className="absolute left-2 top-2 z-10 pointer-events-none select-none">
      <div className="px-2 py-1 rounded-md border text-xs font-medium" style={{ background: "rgba(0,0,0,0.4)", borderColor: "rgba(255,255,255,0.1)" }}>
        <div className="text-slate-100">fps {formatF(minFps)} / {formatF(avgFps)} / {formatF(maxFps)}</div>
        <svg width={width} height={height} className="mt-1 block">
          {counts.map((c, i) => {
            const h = Math.round((c / maxCount) * (height - pad * 2));
            const x = pad + i * barW;
            const y = height - pad - h;
            const k = `${edges[i]}-${edges[i+1]}`;
            return <rect key={k} x={x} y={y} width={barW - 1} height={h} fill="rgb(148,163,184)" />;
          })}
        </svg>
        <div className="mt-1 text-[10px] text-slate-300">dropped {windowDropped} · total {lifetimeDropped}</div>
      </div>
    </div>
  );
}
