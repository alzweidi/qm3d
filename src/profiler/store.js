import { createMonitor } from './rafMonitor.js';
import React from 'react';

const KEY = 'hud:profilerEnabled';

function readPersisted() {
  try {
    if (import.meta?.env?.MODE === 'test') return null;
    if (typeof localStorage === 'undefined') return null;
    const v = localStorage.getItem(KEY);
    if (v == null) return null;
    return v === '1';
  } catch { return null; }
}

function writePersisted(enabled) {
  try {
    if (import.meta?.env?.MODE === 'test') return;
    if (typeof localStorage === 'undefined') return;
    localStorage.setItem(KEY, enabled ? '1' : '0');
  } catch { void 0; }
}

const envDefault = () => {
  if (import.meta?.env?.MODE === 'test') return false;
  const persisted = readPersisted();
  if (persisted != null) return persisted;
  const flag = import.meta?.env?.VITE_ENABLE_PROFILER;
  if (flag === 'true') return true;
  return false;
};

const monitor = createMonitor();

const state = {
  enabled: envDefault(),
  snapshot: null,
};

const listeners = new Set();

monitor.setOnSnapshot((snap) => {
  state.snapshot = snap;
  for (const l of listeners) l();
});

function startStop() {
  if (state.enabled) monitor.start(); else monitor.stop();
}

startStop();

export const ProfilerStore = {
  subscribe(listener) { listeners.add(listener); return () => listeners.delete(listener); },
  getSnapshot() { return state.snapshot; },
  getEnabled() { return state.enabled; },
  setEnabled(v) {
    const next = !!v;
    if (state.enabled === next) return;
    state.enabled = next;
    writePersisted(state.enabled);
    startStop();
    for (const l of listeners) l();
  },
  toggle() { this.setEnabled(!state.enabled); },
};

export function useProfilerSnapshot() {
  return React.useSyncExternalStore(ProfilerStore.subscribe, ProfilerStore.getSnapshot, ProfilerStore.getSnapshot);
}

export function useProfilerEnabled() {
  return React.useSyncExternalStore(ProfilerStore.subscribe, ProfilerStore.getEnabled, ProfilerStore.getEnabled);
}

export function useProfilerControls() {
  const enabled = useProfilerEnabled();
  const toggle = React.useCallback(() => ProfilerStore.toggle(), []);
  const setEnabled = React.useCallback((v) => ProfilerStore.setEnabled(!!v), []);
  return { enabled, toggle, setEnabled };
}
