import React from 'react';
import { ProfilerStore } from './store.js';

export default function ProfilerProvider({ children }) {
  React.useEffect(() => {
    const onKey = (e) => {
      if (e.defaultPrevented) return;
      if (e.metaKey || e.ctrlKey || e.altKey) return;
      const tag = e.target && e.target.tagName;
      if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT') return;
      if (e.isComposing) return;
      if (e.key === 'p' || e.key === 'P') {
        ProfilerStore.toggle();
      }
    };
    if (typeof window !== 'undefined') window.addEventListener('keydown', onKey);
    return () => { if (typeof window !== 'undefined') window.removeEventListener('keydown', onKey); };
  }, []);

  return children;
}
