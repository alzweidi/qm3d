import React from 'react'
import ReactDOM from 'react-dom/client'
import './index.css'
import QuantumWaveEngine from './components/engine.jsx'
import ProfilerProvider from './profiler/provider.jsx'

ReactDOM.createRoot(document.getElementById('root')).render(
  <React.StrictMode>
    <ProfilerProvider>
      <QuantumWaveEngine />
    </ProfilerProvider>
  </React.StrictMode>,
)