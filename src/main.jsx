import React from 'react'
import ReactDOM from 'react-dom/client'
import './index.css'
import QuantumWaveSimulator3D from '../3d_quantum_wave_engine.jsx'

ReactDOM.createRoot(document.getElementById('root')).render(
  <React.StrictMode>
    <QuantumWaveSimulator3D />
  </React.StrictMode>,
)