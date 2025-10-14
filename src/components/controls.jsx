// UI controls component for quantum wave simulation parameters
// handles all user interface elements and parameter management

import React from 'react';
import PropTypes from 'prop-types';

/**
 * reusable control component with label and input
 */
const Control = ({ label, children }) => (
  <div className="flex items-center justify-between gap-3 py-1.5">
    <div className="text-sm text-slate-200/90 w-48 font-medium">{label}</div>
    <div className="flex-1">{children}</div>
  </div>
);

Control.propTypes = {
  label: PropTypes.string.isRequired,
  children: PropTypes.node.isRequired,
};

/**
 * main controls panel component
 */
export default function Controls({
  // simulation parameters
  N, setN,
  L, setL,
  dtScale, setDtScale,
  dt,
  stepsPerFrame, setStepsPerFrame,
  
  // absorbing boundaries
  absorbFrac, setAbsorbFrac,
  absorbStrength, setAbsorbStrength,
  
  // wave packet parameters
  sigma, setSigma,
  k0x, setK0x,
  k0y, setK0y,
  k0z, setK0z,
  amp, setAmp,
  x0, setX0,
  y0, setY0,
  z0, setZ0,
  
  // visualisation
  densityScale, setDensityScale,
  showPhase, setShowPhase,
  
  // control actions
  running, setRunning,
  onAddPacket,
  onResetPsi,
  onPresetFree,
  onPresetPlaneBarrier,
  onPresetBoxWell,
  onPresetSphere,
  onPresetHarmonic
}) {
  return (
    <div className="card p-4">
      <h2 className="section-title text-base mb-2">controls</h2>

      {/* simulation parameters */}
      <Control label={`grid n = ${N}`}>
        <input 
          type="range" 
          min={32} 
          max={64} 
          step={32} 
          value={N} 
          onChange={e => setN(parseInt(e.target.value))} 
          className="w-full"
        />
      </Control>
      
      <Control label={`domain l = ${L.toFixed(1)}`}>
        <input 
          type="range" 
          min={6} 
          max={16} 
          step={0.5} 
          value={L} 
          onChange={e => setL(parseFloat(e.target.value))} 
          className="w-full"
        />
      </Control>
      
      <Control label={`dt scale = ${dtScale.toFixed(3)} (dt≈${dt.toExponential(1)})`}>
        <input 
          type="range" 
          min={0.02} 
          max={0.15} 
          step={0.001} 
          value={dtScale} 
          onChange={e => setDtScale(parseFloat(e.target.value))} 
          className="w-full"
        />
      </Control>
      
      <Control label={`Steps/frame = ${stepsPerFrame}`}>
        <input 
          type="range" 
          min={1} 
          max={4} 
          step={1} 
          value={stepsPerFrame} 
          onChange={e => setStepsPerFrame(parseInt(e.target.value))} 
          className="w-full"
        />
      </Control>

      <div className="divider my-3" />

      {/* wave packet parameters */}
      <h3 className="section-title text-sm mb-1">packet</h3>
      
      <Control label={`σ = ${sigma.toFixed(2)}`}>
        <input 
          type="range" 
          min={0.3} 
          max={1.5} 
          step={0.05} 
          value={sigma} 
          onChange={e => setSigma(parseFloat(e.target.value))} 
          className="w-full"
        />
      </Control>
      
      <Control label={`centre x = ${x0.toFixed(2)}`}>
        <input 
          type="range" 
          min={-L/2} 
          max={L/2} 
          step={0.1} 
          value={x0} 
          onChange={e => setX0(parseFloat(e.target.value))} 
          className="w-full"
        />
      </Control>
      
      <Control label={`centre y = ${y0.toFixed(2)}`}>
        <input 
          type="range" 
          min={-L/2} 
          max={L/2} 
          step={0.1} 
          value={y0} 
          onChange={e => setY0(parseFloat(e.target.value))} 
          className="w-full"
        />
      </Control>
      
      <Control label={`centre z = ${z0.toFixed(2)}`}>
        <input 
          type="range" 
          min={-L/2} 
          max={L/2} 
          step={0.1} 
          value={z0} 
          onChange={e => setZ0(parseFloat(e.target.value))} 
          className="w-full"
        />
      </Control>
      
      <Control label={`k0x = ${k0x.toFixed(1)}`}>
        <input 
          type="range" 
          min={-20} 
          max={20} 
          step={0.5} 
          value={k0x} 
          onChange={e => setK0x(parseFloat(e.target.value))} 
          className="w-full"
        />
      </Control>
      
      <Control label={`k0y = ${k0y.toFixed(1)}`}>
        <input 
          type="range" 
          min={-20} 
          max={20} 
          step={0.5} 
          value={k0y} 
          onChange={e => setK0y(parseFloat(e.target.value))} 
          className="w-full"
        />
      </Control>
      
      <Control label={`k0z = ${k0z.toFixed(1)}`}>
        <input 
          type="range" 
          min={-20} 
          max={20} 
          step={0.5} 
          value={k0z} 
          onChange={e => setK0z(parseFloat(e.target.value))} 
          className="w-full"
        />
      </Control>
      
      <Control label={`amplitude = ${amp.toFixed(2)}`}>
        <input 
          type="range" 
          min={0.5} 
          max={2.0} 
          step={0.1} 
          value={amp} 
          onChange={e => setAmp(parseFloat(e.target.value))} 
          className="w-full"
        />
      </Control>

      <div className="divider my-3" />

      {/* absorbing boundaries */}
      <h3 className="section-title text-sm mb-1">absorbing boundaries (cap)</h3>
      
      <Control label={`width = ${(absorbFrac*100).toFixed(0)}%`}>
        <input 
          type="range" 
          min={0.05} 
          max={0.3} 
          step={0.01} 
          value={absorbFrac} 
          onChange={e => setAbsorbFrac(parseFloat(e.target.value))} 
          className="w-full"
        />
      </Control>
      
      <Control label={`strength = ${absorbStrength.toFixed(2)}`}>
        <input 
          type="range" 
          min={1} 
          max={8} 
          step={0.1} 
          value={absorbStrength} 
          onChange={e => setAbsorbStrength(parseFloat(e.target.value))} 
          className="w-full"
        />
      </Control>

      <div className="divider my-3" />

      {/* visualisation */}
      <h3 className="section-title text-sm mb-1">visualisation</h3>
      
      <Control label={`density scale = ${densityScale.toFixed(2)}`}>
        <input 
          type="range" 
          min={0.4} 
          max={2.0} 
          step={0.05} 
          value={densityScale} 
          onChange={e => setDensityScale(parseFloat(e.target.value))} 
          className="w-full"
        />
      </Control>
      
      <Control label={`phase hue = ${showPhase ? "on" : "off"}`}>
        <input 
          type="checkbox" 
          checked={showPhase} 
          onChange={e => setShowPhase(e.target.checked)} 
        />
      </Control>

      <div className="divider my-3" />

      {/* action buttons */}
      <div className="space-y-2">
        <div className="flex gap-2 flex-wrap">
          <button 
            className="btn btn--secondary" 
            onClick={() => setRunning(r => !r)}
          >
            {running ? "pause" : "run"}
          </button>
          
          <button 
            className="btn btn--secondary" 
            onClick={onResetPsi}
          >
            reset ψ
          </button>
          
          <button 
            className="btn btn--primary" 
            onClick={onAddPacket}
          >
            add packet
          </button>
        </div>
        
        <div className="flex gap-2 flex-wrap">
          <button 
            className="btn btn--ghost" 
            onClick={onPresetFree}
          >
            free
          </button>
          
          <button 
            className="btn btn--ghost" 
            onClick={onPresetPlaneBarrier}
          >
            plane barrier
          </button>
          
          <button 
            className="btn btn--ghost" 
            onClick={onPresetBoxWell}
          >
            box well
          </button>
          
          <button 
            className="btn btn--ghost" 
            onClick={onPresetSphere}
          >
            spherical well
          </button>
          
          <button 
            className="btn btn--ghost" 
            onClick={onPresetHarmonic}
          >
            harmonic
          </button>
        </div>
      </div>
    </div>
  );
}

Controls.propTypes = {
  // simulation parameters
  N: PropTypes.number.isRequired,
  setN: PropTypes.func.isRequired,
  L: PropTypes.number.isRequired,
  setL: PropTypes.func.isRequired,
  dtScale: PropTypes.number.isRequired,
  setDtScale: PropTypes.func.isRequired,
  dt: PropTypes.number.isRequired,
  stepsPerFrame: PropTypes.number.isRequired,
  setStepsPerFrame: PropTypes.func.isRequired,

  // absorbing boundaries
  absorbFrac: PropTypes.number.isRequired,
  setAbsorbFrac: PropTypes.func.isRequired,
  absorbStrength: PropTypes.number.isRequired,
  setAbsorbStrength: PropTypes.func.isRequired,

  // wave packet parameters
  sigma: PropTypes.number.isRequired,
  setSigma: PropTypes.func.isRequired,
  k0x: PropTypes.number.isRequired,
  setK0x: PropTypes.func.isRequired,
  k0y: PropTypes.number.isRequired,
  setK0y: PropTypes.func.isRequired,
  k0z: PropTypes.number.isRequired,
  setK0z: PropTypes.func.isRequired,
  amp: PropTypes.number.isRequired,
  setAmp: PropTypes.func.isRequired,
  x0: PropTypes.number.isRequired,
  setX0: PropTypes.func.isRequired,
  y0: PropTypes.number.isRequired,
  setY0: PropTypes.func.isRequired,
  z0: PropTypes.number.isRequired,
  setZ0: PropTypes.func.isRequired,

  // visualisation
  densityScale: PropTypes.number.isRequired,
  setDensityScale: PropTypes.func.isRequired,
};