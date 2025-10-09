// @vitest-environment jsdom
import React from 'react';
import { describe, it, expect, vi } from 'vitest';
import { render, screen, fireEvent } from '@testing-library/react';
import Controls from '../../src/components/controls.jsx';

function renderControls(overrides = {}) {
  const props = {
    // sim params
    N: 32, setN: vi.fn(),
    L: 10, setL: vi.fn(),
    dtScale: 0.08, setDtScale: vi.fn(),
    dt: 0.08 * (10/32) * (10/32),
    stepsPerFrame: 1, setStepsPerFrame: vi.fn(),
    // absorbing
    absorbFrac: 0.15, setAbsorbFrac: vi.fn(),
    absorbStrength: 3.0, setAbsorbStrength: vi.fn(),
    // packet
    sigma: 0.6, setSigma: vi.fn(),
    k0x: 8.0, setK0x: vi.fn(),
    k0y: 0.0, setK0y: vi.fn(),
    k0z: 0.0, setK0z: vi.fn(),
    amp: 1.0, setAmp: vi.fn(),
    x0: -5, setX0: vi.fn(),
    y0: 0, setY0: vi.fn(),
    z0: 0, setZ0: vi.fn(),
    // viz
    densityScale: 1.2, setDensityScale: vi.fn(),
    showPhase: true, setShowPhase: vi.fn(),
    // actions
    running: true, setRunning: vi.fn(),
    onAddPacket: vi.fn(),
    onResetPsi: vi.fn(),
    onPresetFree: vi.fn(),
    onPresetPlaneBarrier: vi.fn(),
    onPresetBoxWell: vi.fn(),
    onPresetSphere: vi.fn(),
    onPresetHarmonic: vi.fn(),
    ...overrides,
  };
  render(<Controls {...props} />);
  return props;
}

describe('Controls component', () => {
  it('invokes handlers for buttons', () => {
    const p = renderControls();
    fireEvent.click(screen.getByText(/add packet/i));
    expect(p.onAddPacket).toHaveBeenCalledTimes(1);
    fireEvent.click(screen.getByText(/reset ψ/i));
    expect(p.onResetPsi).toHaveBeenCalledTimes(1);
    fireEvent.click(screen.getByText(/^free$/i));
    expect(p.onPresetFree).toHaveBeenCalledTimes(1);
    fireEvent.click(screen.getByText(/plane barrier/i));
    expect(p.onPresetPlaneBarrier).toHaveBeenCalledTimes(1);
    fireEvent.click(screen.getByText(/box well/i));
    expect(p.onPresetBoxWell).toHaveBeenCalledTimes(1);
    fireEvent.click(screen.getByText(/spherical well/i));
    expect(p.onPresetSphere).toHaveBeenCalledTimes(1);
    fireEvent.click(screen.getByText(/harmonic/i));
    expect(p.onPresetHarmonic).toHaveBeenCalledTimes(1);
  });

  it('toggles showPhase checkbox', () => {
    const p = renderControls({ showPhase: false });
    const label = screen.getByText(/phase hue/i);
    // The checkbox lives in the sibling container of the label within the same row
    const row = label.closest('div');
    const inputContainer = row?.nextElementSibling;
    const checkbox = inputContainer?.querySelector('input[type="checkbox"]');
    expect(checkbox).toBeTruthy();
    fireEvent.click(checkbox);
    expect(p.setShowPhase).toHaveBeenCalledWith(true);
  });

  it('parses numeric values from sliders (amplitude)', () => {
    const p = renderControls({ amp: 1.0 });
    // find the label text node for amplitude and then the sibling input
    const label = screen.getByText(/amplitude\s*=/i);
    const container = label.closest('div');
    // the input is inside the sibling div
    const wrapper = container?.nextElementSibling;
    const input = wrapper?.querySelector('input[type="range"]');
    expect(input).toBeTruthy();
    fireEvent.change(input, { target: { value: '1.3' } });
    expect(p.setAmp).toHaveBeenCalledWith(1.3);
  });
});
