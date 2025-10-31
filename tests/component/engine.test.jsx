// @vitest-environment jsdom
import React from 'react';
import { beforeAll, describe, it, expect, vi } from 'vitest';
import { render, screen, fireEvent, within, waitFor } from '@testing-library/react';

// Mock the visualisation layer to avoid WebGL/Three
vi.mock('../../src/rendering/visualisation.js', () => {
  const updatePointCloud = vi.fn();
  let lastPointCloud = null;
  function initialiseThreeJS() {
    // minimal stubs
    const renderer = {
      setPixelRatio: vi.fn(),
      setSize: vi.fn(),
      getPixelRatio: () => 1,
      getContext: () => ({ getParameter: () => [1, 64], ALIASED_POINT_SIZE_RANGE: 'ALIASED_POINT_SIZE_RANGE' }),
      render: vi.fn(),
      dispose: vi.fn(),
      domElement: (typeof document === 'undefined') ? null : document.createElement('canvas'),
    };
    const camera = { fov: 45, updateProjectionMatrix: vi.fn(), aspect: 1 };
    const scene = { add: vi.fn() };
    const controls = { update: vi.fn(), addEventListener: vi.fn(), removeEventListener: vi.fn(), dispose: vi.fn() };
    return { scene, camera, renderer, controls, maxPointSize: 64 };
  }

  function createQuantumPointCloud(_coord, N, maxPointSize, showPhase) {
    const size = N * N * N;
    const ampAttr = { array: new Float32Array(size), needsUpdate: false };
    const phaseAttr = { array: new Float32Array(size), needsUpdate: false };
    const material = { uniforms: { uShowPhase: { value: showPhase ? 1 : 0 }, uSizeScale: { value: 120 }, uMaxPointSize: { value: maxPointSize } } };
    const geometry = { attributes: { aAmp: ampAttr, aPhase: phaseAttr }, dispose: vi.fn() };
    const points = { material, frustumCulled: false };
    lastPointCloud = { points, geometry, material, ampAttr, phaseAttr };
    return { points, geometry, material, ampAttr, phaseAttr };
  }

  function createBoundingBox() { return { geometry: { dispose: vi.fn() }, material: { dispose: vi.fn() } }; }
  const handleResize = vi.fn();
  const renderScene = vi.fn();
  const disposeThreeJS = vi.fn();

  return {
    initialiseThreeJS,
    createQuantumPointCloud,
    createBoundingBox,
    updatePointCloud,
    handleResize,
    renderScene,
    disposeThreeJS,
    __getLastPointCloud: () => lastPointCloud,
  };
});

// Ensure RAF exists in JSDOM
beforeAll(() => {
  if (globalThis.requestAnimationFrame == null) {
    globalThis.requestAnimationFrame = (cb) => setTimeout(cb, 0);
  }
  if (globalThis.cancelAnimationFrame == null) {
    globalThis.cancelAnimationFrame = (id) => clearTimeout(id);
  }
  Object.defineProperty(document, 'hidden', { value: false, configurable: true });
});

import QuantumWaveEngine from '../../src/components/engine.jsx';
import * as Vis from '../../src/rendering/visualisation.js';
import * as Quantum from '../../src/physics/quantum.js';

describe('QuantumWaveEngine (mocked visualisation)', () => {
  it('toggling phase hue updates shader uniform', async () => {
    render(<QuantumWaveEngine />);

    await waitFor(() => expect(Vis.__getLastPointCloud()).toBeTruthy());
    const cloud = Vis.__getLastPointCloud();
    expect(cloud.material.uniforms.uShowPhase.value).toBe(1);

    // Toggle the specific 'phase hue' checkbox inside Controls card
    const controlsHeading = screen.getByRole('heading', { name: /controls/i });
    const controlsCard = controlsHeading.closest('.card');
    expect(controlsCard).toBeTruthy();
    const wc = within(controlsCard);
    const phaseLabel = wc.getByText(/phase hue/i);
    const row = phaseLabel.closest('div');
    const inputContainer = row?.nextElementSibling;
    const checkbox = inputContainer?.querySelector('input[type="checkbox"]');
    expect(checkbox).toBeTruthy();
    fireEvent.click(checkbox);

    // uShowPhase should switch to 0 (via effect)
    await waitFor(() => expect(cloud.material.uniforms.uShowPhase.value).toBe(0));
  });

  it('renders an FPS overlay in the main view', async () => {
    render(<QuantumWaveEngine />);

    await waitFor(() => expect(Vis.__getLastPointCloud()).toBeTruthy());

    const mainCard = document.querySelector(String.raw`.lg\:col-span-2.card.p-4`);
    expect(mainCard).toBeTruthy();
    // Should show either an FPS readout or a mode label
    const hasFps = !!within(mainCard).queryByText(/\bfps\s+\d+\.\d\b/i);
    const hasMode = !!within(mainCard).queryByText(/\b(paused|idle)\b/i);
    expect(hasFps || hasMode).toBe(true);
  });

  it('overlay switches to paused when toggled and to idle when tab is hidden', async () => {
    render(<QuantumWaveEngine />);

    await waitFor(() => expect(Vis.__getLastPointCloud()).toBeTruthy());

    const mainCard = document.querySelector(String.raw`.lg\:col-span-2.card.p-4`);
    const utils = within(mainCard);

    // Pause
    const toggleBtn = utils.getByText(/pause|run/i);
    expect(toggleBtn).toBeTruthy();
    fireEvent.click(toggleBtn); // running -> paused

    await waitFor(() => expect(utils.getByText(/\bpaused\b/i)).toBeTruthy());

    // Simulate tab hidden
    Object.defineProperty(document, 'hidden', { value: true, configurable: true });
    document.dispatchEvent(new Event('visibilitychange'));

    await waitFor(() => expect(utils.getByText(/\bidle\b/i)).toBeTruthy());

    // Back to visible and run again
    Object.defineProperty(document, 'hidden', { value: false, configurable: true });
    document.dispatchEvent(new Event('visibilitychange'));
    fireEvent.click(utils.getByText(/pause|run/i)); // paused -> running

    await waitFor(() => {
      const maybeFps = utils.queryByText(/\bfps\s+\d+\.\d\b/i);
      expect(maybeFps).toBeTruthy();
    });
  });

  it('default k0x ensures >=8 points-per-wavelength for defaults', async () => {
    const spy = vi.spyOn(Quantum, 'addPacket3D');
    render(<QuantumWaveEngine />);

    await waitFor(() => expect(Vis.__getLastPointCloud()).toBeTruthy());

    const mainCard = document.querySelector(String.raw`.lg\:col-span-2.card.p-4`);
    const utils = within(mainCard);
    fireEvent.click(utils.getByText(/add packet/i));

    await waitFor(() => expect(spy).toHaveBeenCalled());
    const call = spy.mock.calls.at(-1);
    const kxArg = call[10];
    const N = 32, L = 10;
    const dx = L / N;
    const ppw = (2 * Math.PI / Math.max(1e-12, Math.abs(kxArg))) / dx;
    expect(ppw).toBeGreaterThanOrEqual(7.9);
  });

  it('add packet triggers updatePointCloud', async () => {
    Vis.updatePointCloud.mockClear();
    render(<QuantumWaveEngine />);

    await waitFor(() => expect(Vis.__getLastPointCloud()).toBeTruthy());

    // Click the 'add packet' button in the main left panel (not the Controls panel)
    const mainCard = document.querySelector(String.raw`.lg\:col-span-2.card.p-4`);
    expect(mainCard).toBeTruthy();
    const utils = within(mainCard);
    fireEvent.click(utils.getByText(/add packet/i));

    expect(Vis.updatePointCloud).toHaveBeenCalled();
  });

  it('clamps k0 to kMax before addPacket', async () => {
    const spy = vi.spyOn(Quantum, 'addPacket3D');
    render(<QuantumWaveEngine />);

    await waitFor(() => expect(Vis.__getLastPointCloud()).toBeTruthy());

    // Find Controls card
    const controlsHeading = screen.getByRole('heading', { name: /controls/i });
    const controlsCard = controlsHeading.closest('.card');
    const withinControls = within(controlsCard);

    // Locate k0x slider and set to a very large value to trigger clamp
    const k0xLabel = withinControls.getByText(/k0x\s*=/i);
    const k0xRow = k0xLabel.closest('div');
    const k0xInput = k0xRow?.nextElementSibling?.querySelector('input[type="range"]');
    expect(k0xInput).toBeTruthy();
    // Change to an out-of-range value
    fireEvent.change(k0xInput, { target: { value: '9999' } });

    // Click the main panel 'add packet'
    const mainCard = document.querySelector(String.raw`.lg\:col-span-2.card.p-4`);
    const utils = within(mainCard);
    fireEvent.click(utils.getByText(/add packet/i));

    // Assert addPacket3D was called with clamped kx argument
    await waitFor(() => expect(spy).toHaveBeenCalled());
    const call = spy.mock.calls.at(-1);
    // args: [..., sigma,sigma,sigma, kx, ky, kz, ...]
    const kxArg = call[10];
    const N = 32, L = 10; // engine defaults
    const dx = L / N;
    const kMax = Math.PI / dx * 0.9;
    expect(Math.abs(kxArg)).toBeLessThanOrEqual(kMax + 1e-12);
  });

  it('re-clamps k0 when grid spacing increases (L larger)', async () => {
    const spy = vi.spyOn(Quantum, 'addPacket3D');
    render(<QuantumWaveEngine />);

    await waitFor(() => expect(Vis.__getLastPointCloud()).toBeTruthy());

    // Increase L via Controls to reduce kMax
    const controlsHeading = screen.getByRole('heading', { name: /controls/i });
    const controlsCard = controlsHeading.closest('.card');
    const { getByText } = within(controlsCard);

    // Find domain slider row by label containing 'domain l'
    const lLabel = getByText(/domain l/i);
    const lRow = lLabel.closest('div');
    const lInput = lRow?.nextElementSibling?.querySelector('input[type="range"]');
    expect(lInput).toBeTruthy();
    // Set L to 16 (from min 6 to max 16 per component)
    fireEvent.change(lInput, { target: { value: '16' } });

    // Now click add packet, engine will clamp again
    const mainCard = document.querySelector(String.raw`.lg\:col-span-2.card.p-4`);
    const utils = within(mainCard);
    fireEvent.click(utils.getByText(/add packet/i));

    await waitFor(() => expect(spy).toHaveBeenCalled());
    const call = spy.mock.calls.at(-1);
    const kxArg = call[10];
    const N = 32, L = 16;
    const dx = L / N;
    const kMax = Math.PI / dx * 0.9; // ~5.6549
    expect(Math.abs(kxArg)).toBeLessThanOrEqual(kMax + 1e-12);
  });
});
