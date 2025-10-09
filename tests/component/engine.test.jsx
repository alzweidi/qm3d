// @vitest-environment jsdom
import React from 'react';
import { beforeAll, describe, it, expect, vi } from 'vitest';
import { render, screen, fireEvent, within, waitFor } from '@testing-library/react';

// Mock the visualisation layer to avoid WebGL/Three
vi.mock('../../src/rendering/visualisation.js', () => {
  const updatePointCloud = vi.fn();
  let lastPointCloud = null;
  function initialiseThreeJS(mountElement, L) {
    // minimal stubs
    const renderer = {
      setPixelRatio: vi.fn(),
      setSize: vi.fn(),
      getPixelRatio: () => 1,
      getContext: () => ({ getParameter: () => [1, 64], ALIASED_POINT_SIZE_RANGE: 'ALIASED_POINT_SIZE_RANGE' }),
      render: vi.fn(),
      dispose: vi.fn(),
      domElement: (typeof document !== 'undefined') ? document.createElement('canvas') : null,
    };
    const camera = { fov: 45, updateProjectionMatrix: vi.fn(), aspect: 1 };
    const scene = { add: vi.fn() };
    const controls = { update: vi.fn(), addEventListener: vi.fn(), removeEventListener: vi.fn(), dispose: vi.fn() };
    return { scene, camera, renderer, controls, maxPointSize: 64 };
  }

  function createQuantumPointCloud(coord, N, maxPointSize, showPhase) {
    const size = N * N * N;
    const ampAttr = { array: new Float32Array(size), needsUpdate: false };
    const phaseAttr = { array: new Float32Array(size), needsUpdate: false };
    const material = { uniforms: { uShowPhase: { value: showPhase ? 1.0 : 0.0 }, uSizeScale: { value: 120.0 }, uMaxPointSize: { value: maxPointSize } } };
    const geometry = { attributes: { aAmp: ampAttr, aPhase: phaseAttr }, dispose: vi.fn() };
    const points = { material, frustumCulled: false };
    lastPointCloud = { points, geometry, material, ampAttr, phaseAttr };
    return { points, geometry, material, ampAttr, phaseAttr };
  }

  function createBoundingBox(L) { return { geometry: { dispose: vi.fn() }, material: { dispose: vi.fn() } }; }
  function handleResize() {}
  function renderScene() {}
  function disposeThreeJS() {}

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
  if (!globalThis.requestAnimationFrame) {
    globalThis.requestAnimationFrame = (cb) => setTimeout(cb, 0);
  }
  if (!globalThis.cancelAnimationFrame) {
    globalThis.cancelAnimationFrame = (id) => clearTimeout(id);
  }
  Object.defineProperty(document, 'hidden', { value: false, configurable: true });
});

import QuantumWaveEngine from '../../src/components/engine.jsx';
import * as Vis from '../../src/rendering/visualisation.js';

describe('QuantumWaveEngine (mocked visualisation)', () => {
  it('toggling phase hue updates shader uniform', async () => {
    render(<QuantumWaveEngine />);

    await waitFor(() => expect(Vis.__getLastPointCloud()).toBeTruthy());
    const cloud = Vis.__getLastPointCloud();
    expect(cloud.material.uniforms.uShowPhase.value).toBe(1);

    // Toggle checkbox via Controls (scope to Controls card and match the label with '=')
    const controlsHeading = screen.getByRole('heading', { name: /controls/i });
    const controlsCard = controlsHeading.closest('.card');
    expect(controlsCard).toBeTruthy();
    const { getByText } = within(controlsCard);
    const phaseLabelNode = getByText(/^\s*phase hue\s*=/i);
    const row = phaseLabelNode.closest('div');
    const checkbox = row?.nextElementSibling?.querySelector('input[type="checkbox"]');
    expect(checkbox).toBeTruthy();
    fireEvent.click(checkbox);

    // uShowPhase should switch to 0 (via effect)
    await waitFor(() => expect(cloud.material.uniforms.uShowPhase.value).toBe(0));
  });

  it('add packet triggers updatePointCloud', async () => {
    Vis.updatePointCloud.mockClear();
    render(<QuantumWaveEngine />);

    await waitFor(() => expect(Vis.__getLastPointCloud()).toBeTruthy());

    // Click the 'add packet' button in the main left panel (not the Controls panel)
    const mainCard = document.querySelector('.lg\\:col-span-2.card.p-4');
    expect(mainCard).toBeTruthy();
    const utils = within(mainCard);
    fireEvent.click(utils.getByText(/add packet/i));

    expect(Vis.updatePointCloud).toHaveBeenCalled();
  });
});
