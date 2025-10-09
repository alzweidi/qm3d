# 3D Quantum Wave Engine

An interactive 3D simulator for quantum wavefunctions in the browser.
Built with React, Three.js, and a split-step FFT solver for time
evolution under common potentials.

## Features

- Gaussian wave packet injection with adjustable parameters
- Presets: free space, plane barrier, box well, spherical well, harmonic oscillator
- Real-time visualization with phase coloring and orbit controls
- Absorbing boundaries (CAP) to mitigate periodic FFT wrap-around
- Adjustable grid size, time step, and rendering options

## Quick start

Prerequisites: Node.js 16+ (includes npm 7+).

```bash
npm install
npm run dev
```

Open [http://localhost:5173](http://localhost:5173) in your browser.

Or

Visit the site

<https://qm3d.netlify.app/>

## Usage & docs

- See the Usage Guide: [usage.md](./usage.md)
- Concepts & Reference: [usage/concepts.md](./usage/concepts.md)
- Troubleshooting lives in the guide: [usage.md#troubleshooting](./usage.md#troubleshooting)

## Tech stack

React • Three.js • Vite • Tailwind CSS

## Why

Built as a summer project to explore split-step spectral solvers and WebGL.
There’s also a 2D version in another repo.
