# 3D Quantum Wave Engine Usage Guide

this guide provides step-by-step instructions for installing dependencies, running the development server, and using the 3D Quantum Wave Engine locally.

## Prerequisites

- **Node.js**: version 16 or higher (required for Vite and React). Download from [nodejs.org](https://nodejs.org/).
- **npm**: comes bundled with Node.js. version 7 or higher is recommended.
- **Modern Web Browser**: Chrome, Firefox, Safari, or Edge with WebGL support enabled (required for Three.js rendering).

to verify your Node.js and npm versions, run:

```bash
node --version
npm --version
```

## Installation

1. **Clone or Download the Project**:
   - if using Git: `git clone https://github.com/alzweidi/qm3d.git`
   - or download the project ZIP and extract it to your desired location.

2. **Navigate to the Project Directory**:

   ```bash
   cd /path/to/qm3d
   ```

3. **Install Dependencies**:
   run the following command to install all required packages:

   ```bash
   npm install
   ```

   this will install React, ReactDOM, Three.js, Vite, Tailwind CSS, and other development dependencies as specified in `package.json`.

## Running the Development Server

1. **Start the Development Server**:

   ```bash
   npm run dev
   ```

2. **Access the Application**:
   - open your web browser and navigate to `http://localhost:5173` (the default Vite development server URL).
   - the 3D Quantum Wave Engine should load automatically.

3. **Stop the Server**:
   - press `Ctrl + C` (or `Cmd + C` on macOS) in the terminal to stop the development server.

## Using the 3D Quantum Wave Engine

the application consists of a 3D visualization canvas and a control panel. here's how to use it:

### Basic Operation

1. **Add a Wave Packet**:
   - use the "Packet" section in the controls to configure:
     - **σ (Width)**: controls the spread of the Gaussian packet (0.3 to 1.5).
     - **Center (x, y, z)**: position of the packet center (-L/2 to L/2).
     - **k0 (x, y, z)**: wave vector components (-20 to 20).
     - **Amplitude**: overall strength of the packet (0.5 to 2.0).
   - click "Add Packet" to inject the wave packet into the simulation.

2. **Choose a Potential**:
   - click one of the preset buttons:
     - **Free**: no potential (flat landscape).
     - **Plane Barrier**: a potential wall at x > 0.
     - **Box Well**: a cubic potential well in the center.
     - **Spherical Well**: a spherical potential well.
     - **Harmonic**: harmonic oscillator potential (V ∝ r²).

3. **Control the Simulation**:
   - **Run/Pause**: toggle simulation playback.
   - **Reset ψ**: clear the wave function to zero.
   - **Grid N**: simulation resolution (32 or 64). keep at 32 for smooth performance.
   - **Domain L**: physical size of the simulation box (6 to 16).
   - **dt Scale**: time step size (affects accuracy vs. speed).
   - **Steps/Frame**: number of simulation steps per animation frame (1-4).

4. **Visualization Controls**:
   - **Density Scale**: adjusts point sizes (0.4 to 2.0).
   - **Phase Hue**: toggle color coding by wave phase (on/off).
   - **Orbit Controls**: click and drag to rotate the view, scroll to zoom.

### Absorbing Boundaries (CAP)

- **Width**: fraction of the domain edge used for absorption (5% to 30%).
- **Strength**: how strongly waves are damped at boundaries (1 to 8).

this prevents unwanted reflections from the periodic boundaries imposed by the FFT.

### Tips for Effective Use

- start with N=32 for responsive interaction.
- use smaller dt values for more accurate simulations.
- the harmonic and spherical well potentials are good for testing bound states.
- points represent |ψ|² (probability density) by size and arg(ψ) (phase) by color when phase hue is enabled.

## Troubleshooting

### Common Issues

1. **Development Server Won't Start**:
   - **Cause**: port 5173 may be in use or Node.js version too old.
   - **Solution**: try a different port with `npm run dev -- --port 3000`, or update Node.js.

2. **Dependencies Installation Fails**:
   - **Cause**: network issues or permission problems.
   - **Solution**: clear npm cache (`npm cache clean --force`) and retry. use `sudo` if on Linux/macOS and permission denied.

3. **Application Doesn't Load or Shows Errors**:
   - **Cause**: missing dependencies or browser compatibility.
   - **Solution**: ensure all dependencies are installed. check browser console for errors. verify WebGL is enabled.

4. **Poor Performance or Lag**:
   - **Cause**: high grid resolution or too many steps per frame.
   - **Solution**: reduce N to 32, lower steps/frame, or increase dt scale.

5. **Visualization Issues**:
   - **Cause**: WebGL not supported or disabled.
   - **Solution**: update browser or enable WebGL in browser settings.

6. **Linting Errors**:
   - run `npm run lint` to check for code style issues.

### Additional Resources

- **Three.js Documentation**: [threejs.org/docs](https://threejs.org/docs/)
- **React Documentation**: [react.dev](https://react.dev/)
- **Vite Documentation**: [vitejs.dev](https://vitejs.dev/)

if you encounter issues not covered here, check the browser console for error messages and ensure your environment meets the prerequisites.
