# 3D Quantum Wave Engine — Usage Guide

This guide provides step-by-step instructions for installing dependencies,
running the development server, and using the app locally.

See Concepts & Reference: [usage/concepts.md](./usage/concepts.md)

## prerequisites

- **node.js**: version 16 or higher (required for vite and react). download from [nodejs.org](https://nodejs.org/).
- **npm**: comes bundled with node.js. version 7 or higher is recommended.
- **modern web browser**: chrome, firefox, safari, or edge with webgl support enabled (required for three.js rendering).

to verify your node.js and npm versions, run:

```bash
node --version
npm --version
```

## installation

1. **clone or download the project**:
   - if using git: `git clone https://github.com/alzweidi/qm3d.git`
   - or download the project zip and extract it to your desired location.

2. **navigate to the project directory**:

   ```bash
   cd /path/to/qm3d
   ```

3. **install dependencies**:
   run the following command to install all required packages:

   ```bash
   npm install
   ```

   this will install react, reactdom, three.js, vite, tailwind css, and other development dependencies as specified in `package.json`.

## running the development server

1. **start the development server**:

   ```bash
   npm run dev
   ```

2. **access the application**:
   - open your web browser and navigate to `http://localhost:5173` (the default vite development server url).
   - the 3d quantum wave engine should load automatically.

3. **stop the server**:
   - press `Ctrl + C` (or `Cmd + C` on macos) in the terminal to stop the development server.

## Using the app

the application consists of a 3d visualisation canvas and a control panel. here's how to use it:

### basic operation

1. **add a wave packet**:
   - use the "packet" section in the controls to configure:
     - **σ (width)**: controls the spread of the gaussian packet (0.3 to 1.5).
     - **centre (x, y, z)**: position of the packet centre (-l/2 to l/2).
     - **k0 (x, y, z)**: wave vector components (-20 to 20).
     - **amplitude**: overall strength of the packet (0.5 to 2.0).
   - click "add packet" to inject the wave packet into the simulation.

2. **choose a potential**:
   - click one of the preset buttons:
     - **free**: no potential (flat landscape).
     - **plane barrier**: a potential wall at x > 0.
     - **box well**: a cubic potential well in the centre.
     - **spherical well**: a spherical potential well.
     - **harmonic**: harmonic oscillator potential (v ∝ r²).

3. **control the simulation**:
   - **run/pause**: toggle simulation playback.
   - **reset ψ**: clear the wave function to zero.
   - **grid n**: simulation resolution (32 or 64). keep at 32 for smooth performance.
   - **domain l**: physical size of the simulation box (6 to 16).
   - **dt scale**: time step size (affects accuracy vs. speed).
   - **steps/frame**: number of simulation steps per animation frame (1-4).

4. **visualisation controls**:
   - **density scale**: adjusts point sizes (0.4 to 2.0).
   - **phase hue**: toggle colour coding by wave phase (on/off).
   - **orbit controls**: click and drag to rotate the view, scroll to zoom.

### absorbing boundaries (cap)

- **width**: fraction of the domain edge used for absorption (5% to 30%).
- **strength**: how strongly waves are damped at boundaries (1 to 8).

this prevents unwanted reflections from the periodic boundaries imposed by the fft.

### tips for effective use

- start with n=32 for responsive interaction.
- use smaller dt values for more accurate simulations.
- the harmonic and spherical well potentials are good for testing bound states.
- points represent |ψ|² (probability density) by size and arg(ψ) (phase) by colour when phase hue is enabled.

## troubleshooting

### common issues

1. **development server won't start**:
   - **cause**: port 5173 may be in use or node.js version too old.
   - **solution**: try a different port with `npm run dev -- --port 3000`, or update node.js.

2. **dependencies installation fails**:
   - **cause**: network issues or permission problems.
   - **solution**: clear npm cache (`npm cache clean --force`) and retry. use `sudo` if on linux/macos and permission denied.

3. **application doesn't load or shows errors**:
   - **cause**: missing dependencies or browser compatibility.
   - **solution**: ensure all dependencies are installed. check browser console for errors. verify webgl is enabled.

4. **poor performance or lag**:
   - **cause**: high grid resolution or too many steps per frame.
   - **solution**: reduce n to 32, lower steps/frame, or increase dt scale.

5. **visualisation issues**:
   - **cause**: webgl not supported or disabled.
   - **solution**: update browser or enable webgl in browser settings.

6. **linting errors**:
   - run `npm run lint` to check for code style issues.

### additional resources

- **three.js documentation**: [threejs.org/docs](https://threejs.org/docs/)
- **react documentation**: [react.dev](https://react.dev/)
- **vite documentation**: [vitejs.dev](https://vitejs.dev/)

if you encounter issues not covered here, check the browser console for error messages and ensure your environment meets the prerequisites.
