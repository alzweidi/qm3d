# Contributing

Thanks for your interest in contributing! This document explains how to work on the project and submit changes.

## Code of Conduct

By participating, you agree to abide by our Code of Conduct (see `CODE_OF_CONDUCT.md`).

## Prerequisites

- Node.js 20.19+ recommended
- npm (or pnpm/yarn). The repo uses npm scripts; pnpm is preferred if you use it locally.

## Getting started

```bash
# install deps
npm install

# start dev server
npm run dev

# run linters
npm run lint

# run unit + component tests
npm run test

# run just component tests (jsdom)
npm run test:component

# run E2E tests (Playwright)
npm run test:e2e
```

## Project layout (high‑level)

- `src/` React UI and Three.js rendering/physics glue
- `src/physics/` FFT + quantum simulation code
- `tests/` unit, component, and E2E tests

## Style & quality

- Follow existing patterns in the codebase.
- ESLint must pass: `npm run lint`.
- Tests required for new/changed behavior. Aim for high coverage; add unit and/or component tests. For visual changes, prefer an E2E screenshot test.
- Prefer small, focused PRs.

## Commit messages

- Clear and descriptive. Conventional Commits are welcome (e.g. `feat:`, `fix:`, `test:`) but not required.

## Pull requests

- Link the related issue (e.g. `Closes #123`).
- Include:
  - Summary of change and rationale.
  - Screenshots/GIFs for UI changes.
  - Notes on performance/UX impact when relevant.
- Checks to pass:
  - `npm run lint`
  - `npm run test` (and `npm run test:e2e` when applicable)

## Issue triage

- Use labels: `bug`, `feature`, `documentation`, `help wanted`, etc.
- Provide a minimal reproduction when filing bugs.

## Security

Please do not report vulnerabilities via public issues. See `SECURITY.md` for private reporting instructions.
