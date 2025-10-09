import { defineConfig } from '@playwright/test';

export default defineConfig({
  testDir: 'tests/e2e',
  fullyParallel: false,
  retries: process.env.CI ? 1 : 0,
  timeout: 60_000,
  expect: {
    toHaveScreenshot: {
      // Use a cross-OS snapshot path to avoid -linux/-darwin suffix divergence
      // Resulting path for toHaveScreenshot('free-packet.png'):
      // tests/e2e/app.spec.ts-snapshots/free-packet.png
      pathTemplate: '{testDir}/{testFileDir}/{testFileName}-snapshots/{arg}{ext}',
    },
  },
  use: {
    baseURL: 'http://localhost:5173',
    viewport: { width: 1280, height: 800 },
    deviceScaleFactor: 1,
    colorScheme: 'dark',
    screenshot: 'only-on-failure',
    video: 'retain-on-failure',
  },
  webServer: [
    {
      command: 'npm run dev -- --port=5173 --strictPort',
      url: 'http://localhost:5173',
      reuseExistingServer: !process.env.CI,
      timeout: 120_000,
    },
  ],
});
