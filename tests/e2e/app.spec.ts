import { test, expect } from '@playwright/test';
import type { Page, Locator } from '@playwright/test';

async function getMainCard(page: Page): Promise<Locator> {
  // left-hand main card container
  return page.locator(String.raw`div.lg\:col-span-2.card.p-4`).first();
}

async function getCanvas(page: Page): Promise<Locator> {
  return page.locator('canvas').first();
}

// For first-run, create baselines with: npx playwright test --update-snapshots

test.describe('App E2E', () => {
  test('loads app and renders canvas without console errors', async ({ page }) => {
    const errors: string[] = [];
    page.on('console', (msg) => {
      if (msg.type() === 'error') errors.push(msg.text());
    });

    await page.goto('/');
    await expect(await getCanvas(page)).toBeVisible();
    // Allow a small idle to avoid transient dev server warnings
    await page.waitForTimeout(200);
    expect(errors, `Console errors: ${errors.join('\n')}`).toHaveLength(0);
  });

  test('free space: add packet then pause and snapshot', async ({ page }) => {
    await page.goto('/');

    const main = await getMainCard(page);
    await expect(main).toBeVisible();

    // add packet
    await main.getByRole('button', { name: /add packet/i }).click();
    // let it evolve for a short while, then pause for stable capture
    await page.waitForTimeout(400);
    await main.getByRole('button', { name: /pause/i }).click();

    const canvas = await getCanvas(page);
    await expect(canvas).toHaveScreenshot('free-packet.png', {
      maxDiffPixelRatio: 0.02,
    });
  });

  test('plane barrier: add packet then pause and snapshot', async ({ page }) => {
    await page.goto('/');

    const main = await getMainCard(page);
    await expect(main).toBeVisible();

    await main.getByRole('button', { name: /reset ψ/i }).click();
    await main.getByRole('button', { name: /plane barrier/i }).click();
    await main.getByRole('button', { name: /add packet/i }).click();

    await page.waitForTimeout(500);
    await main.getByRole('button', { name: /pause/i }).click();

    const canvas = await getCanvas(page);
    await expect(canvas).toHaveScreenshot('barrier.png', {
      maxDiffPixelRatio: 0.02,
    });
  });

  test('harmonic: add packet then pause and snapshot', async ({ page }) => {
    await page.goto('/');

    const main = await getMainCard(page);
    await expect(main).toBeVisible();

    await main.getByRole('button', { name: /reset ψ/i }).click();
    await main.getByRole('button', { name: /harmonic/i }).click();
    await main.getByRole('button', { name: /add packet/i }).click();

    await page.waitForTimeout(500);
    await main.getByRole('button', { name: /pause/i }).click();

    const canvas = await getCanvas(page);
    await expect(canvas).toHaveScreenshot('harmonic.png', {
      maxDiffPixelRatio: 0.02,
    });
  });
});
