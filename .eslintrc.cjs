/* ESLint configuration for React + Vite + Vitest */
module.exports = {
  root: true,
  env: {
    browser: true,
    node: true,
    es2021: true,
  },
  settings: {
    react: { version: 'detect' },
  },
  parserOptions: {
    ecmaVersion: 'latest',
    sourceType: 'module',
    ecmaFeatures: { jsx: true },
  },
  plugins: ['react', 'react-hooks', 'react-refresh'],
  extends: [
    'eslint:recommended',
    'plugin:react/recommended',
    'plugin:react-hooks/recommended',
    'plugin:react-refresh/recommended',
  ],
  rules: {
    // using the new JSX transform, React in scope is not required
    'react/react-in-jsx-scope': 'off',
    // do not use prop-types in this project
    'react/prop-types': 'off',
  },
  overrides: [
    {
      files: ['tests/**/*.test.{js,jsx}', 'tests/setup.js'],
      env: { browser: true, node: true, jest: true },
      globals: { vi: 'readonly' },
    },
  ],
};
