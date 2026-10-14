import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// https://vite.dev/config/
export default defineConfig({
  plugins: [react()],
  server: {
    port: 5173,
    proxy: {
      '/api': {
        target: 'http://localhost:8550',
        changeOrigin: true,
      },
      // Served by the CoolSolve binary too: the documentation and the model library
      '/docs': {
        target: 'http://localhost:8550',
        changeOrigin: true,
      },
      '/library': {
        target: 'http://localhost:8550',
        changeOrigin: true,
      },
    },
  },
})
