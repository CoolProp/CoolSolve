import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// https://vite.dev/config/
export default defineConfig({
  plugins: [react()],
  // 127.0.0.1, not localhost: the CoolSolve server listens on IPv4 only, and Node 18 tries ::1 first
  server: {
    port: 5173,
    proxy: {
      '/api': {
        target: 'http://127.0.0.1:8550',
        changeOrigin: true,
      },
      // Served by the CoolSolve binary too: the documentation and the model library
      '/docs': {
        target: 'http://127.0.0.1:8550',
        changeOrigin: true,
      },
      '/library': {
        target: 'http://127.0.0.1:8550',
        changeOrigin: true,
      },
    },
  },
})
