import { defineConfig, type ProxyOptions } from 'vite'
import react from '@vitejs/plugin-react'

// changeOrigin rewrites Host to localhost:8550. The CoolSolve server refuses a POST, PUT or
// DELETE whose Origin (here http://localhost:5173) names neither Host nor X-Forwarded-Host
// (cross-site guard in src/server.cpp), so the host the browser used is passed on.
const forwardHost: ProxyOptions['configure'] = (proxy) => {
  proxy.on('proxyReq', (proxyReq, req) => {
    if (req.headers.host) proxyReq.setHeader('X-Forwarded-Host', req.headers.host)
  })
}

// https://vite.dev/config/
export default defineConfig({
  plugins: [react()],
  server: {
    port: 5173,
    proxy: {
      '/api': {
        target: 'http://localhost:8550',
        changeOrigin: true,
        configure: forwardHost,
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
