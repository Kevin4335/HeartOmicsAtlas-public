import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import path from "node:path";

export default defineConfig({
  plugins: [react()],
  // Load .env from project root (same file as backend uses)
  envDir: path.resolve(__dirname, ".."),
  // Proxy /chat to backend so same-origin fetch works in dev (no CORS). Backend should be on 8000.
  server: {
    proxy: {
      "/chat": "http://localhost:8000",
    },
  },
  resolve: {
    alias: {
      "@": path.resolve(__dirname, "./src"),
    },
  },
});
