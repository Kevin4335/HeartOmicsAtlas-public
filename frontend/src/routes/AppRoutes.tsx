import { Routes, Route, Navigate } from "react-router-dom";

import Home from "../pages/Home";
import AIChat from "../pages/AIChat";
import SpatialTranscriptomics from "../pages/SpatialTranscriptomics";
import Multiomics from "../pages/Multiomics";
import ScRNA from "../pages/ScRNA";
import NotFound from "../pages/NotFound";

export default function AppRoutes() {
  return (
    <Routes>
      <Route path="/" element={<Home />} />
      <Route path="/chat" element={<AIChat />} />
      <Route path="/st" element={<SpatialTranscriptomics />} />
      <Route path="/multiomics" element={<Multiomics />} />
      <Route path="/scrna" element={<ScRNA />} />

      <Route path="/404" element={<NotFound />} />
      <Route path="*" element={<Navigate to="/404" replace />} />
    </Routes>
  );
}
