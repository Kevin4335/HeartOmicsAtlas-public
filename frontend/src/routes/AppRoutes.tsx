import { Routes, Route, Navigate } from "react-router-dom";

import Home from "../pages/Home";
import AIChat from "../pages/AIChat";
import SpatialTranscriptomics from "../pages/SpatialTranscriptomics";
import Multiomics from "../pages/Multiomics";

import ScRNALayout from "../pages/scrna/ScRNALayout";
import ScrnaAcmVcmSan from "../pages/scrna/ScrnaAcmVcmSan";
import ScrnaSanPco from "../pages/scrna/ScrnaSanPco";
import ScrnaMiniHeart from "../pages/scrna/ScrnaMiniHeart";

import NotFound from "../pages/NotFound";

export default function AppRoutes() {
  return (
    <Routes>
      <Route path="/" element={<Home />} />
      <Route path="/chat" element={<AIChat />} />
      <Route path="/st" element={<SpatialTranscriptomics />} />
      <Route path="/multiomics" element={<Multiomics />} />

      {/* scRNA hub + subpages */}
      <Route path="/scrna" element={<ScRNALayout />}>
        <Route index element={<Navigate to="acm-vcm-san" replace />} />
        <Route path="acm-vcm-san" element={<ScrnaAcmVcmSan />} />
        <Route path="san-pco" element={<ScrnaSanPco />} />
        <Route path="mini-heart" element={<ScrnaMiniHeart />} />
      </Route>

      <Route path="/404" element={<NotFound />} />
      <Route path="*" element={<Navigate to="/404" replace />} />
    </Routes>
  );
}
