import { Routes, Route, Navigate } from "react-router-dom";

import Home from "../pages/Home";
import AIChat from "../pages/AIChat";
import AILanding from "../pages/AILanding";
import ExploreAtlas from "../pages/ExploreAtlas";
import ExploreResults from "../pages/ExploreResults";
import NotFound from "../pages/NotFound";

export default function AppRoutes() {
  return (
    <Routes>
      <Route path="/" element={<Home />} />
      <Route path="/explore" element={<ExploreAtlas />} />
      <Route path="/explore/results" element={<ExploreResults />} />
      <Route path="/chat" element={<AILanding />} />
      <Route path="/chat/conversation" element={<AIChat />} />

      <Route path="/404" element={<NotFound />} />
      <Route path="*" element={<Navigate to="/404" replace />} />
    </Routes>
  );
}
