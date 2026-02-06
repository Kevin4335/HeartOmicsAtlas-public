import type { PropsWithChildren } from "react";
import { Box } from "@mui/material";

/**
 * Standard page wrapper that provides consistent margins (7.7% horizontal padding)
 * and vertical padding for all pages except Home which has its own layout.
 */
export default function PageWrapper({ children }: PropsWithChildren) {
  return (
    <Box sx={{ px: "7.7%", py: 4 }}>
      {children}
    </Box>
  );
}
