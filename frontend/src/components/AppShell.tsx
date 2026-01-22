import type { PropsWithChildren } from "react";
import { Link as RouterLink } from "react-router-dom";

import {
  AppBar,
  Box,
  Container,
  Toolbar,
  Typography,
  Link,
  Stack,
} from "@mui/material";

function NavLink(props: { to: string; label: string }) {
  return (
    <Link
      component={RouterLink}
      to={props.to}
      underline="none"
      sx={{
        color: "inherit",
        fontWeight: 500,
        fontSize: 15,
        px: 1.5,
        py: 0.5,
        borderRadius: 1,
        transition: "background-color 0.2s",
        "&:hover": {
          backgroundColor: "rgba(255,255,255,0.12)",
        },
      }}
    >
      {props.label}
    </Link>
  );
}

export default function AppShell({ children }: PropsWithChildren) {
  return (
    <Box sx={{ minHeight: "100vh", display: "flex", flexDirection: "column" }}>
      <AppBar position="sticky" elevation={0}>
        <Container maxWidth="lg">
          <Toolbar disableGutters sx={{ height: 64 }}>
            
            {/* Left: Logo + Title */}
            <Stack direction="row" alignItems="center" spacing={1.5} sx={{ flexGrow: 1 }}>
              
              {/* Simple logo circle. Replace later with real SVG if desired */}
              <Box
                sx={{
                  width: 28,
                  height: 28,
                  borderRadius: "50%",
                  backgroundColor: "white",
                  opacity: 0.9,
                }}
              />

              <Typography
                variant="h6"
                sx={{
                  fontWeight: 700,
                  letterSpacing: 0.3,
                }}
              >
                HeartOmicsAtlas
              </Typography>
            </Stack>

            {/* Right: Navigation */}
            <Stack direction="row" spacing={1}>
              <NavLink to="/" label="Home" />
              <NavLink to="/chat" label="Chat with AI" />
              <NavLink to="/st" label="Spatial Transcriptomics" />
              <NavLink to="/multiomics" label="Multiomics" />
              <NavLink to="/scrna" label="scRNA" />
            </Stack>

          </Toolbar>
        </Container>
      </AppBar>

      {/* Main content */}
      <Box component="main" sx={{ flex: 1, py: 4 }}>
        <Container maxWidth="lg">{children}</Container>
      </Box>

      {/* Footer */}
      <Box component="footer" sx={{ py: 3 }}>
        <Container maxWidth="lg">
          <Typography variant="body2" color="text.secondary">
            Copyright © Chen lab at Weill Cornell Medicine {new Date().getFullYear()} All rights reserved.
          </Typography>
        </Container>
      </Box>
    </Box>
  );
}
