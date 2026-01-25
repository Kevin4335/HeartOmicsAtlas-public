import type { PropsWithChildren } from "react";
import { Link as RouterLink } from "react-router-dom";
import { useState } from "react";
import logo from "../assets/heart_logo_1.png";

import {
  AppBar,
  Box,
  Container,
  Link,
  Stack,
  Toolbar,
  Typography,
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
        borderRadius: 0,
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
  const [scrnaOpen, setScrnaOpen] = useState(false);

  return (
    <Box sx={{ minHeight: "100vh", display: "flex", flexDirection: "column" }}>
      <AppBar position="sticky" elevation={0}>
        <Container maxWidth="lg">
          <Toolbar disableGutters sx={{ height: 64 }}>
            {/* Left */}
            <Stack direction="row" alignItems="center" spacing={1.5} sx={{ flexGrow: 1 }}>
              <Box
                component="img"
                src={logo}
                alt="HeartOmicsAtlas logo"
                sx={{
                  width: 64,
                  height: 64,
                  display: "block",
                }}
              />
              <Typography variant="h6" sx={{ fontWeight: 700, letterSpacing: 0.3 }}>
                HeartOmicsAtlas
              </Typography>
            </Stack>

            {/* Right */}
            <Stack direction="row" spacing={1} alignItems="center">
              <NavLink to="/" label="Home" />
              <NavLink to="/chat" label="Chat with AI" />
              <NavLink to="/st" label="Spatial Transcriptomics" />
              <NavLink to="/multiomics" label="Multiomics" />

              {/* scRNA hover dropdown (pure CSS-style hover zone, no portal) */}
              <Box
                onMouseEnter={() => setScrnaOpen(true)}
                onMouseLeave={() => setScrnaOpen(false)}
                sx={{
                  position: "relative",
                  display: "inline-flex",
                  alignItems: "center",
                }}
              >
                <NavLink to="/scrna" label="scRNA" />

                {/* Hover bridge + dropdown */}
                {scrnaOpen ? (
                  <Box
                    sx={{
                      position: "absolute",
                      top: "100%",
                      left: 0,
                      pt: 0.5, // creates a small bridge area under the link
                      zIndex: 2000,
                    }}
                  >
                    <Box
                      sx={{
                        width: 220,
                        borderRadius: 0, // perfect rectangle
                        backgroundColor: "#8B0000",
                        color: "#ffffff",
                        border: "1px solid rgba(255,255,255,0.15)",
                        overflow: "hidden",
                      }}
                    >
                      <Link
                        component={RouterLink}
                        to="/scrna/acm-vcm-san"
                        underline="none"
                        sx={{
                          display: "block",
                          color: "inherit",
                          px: 2,
                          py: 1,
                          fontSize: 14,
                          "&:hover": { backgroundColor: "rgba(255,255,255,0.10)" },
                        }}
                      >
                        ACM_VCM_SAN
                      </Link>

                      <Link
                        component={RouterLink}
                        to="/scrna/san-pco"
                        underline="none"
                        sx={{
                          display: "block",
                          color: "inherit",
                          px: 2,
                          py: 1,
                          fontSize: 14,
                          "&:hover": { backgroundColor: "rgba(255,255,255,0.10)" },
                        }}
                      >
                        SAN-PCO
                      </Link>

                      <Link
                        component={RouterLink}
                        to="/scrna/mini-heart"
                        underline="none"
                        sx={{
                          display: "block",
                          color: "inherit",
                          px: 2,
                          py: 1,
                          fontSize: 14,
                          "&:hover": { backgroundColor: "rgba(255,255,255,0.10)" },
                        }}
                      >
                        Mini-heart
                      </Link>
                    </Box>
                  </Box>
                ) : null}
              </Box>
            </Stack>
          </Toolbar>
        </Container>
      </AppBar>

      {/* Main */}
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
