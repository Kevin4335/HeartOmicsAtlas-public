import type { PropsWithChildren } from "react";
import { Link as RouterLink, useLocation } from "react-router-dom";
import logo from "../assets/heart_logo_1.svg";
import EmailOutlined from "@mui/icons-material/EmailOutlined";
import PhoneOutlined from "@mui/icons-material/PhoneOutlined";
import LocationOnOutlined from "@mui/icons-material/LocationOnOutlined";

import {
  AppBar,
  Box,
  Link,
  Stack,
  Toolbar,
  Typography,
} from "@mui/material";

function NavLink(props: { to: string; label: string; matchPaths?: string[] }) {
  const location = useLocation();
  // Check if current path matches this nav link (exact match or starts with for nested routes)
  const isActive = props.matchPaths 
    ? props.matchPaths.some(path => location.pathname === path || location.pathname.startsWith(path + "/"))
    : location.pathname === props.to || (props.to !== "/" && location.pathname.startsWith(props.to + "/"));
  
  // Special case for home - only exact match
  const isHomeActive = props.to === "/" && location.pathname === "/";
  const active = props.to === "/" ? isHomeActive : isActive;

  return (
    <Link
      component={RouterLink}
      to={props.to}
      underline="none"
      sx={{
        color: active ? "#C30F1A" : "#000000",
        fontWeight: active ? 700 : 500,
        fontSize: 15,
        px: 1.5,
        py: 0.5,
        borderRadius: 0,
        position: "relative",
        textDecoration: "none",
        display: "inline-flex",
        flexDirection: "column",
        alignItems: "center",
        // Invisible bold text to reserve space and prevent layout shift
        "&::before": {
          content: `"${props.label}"`,
          fontWeight: 700,
          height: 0,
          visibility: "hidden",
          overflow: "hidden",
          userSelect: "none",
          pointerEvents: "none",
        },
        "&::after": active ? {
          content: '""',
          position: "absolute",
          bottom: 0,
          left: "50%",
          transform: "translateX(-50%)",
          width: "80%",
          height: "2px",
          backgroundColor: "#C30F1A",
        } : {},
      }}
    >
      {props.label}
    </Link>
  );
}

export default function AppShell({ children }: PropsWithChildren) {
  const location = useLocation();
  
  // Pages that should hide the full footer (they have their own simplified footer)
  const hideFooter = location.pathname.startsWith("/explore") || location.pathname.startsWith("/chat");

  return (
    <Box sx={{ minHeight: "100vh", display: "flex", flexDirection: "column", backgroundColor: "#ffffff" }}>
      <AppBar 
        position="sticky" 
        elevation={0}
        sx={{
          backgroundColor: "#ffffff",
          boxShadow: "0 1px 3px rgba(0, 0, 0, 0.05)",
        }}
      >
        <Box sx={{ px: "14%" }}>
          <Toolbar disableGutters sx={{ height: "8.82vh" }}>
            {/* Left - 50% width */}
            <Stack direction="row" alignItems="center" spacing={1.5} sx={{ width: "50%" }}>
              <Box
                component="img"
                src={logo}
                alt="HeartOmicsAtlas logo"
                sx={{
                  width: 52,
                  height: 52,
                  display: "block",
                }}
              />
              <Typography variant="h5" component="span" sx={{ fontWeight: 600, letterSpacing: 0.3 }}>
                <Box component="span" sx={{ color: "#000000" }}>Heart</Box>
                <Box component="span" sx={{ color: "#BE1B23" }}>Omics</Box>
                <Box component="span" sx={{ color: "#000000" }}>Atlas</Box>
              </Typography>
            </Stack>

            {/* Right - 50% width */}
            <Stack direction="row" spacing={2} alignItems="center" justifyContent="flex-end" sx={{ width: "50%" }}>
              <NavLink to="/" label="Home" />
              <NavLink to="/explore" label="Explore Atlas" />
              <NavLink to="/chat" label="AI Assistant" />

              <Link
                href="https://www.chen-stemcell-lab.com/"
                target="_blank"
                rel="noopener noreferrer"
                underline="none"
                sx={{
                  color: "#000000",
                  fontWeight: 500,
                  fontSize: 15,
                  px: 1.5,
                  py: 0.5,
                  textDecoration: "none",
                }}
              >
                Chen Lab
              </Link>
            </Stack>
          </Toolbar>
        </Box>
      </AppBar>

      {/* Main */}
      <Box component="main" sx={{ flex: 1, backgroundColor: "#ffffff" }}>
        {children}
      </Box>

      {/* Footer - hidden on pages with their own footer */}
      {!hideFooter && (
        <Box
          component="footer"
          sx={{
            width: "100%",
            backgroundColor: "#980C10",
            color: "white",
            py: { xs: 4, md: 5 },
          }}
        >
          {/* Working area */}
          <Box sx={{ px: "14%" }}>
            {/* Top section: Contact */}
            <Box
              sx={{
                display: "grid",
                gridTemplateColumns: { xs: "1fr", md: "220px 1fr" },
                alignItems: "start",
                columnGap: { xs: 2, md: 8 },
                rowGap: { xs: 2, md: 0 },
              }}
            >
              <Typography sx={{ fontWeight: 800, fontSize: 16, lineHeight: 1.2 }}>
                Contact
              </Typography>

              {/* Push the contact details toward center on desktop */}
              <Box
                sx={{
                  justifySelf: { xs: "start", md: "start" },
                  pl: { xs: 0, md: 8 }, // adjust this to move the block more right/left
                }}
              >
                <Stack spacing={1.1}>
                  <Stack direction="row" spacing={1.1} alignItems="center">
                    <EmailOutlined sx={{ fontSize: 18, color: "white" }} />
                    <Typography sx={{ fontSize: 13, lineHeight: 1.3, color: "white" }}>
                      shc2034@med.cornell.edu
                    </Typography>
                  </Stack>

                  <Stack direction="row" spacing={1.1} alignItems="center">
                    <PhoneOutlined sx={{ fontSize: 18, color: "white" }} />
                    <Typography sx={{ fontSize: 13, lineHeight: 1.3, color: "white" }}>
                      +1 (212) 746-5431
                    </Typography>
                  </Stack>

                  <Stack direction="row" spacing={1.1} alignItems="center">
                    <LocationOnOutlined sx={{ fontSize: 18, color: "white" }} />
                    <Typography sx={{ fontSize: 13, lineHeight: 1.3, color: "white" }}>
                      A827B, 1300 York Ave, New York, NY, 10065
                    </Typography>
                  </Stack>
                </Stack>
              </Box>
            </Box>

            {/* Bottom section: Copyright */}
            <Box
              sx={{
                mt: { xs: 4, md: 5 },
                pt: 2,
                textAlign: "center",
              }}
            >
              <Typography sx={{ fontSize: 12, color: "rgba(255,255,255,0.75)" }}>
                Copyright © 2026 - 2030 Chen Lab @ Weil Cornell Medicine. All rights reserved.
              </Typography>
            </Box>
          </Box>
        </Box>
      )}

    </Box>
  );
}
