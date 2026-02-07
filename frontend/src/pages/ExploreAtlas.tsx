import { useState } from "react";
import { useNavigate } from "react-router-dom";
import {
  Box,
  Button,
  InputBase,
  Stack,
  Typography,
} from "@mui/material";
import SearchIcon from "@mui/icons-material/Search";
import mainBackground from "../assets/main_background.svg";

// Example genes for the buttons
const exampleGenes = ["TNNT2", "MYH7", "SCN5A"];

export default function ExploreAtlas() {
  const navigate = useNavigate();
  const [searchValue, setSearchValue] = useState("");

  const handleSearch = (gene?: string) => {
    const query = gene || searchValue.trim();
    if (query) {
      // Navigate to results page with gene
      navigate(`/explore/results?gene=${encodeURIComponent(query)}`);
    } else {
      // Navigate to results page without gene (shows defaults)
      navigate("/explore/results");
    }
  };

  return (
    <Box
      sx={{
        minHeight: "calc(100vh - 8.82vh)", // Full height minus navbar
        display: "flex",
        flexDirection: "column",
        backgroundImage: `url(${mainBackground})`,
        backgroundSize: "130% auto", // Slightly zoomed out
        backgroundPosition: "30% 30%", // Offset towards top-left
        backgroundRepeat: "no-repeat",
      }}
    >
      {/* Main content - centered, occupies middle 42% of screen */}
      <Box
        sx={{
          flex: 1,
          display: "flex",
          flexDirection: "column",
          justifyContent: "center",
          alignItems: "center",
          py: 6,
        }}
      >
        {/* Content container - exactly 42% of screen width */}
        <Box
          sx={{
            width: "42%",
            display: "flex",
            flexDirection: "column",
            alignItems: "center",
          }}
        >
          {/* Title - semi bold, ~27% of screen width */}
          <Typography
            sx={{
              fontWeight: 600,
              fontSize: "2.4vw", // Scales so title spans ~27% of viewport width
              mb: "40px",
              textAlign: "center",
              whiteSpace: "nowrap",
            }}
          >
            <Box component="span" sx={{ color: "#000000" }}>Explore Heart</Box>
            <Box component="span" sx={{ color: "#BE1B23" }}>Omics</Box>
            <Box component="span" sx={{ color: "#000000" }}>Atlas</Box>
          </Typography>

          {/* Search bar container */}
          <Box
            sx={{
              width: "100%",
              mb: 3,
            }}
          >
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                backgroundColor: "#ffffff",
                borderRadius: "8px",
                boxShadow: "0 2px 12px rgba(0, 0, 0, 0.1)",
                overflow: "hidden",
              }}
            >
              <InputBase
                placeholder="Enter a gene name"
                value={searchValue}
                onChange={(e) => setSearchValue(e.target.value)}
                onKeyDown={(e) => {
                  if (e.key === "Enter") handleSearch();
                }}
                sx={{
                  flex: 1,
                  px: 2,
                  py: 1,
                  fontSize: "1rem",
                }}
              />
              <Button
                variant="contained"
                onClick={() => handleSearch()}
                disableElevation
                sx={{
                  borderRadius: "50px",
                  px: 2,
                  py: 0.4,
                  mr: 2,
                  minWidth: 80,
                  fontSize: "0.9rem",
                  backgroundColor: "#C30F1A",
                  boxShadow: "none",
                  "&:hover": {
                    backgroundColor: "#a00d16",
                    boxShadow: "none",
                  },
                }}
              >
                Search
              </Button>
            </Box>
          </Box>

          {/* "For example:" section */}
          <Box
            sx={{
              width: "100%",
            }}
          >
          <Typography
            sx={{
              color: "#000000",
              fontWeight: 600,
              fontSize: "0.9rem",
              mb: 0.8,
            }}
          >
            For example:
          </Typography>

          {/* Example gene buttons */}
          <Stack direction="row" spacing={1.5}>
            {exampleGenes.map((gene) => (
              <Button
                key={gene}
                onClick={() => handleSearch(gene)}
                sx={{
                  flex: 1,
                  backgroundColor: "#ffffff",
                  color: "#2C2C2B",
                  borderRadius: "6px",
                  boxShadow: "0 2px 8px rgba(0, 0, 0, 0.08)",
                  py: 0.8,
                  px: 1.2,
                  textTransform: "none",
                  fontWeight: 400,
                  fontSize: "0.9rem",
                  justifyContent: "space-between",
                  "&:hover": {
                    backgroundColor: "#fafafa",
                    boxShadow: "0 3px 12px rgba(0, 0, 0, 0.12)",
                  },
                }}
              >
                {gene}
                <SearchIcon sx={{ color: "#C30F1A", fontSize: 18 }} />
              </Button>
            ))}
          </Stack>
          </Box>
        </Box>
      </Box>

      {/* Simplified footer */}
      <Box
        component="footer"
        sx={{
          py: 2,
          textAlign: "center",
        }}
      >
        <Typography
          sx={{
            fontSize: "0.85rem",
            color: "#888888",
          }}
        >
          Copyright © Chen lab at Weill Cornell Medicine {new Date().getFullYear()} All rights reserved.
        </Typography>
      </Box>
    </Box>
  );
}
