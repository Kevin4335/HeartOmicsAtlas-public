import { useNavigate } from "react-router-dom";
import {
  Box,
  Button,
  Stack,
  Typography,
} from "@mui/material";
import mainBackground from "../assets/main_background.svg";

// Placeholder for the 3 feature items - replace with actual images later
const featureItems = [
  { image: "/placeholder1.png", text: "Spatial Transcriptomics" },
  { image: "/placeholder2.png", text: "Single-cell RNA Sequencing" },
  { image: "/placeholder3.png", text: "Multi-omics Integration" },
];

// Placeholder for the 3 "How Researchers Use" cards
const researcherSteps = [
  {
    image: "/step1.png",
    title: "Explore Data",
    subtitle: "Browse spatial and single-cell datasets",
  },
  {
    image: "/step2.png",
    title: "Analyze Genes",
    subtitle: "Visualize gene expression patterns",
  },
  {
    image: "/step3.png",
    title: "Get Insights",
    subtitle: "AI-powered analysis and interpretation",
  },
];

export default function Home() {
  const navigate = useNavigate();

  return (
    <Box>
      {/* Hero Section - Full width background */}
      <Box
        sx={{
          width: "100%",
          backgroundImage: `url(${mainBackground})`,
          backgroundSize: "100% auto",
          backgroundPosition: "center center",
          backgroundRepeat: "no-repeat",
          py: { xs: 6, md: 8 },
        }}
      >
        {/* Content container with margins */}
        <Box sx={{ px: "6.7%" }}>
          {/* Content limited to 66% of the work area, left aligned */}
          <Box sx={{ width: "66%", maxWidth: "66%" }}>
            <Typography
              variant="h3"
              sx={{
                fontWeight: 800,
                color: "#000000",
                mb: 2,
                fontSize: { xs: "2rem", md: "2.5rem", lg: "3rem" },
              }}
            >
              Welcome to HeartOmicsAtlas
            </Typography>

            <Typography
              variant="body1"
              sx={{
                color: "#000000",
                lineHeight: 1.8,
                mb: 4,
                fontSize: { xs: "0.95rem", md: "1.05rem" },
              }}
            >
              HeartOmicsAtlas is an AI-powered, user-friendly, open-access platform for analyzing
              spatial and single-nucleus (sn-)multiomics data. It features spatial transcriptomics
              and sn-multiomics datasets from fetal heart samples, as well as scRNA-seq data split
              into three subchannels.
            </Typography>

            {/* Two side-by-side buttons */}
            <Stack direction="row" spacing={2}>
              <Button
                variant="contained"
                sx={{
                  borderRadius: "6px",
                  px: 4,
                  py: 1.5,
                  fontSize: "1rem",
                  fontWeight: 600,
                  textTransform: "none",
                  minWidth: 180,
                }}
                onClick={() => navigate("/chat")}
              >
                Chat with AI
              </Button>
              <Button
                variant="outlined"
                sx={{
                  borderRadius: "6px",
                  px: 4,
                  py: 1.5,
                  fontSize: "1rem",
                  fontWeight: 600,
                  textTransform: "none",
                  minWidth: 180,
                  borderColor: "#000000",
                  color: "#000000",
                  "&:hover": {
                    borderColor: "#000000",
                    backgroundColor: "rgba(0, 0, 0, 0.04)",
                  },
                }}
                onClick={() => navigate("/st")}
              >
                Explore Data
              </Button>
            </Stack>
          </Box>
        </Box>
      </Box>

      {/* Second Section - What does HeartOmicsAtlas Include? */}
      <Box sx={{ px: "6.7%", py: 6 }}>
        {/* Section Title */}
        <Typography
          variant="h4"
          sx={{
            fontWeight: 700,
            mb: 4,
            textAlign: "center",
          }}
        >
          <Box component="span" sx={{ color: "#000000" }}>What does Heart</Box>
          <Box component="span" sx={{ color: "#BE1B23" }}>Omics</Box>
          <Box component="span" sx={{ color: "#000000" }}>Atlas Include?</Box>
        </Typography>

        {/* Feature Card */}
        <Box
          sx={{
            backgroundColor: "#FAF8FD",
            borderRadius: "12px",
            boxShadow: "0 4px 16px rgba(0, 0, 0, 0.08)",
            py: 5,
            px: 4,
          }}
        >
          {/* Three items centered horizontally */}
          <Stack
            direction="row"
            justifyContent="center"
            alignItems="center"
            spacing={6}
          >
            {featureItems.map((item, index) => (
              <Stack
                key={index}
                direction="row"
                alignItems="center"
                spacing={2}
              >
                {/* Placeholder image box - replace with actual images */}
                <Box
                  sx={{
                    width: 48,
                    height: 48,
                    backgroundColor: "#e0e0e0",
                    borderRadius: "8px",
                    display: "flex",
                    alignItems: "center",
                    justifyContent: "center",
                  }}
                >
                  {/* Placeholder icon */}
                  <Box
                    component="span"
                    sx={{
                      fontSize: 24,
                      color: "#666",
                    }}
                  >
                    {index + 1}
                  </Box>
                </Box>
                <Typography
                  sx={{
                    color: "#000000",
                    fontWeight: 500,
                    fontSize: "1rem",
                  }}
                >
                  {item.text}
                </Typography>
              </Stack>
            ))}
          </Stack>
        </Box>
      </Box>

      {/* Third Section - How Researchers Use HeartOmicsAtlas */}
      <Box sx={{ px: "6.7%", py: 6 }}>
        {/* Section Title */}
        <Typography
          variant="h4"
          sx={{
            fontWeight: 700,
            mb: 6,
            textAlign: "center",
            color: "#000000",
          }}
        >
          How Researchers Use 
          <Box component="span" sx={{ color: "#000000" }}>What does Heart</Box>
          <Box component="span" sx={{ color: "#BE1B23" }}>Omics</Box>
          <Box component="span" sx={{ color: "#000000" }}>Atlas.</Box>
        </Typography>

        {/* Cards with arrows - spans full working area width */}
        <Stack
          direction="row"
          justifyContent="space-between"
          alignItems="center"
          sx={{ width: "100%", marginBottom: '5vw'}}
        >
          {researcherSteps.map((step, index) => (
            <Box 
              key={index} 
              sx={{ 
                display: "flex", 
                alignItems: "center",
                flex: 1,
              }}
            >
              {/* Card with overlapping red circle */}
              <Box
                sx={{
                  position: "relative",
                  pt: 3, // Space for the overlapping circle
                  flex: 1,
                }}
              >
                {/* Red numbered circle - positioned to overlap top of card */}
                <Box
                  sx={{
                    position: "absolute",
                    top: 0,
                    left: 40, // More to the right within the card
                    width: 48,
                    height: 48,
                    borderRadius: "50%",
                    backgroundColor: "#C30F1A",
                    display: "flex",
                    alignItems: "center",
                    justifyContent: "center",
                    zIndex: 1,
                  }}
                >
                  <Typography
                    sx={{
                      color: "#ffffff",
                      fontWeight: 700,
                      fontSize: "1.25rem",
                    }}
                  >
                    {index + 1}
                  </Typography>
                </Box>

                {/* Card */}
                <Box
                  sx={{
                    minHeight: 340,
                    backgroundColor: "#FAF8FD",
                    borderRadius: "16px",
                    boxShadow: "0 4px 16px rgba(0, 0, 0, 0.08)",
                    p: 3,
                    pt: 5, // Extra top padding for circle overlap area
                    display: "flex",
                    flexDirection: "column",
                    alignItems: "center",
                    justifyContent: "center",
                  }}
                >
                  {/* Placeholder image */}
                  <Box
                    sx={{
                      width: 100,
                      height: 100,
                      backgroundColor: "#e8e4ed",
                      borderRadius: "12px",
                      mb: 3,
                      display: "flex",
                      alignItems: "center",
                      justifyContent: "center",
                    }}
                  >
                    <Typography sx={{ color: "#999", fontSize: 40 }}>
                      📊
                    </Typography>
                  </Box>

                  {/* Title - similar size to section title */}
                  <Typography
                    variant="h5"
                    sx={{
                      fontWeight: 700,
                      color: "#000000",
                      textAlign: "center",
                      mb: 1,
                    }}
                  >
                    {step.title}
                  </Typography>

                  {/* Subtitle */}
                  <Typography
                    sx={{
                      fontSize: "1rem",
                      color: "#666666",
                      textAlign: "center",
                      lineHeight: 1.5,
                    }}
                  >
                    {step.subtitle}
                  </Typography>
                </Box>
              </Box>

              {/* Arrow between cards (not after the last card) */}
              {index < researcherSteps.length - 1 && (
                <Box
                  sx={{
                    mx: 2,
                    display: "flex",
                    alignItems: "center",
                    mt: 3, // Align with card center
                    flexShrink: 0,
                  }}
                >
                  {/* Simple arrow using CSS */}
                  <Box
                    sx={{
                      width: 40,
                      height: 2,
                      backgroundColor: "#C30F1A",
                      position: "relative",
                      "&::after": {
                        content: '""',
                        position: "absolute",
                        right: -1,
                        top: "50%",
                        transform: "translateY(-50%)",
                        width: 0,
                        height: 0,
                        borderTop: "6px solid transparent",
                        borderBottom: "6px solid transparent",
                        borderLeft: "10px solid #C30F1A",
                      },
                    }}
                  />
                </Box>
              )}
            </Box>
          ))}
        </Stack>
      </Box>
    </Box>
  );
}
