import { useNavigate } from "react-router-dom";
import {
  Box,
  Button,
  Stack,
  Typography,
} from "@mui/material";
import DoubleArrowIcon from "@mui/icons-material/DoubleArrow";
import mainBackground from "../assets/main_background.svg";
import logo from "../assets/heart_logo_1.svg";
import scrnaseq from "../assets/scrnaseq.png";
import spatialtrans from "../assets/spatialtrans.png";
import multiomics from "../assets/multiomics.png";
import consultai from "../assets/consultai.png";
import comparemodalities from "../assets/comparemodalities.png";
import selectagene from "../assets/selectagene.png";

//  3 feature items 
const featureItems = [
  { image: scrnaseq, text: "scRNA-seq" },
  { image: spatialtrans, text: "Spatial Transcriptomics" },
  { image: multiomics, text: "scMultiomics" },
];

// Placeholder for the 3 "How Researchers Use" cards
const researcherSteps = [
  {
    image: selectagene,
    title: "Search a gene",
    subtitle: "View expression across cell types and space",
  },
  {
    image: comparemodalities,
    title: "Compare modalities",
    subtitle: "View expression across scRNA, spatial, and scMultiomics",
  },
  {
    image: consultai,
    title: "Ask questions with AI",
    subtitle: "Interpret data patterns, support hypothesis generation",
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
        <Box sx={{ px: "7.7%" }}>
          {/* Flex row */}
          <Box
            sx={{
              display: "flex",
              alignItems: "center",
              justifyContent: "space-between",
              gap: 4,
            }}
          >
            {/* LEFT SIDE — text */}
            <Box sx={{ width: { xs: "100%", md: "75%" } }}>
              <Typography
                variant="h3"
                sx={{
                  fontWeight: 600,
                  color: "#000000",
                  mb: 2,
                  fontSize: { xs: "1.5rem", md: "2em", lg: "3.00rem" },
                }}
              >
                Mapping the molecular architecture of the human sinoatrial node
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
                An interactive atlas for spatial single cell, and scMultiomics data Powered by AI assistant
              </Typography>

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
                    minWidth: 200,
                    background: "linear-gradient(90deg, #C30F1A 0%, #FD441E 100%)",
                    boxShadow: "0 4px 12px rgba(195, 15, 26, 0.3)",
                    border: "none",
                    "&:hover": {
                      background: "linear-gradient(90deg, #a00d16 0%, #e03d1a 100%)",
                      boxShadow: "0 6px 16px rgba(195, 15, 26, 0.4)",
                    },
                  }}
                  onClick={() => navigate("/explore")}
                >
                  Start exploring the atlas
                </Button>
                <Button
                  variant="contained"
                  sx={{
                    borderRadius: "6px",
                    px: 4,
                    py: 1.5,
                    fontSize: "1rem",
                    fontWeight: 600,
                    textTransform: "none",
                    minWidth: 200,
                    backgroundColor: "#ffffff",
                    color: "#000000",
                    boxShadow: "0 4px 12px rgba(0, 0, 0, 0.1)",
                    border: "none",
                    "&:hover": {
                      backgroundColor: "#fafafa",
                      boxShadow: "0 6px 16px rgba(0, 0, 0, 0.15)",
                    },
                  }}
                  onClick={() => navigate("/chat")}
                >
                  Ask the AI assistant
                </Button>
              </Stack>
            </Box>

            {/* RIGHT SIDE — BIG LOGO */}
            <Box
              component="img"
              src={logo}
              alt="HeartOmicsAtlas Logo"
              sx={{
                width: { xs: 0, lg: 204},
                height: "auto",
                objectFit: "contain",
                marginRight: '10vw',
                display: { xs: "none", md: "block" }, // hide on mobile
              }}
            />
          </Box>
        </Box>

      </Box>

      {/* Second Section - What does HeartOmicsAtlas Include? */}
      <Box sx={{ px: "7.7%", py: 6 }}>
        {/* Section Title */}
        <Typography
          variant="h4"
          sx={{
            fontWeight: 400,
            mb: 4,
            textAlign: "center",
          }}
        >
          What data does{" "}
          <Box component="span" sx={{ color: "#000000", fontWeight: 600 }}>Heart</Box>
          <Box component="span" sx={{ color: "#BE1B23", fontWeight: 600 }}>Omics</Box>
          <Box component="span" sx={{ color: "#000000", fontWeight: 600 }}>Atlas</Box>
          {" "}Include?
        </Typography>

        {/* Feature Card */}
        <Box
          sx={{
            backgroundColor: "#FAF8FD",
            borderRadius: "12px",
            boxShadow: "0 4px 16px rgba(0, 0, 0, 0.08)",
            py: 3,
            px: 6,
          }}
        >
          {/* Three items spaced across full width */}
          <Stack
            direction="row"
            justifyContent="space-between"
            alignItems="center"
          >
            {featureItems.map((item, index) => (
              <Stack
                key={index}
                direction="row"
                alignItems="center"
                spacing={1.5}
              >
                {/* Placeholder image box - replace with actual images */}
                <Box
                  component="img"
                  src={item.image}
                  alt={item.text}
                  sx={{
                    width: 78,
                    height: 78,
                    borderRadius: "8px",
                    objectFit: "contain",
                  }}
                />
                <Typography
                  sx={{
                    color: "#000000",
                    fontWeight: 500,
                    fontSize: "1.15rem",
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
      <Box sx={{ px: "7.7%", py: 2 }}>
        {/* Section Title */}
        <Typography
          variant="h4"
          sx={{
            fontWeight: 400,
            mb: 4,
            textAlign: "center",
            color: "#000000",
          }}
        >
          How researchers use{" "}
          <Box component="span" sx={{ color: "#000000", fontWeight: 600 }}>Heart</Box>
          <Box component="span" sx={{ color: "#BE1B23", fontWeight: 600 }}>Omics</Box>
          <Box component="span" sx={{ color: "#000000", fontWeight: 600 }}>Atlas</Box>
        </Typography>

        {/* Cards with arrows - 5 flex items: card, arrow, card, arrow, card */}
        <Box
          sx={{
            display: "flex",
            flexDirection: "row",
            alignItems: "center",
            width: "100%",
            marginBottom: "5vw",
          }}
        >
          {researcherSteps.map((step, index) => (
            <Box 
              key={index} 
              sx={{ 
                display: "contents", // Allows children to participate directly in parent flex
              }}
            >
              {/* Card with overlapping red circle - flex: 3 (larger portion) */}
              <Box
                sx={{
                  position: "relative",
                  pt: 3, // Space for the overlapping circle
                  flex: 3, // Each card takes 3 parts
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
                      fontWeight: 600,
                      fontSize: "1.25rem",
                    }}
                  >
                    {index + 1}
                  </Typography>
                </Box>

                {/* Card - slightly taller than wide */}
                <Box
                  sx={{
                    aspectRatio: "1 / 1.1",
                    backgroundColor: "#FAF8FD",
                    borderRadius: "16px",
                    boxShadow: "0 4px 16px rgba(0, 0, 0, 0.08)",
                    px: 4,
                    py: 1.5,
                    display: "flex",
                    flexDirection: "column",
                    alignItems: "center",
                    justifyContent: "center",
                  }}
                >
                  {/* Image - larger, takes more of the card area */}
                  <Box
                    component="img"
                    src={step.image}
                    alt={step.title}
                    sx={{
                      width: "70%",
                      height: "auto",
                      objectFit: "contain",
                      mb: 1.5,
                    }}
                  />

                  {/* Title - similar size to section title */}
                  <Typography
                    variant="h5"
                    sx={{
                      fontWeight: 600,
                      color: "#000000",
                      textAlign: "center",
                      mb: 2,
                    }}
                  >
                    {step.title}
                  </Typography>

                  {/* Subtitle */}
                  <Typography
                    sx={{
                      fontSize: "1.35rem",
                      color: "#666666",
                      textAlign: "center",
                      lineHeight: 1.6,
                    }}
                  >
                    {step.subtitle}
                  </Typography>
                </Box>
              </Box>

              {/* Arrow between cards (not after the last card) - flex: 1 (smaller portion) */}
              {index < researcherSteps.length - 1 && (
                <Box
                  sx={{
                    flex: 1, // Each arrow takes 1 part
                    display: "flex",
                    alignItems: "center",
                    justifyContent: "center",
                    mt: 3,
                  }}
                >
                  <DoubleArrowIcon
                    sx={{
                      color: "#FFB797",
                      fontSize: 80,
                      transform: "scaleY(2)", // Stretch vertically to make it very tall
                    }}
                  />
                </Box>
              )}
            </Box>
          ))}
        </Box>
      </Box>
    </Box>
  );
}
