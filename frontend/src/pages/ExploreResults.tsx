import { useState, useEffect } from "react";
import { useSearchParams, useNavigate } from "react-router-dom";
import { Box, InputBase, Typography, Collapse, IconButton, keyframes } from "@mui/material";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import ExpandLessIcon from "@mui/icons-material/ExpandLess";
import AddIcon from "@mui/icons-material/Add";
import RemoveIcon from "@mui/icons-material/Remove";
import { TransformWrapper, TransformComponent } from "react-zoom-pan-pinch";

// Animated gradient keyframes - wave-like motion between bottom-left and top-right
const gradientAnimation = keyframes`
  0% {
    background-position: 0% 100%;
  }
  15% {
    background-position: 25% 85%;
  }
  30% {
    background-position: 60% 55%;
  }
  45% {
    background-position: 100% 0%;
  }
  55% {
    background-position: 85% 25%;
  }
  70% {
    background-position: 40% 70%;
  }
  85% {
    background-position: 15% 90%;
  }
  100% {
    background-position: 0% 100%;
  }
`;

// Loading dots animation
const dotsAnimation = keyframes`
  0%, 20% {
    content: ".";
  }
  40% {
    content: "..";
  }
  60% {
    content: "...";
  }
  80% {
    content: "..";
  }
  100% {
    content: ".";
  }
`;

// Images from Home page (for sidebar)
import multiomicsImg from "../assets/multiomics.png";
import spatialImg from "../assets/spatialtrans.png";
import scrnaImg from "../assets/scrnaseq.png";

// Default plot images
import multiomicsDefault from "../assets/imgs/Multiomics_default_plots.png";
import spatialDefault from "../assets/imgs/Spatial_default_plots.png";
import acmVcmSanDefault from "../assets/imgs/ACM_VCM_SAN_default_plots.png";
import sanPcoDefault from "../assets/imgs/SAN_PCO_default_plots.png";
import miniHeartDefault from "../assets/imgs/mini_heart_default_plots.png";

// Map type IDs to default images
const DEFAULT_IMAGES: Record<string, string> = {
  multiomics: multiomicsDefault,
  spatial: spatialDefault,
  acm_vcm_san: acmVcmSanDefault,
  san_pco: sanPcoDefault,
  mini_heart: miniHeartDefault,
};

// Backend URLs for each data type
const BACKEND_URLS = {
  multiomics: "http://128.84.41.80:9026",
  spatial: "http://128.84.41.80:9025",
  acm_vcm_san: "http://128.84.41.80:9027",
  san_pco: "http://128.84.41.80:9028",
  mini_heart: "http://128.84.41.80:9029",
};

// Main data types
const dataTypes = [
  { id: "multiomics", label: "Multiomics", image: multiomicsImg },
  { id: "spatial", label: "Spatial Transcriptomics", image: spatialImg },
  { id: "scrna", label: "scRNA-seq", image: scrnaImg, expandable: true },
];

// scRNA-seq subtypes (no logos)
const scrnaSubtypes = [
  { id: "acm_vcm_san", label: "ACM_VCM_SAN" },
  { id: "san_pco", label: "SAN-PCO" },
  { id: "mini_heart", label: "Mini-heart" },
];

export default function ExploreResults() {
  const [searchParams] = useSearchParams();
  const navigate = useNavigate();
  const gene = searchParams.get("gene") || "";
  const typeParam = searchParams.get("type") || "";
  
  const [searchValue, setSearchValue] = useState(gene);
  const [scrnaExpanded, setScrnaExpanded] = useState(
    // Auto-expand if a scRNA subtype is selected
    scrnaSubtypes.some(s => s.id === typeParam)
  );
  const [selectedType, setSelectedType] = useState<string | null>(typeParam || null);
  const [isLoading, setIsLoading] = useState(false);
  const [hasError, setHasError] = useState(false);

  // Sync state with URL params when navigating back to this page
  useEffect(() => {
    setSearchValue(gene);
    setSelectedType(typeParam || null);
    if (scrnaSubtypes.some(s => s.id === typeParam)) {
      setScrnaExpanded(true);
    }
  }, [gene, typeParam]);

  // Compute the image URL to display
  const getImageUrl = () => {
    if (!selectedType) return null;
    
    // If no gene, show default image
    if (!gene) {
      return DEFAULT_IMAGES[selectedType] || null;
    }
    
    // If gene is present, try to load from backend
    const baseUrl = BACKEND_URLS[selectedType as keyof typeof BACKEND_URLS];
    if (!baseUrl) return null;
    return `${baseUrl}/genes/${encodeURIComponent(gene)}`;
  };

  const imageUrl = getImageUrl();

  // Reset loading/error state when image URL changes
  useEffect(() => {
    if (imageUrl && gene) {
      setIsLoading(true);
      setHasError(false);
    } else {
      setIsLoading(false);
      setHasError(false);
    }
  }, [imageUrl, gene]);

  const handleSearch = () => {
    const query = searchValue.trim();
    if (query) {
      const params = new URLSearchParams();
      params.set("gene", query);
      if (selectedType) params.set("type", selectedType);
      navigate(`/explore/results?${params.toString()}`);
    }
  };

  const handleClear = () => {
    setSearchValue("");
    setSelectedType(null);
    // Clear the gene from URL as well
    navigate("/explore/results");
  };

  const hasSearchValue = searchValue.trim().length > 0;

  const updateSelectedType = (typeId: string) => {
    setSelectedType(typeId);
    // Persist selection in URL
    const params = new URLSearchParams(searchParams);
    params.set("type", typeId);
    navigate(`/explore/results?${params.toString()}`, { replace: true });
  };

  const handleTypeClick = (typeId: string, expandable?: boolean) => {
    if (expandable) {
      setScrnaExpanded(!scrnaExpanded);
    } else {
      updateSelectedType(typeId);
    }
  };

  const handleSubtypeClick = (subtypeId: string) => {
    updateSelectedType(subtypeId);
  };

  return (
    <Box
      sx={{
        minHeight: "100vh",
        display: "flex",
        flexDirection: "column",
        backgroundColor: "#ffffff",
      }}
    >
      {/* Main content area */}
      <Box
        sx={{
          flex: 1,
          display: "flex",
          px: "14%",
          pt: "6vh",
          pb: 4,
          gap: "3.73%",
        }}
      >
        {/* Left side - 18.1% of screen width */}
        <Box
          sx={{
            width: "18.1vw",
            flexShrink: 0,
            display: "flex",
            flexDirection: "column",
          }}
        >
          {/* Search bar - very round, no internal margins */}
          <Box
            sx={{
              display: "flex",
              alignItems: "center",
              backgroundColor: "#ffffff",
              borderRadius: "50px",
              boxShadow: "0 2px 8px rgba(0, 0, 0, 0.1)",
              overflow: "hidden",
              mb: 3,
              border: hasSearchValue ? "none" : "2px solid #C30F1A",
            }}
          >
            <InputBase
              placeholder="Search gene..."
              value={searchValue}
              onChange={(e) => setSearchValue(e.target.value)}
              onKeyDown={(e) => {
                if (e.key === "Enter") handleSearch();
              }}
              sx={{
                flex: 1,
                px: 2,
                py: 1,
                fontSize: "0.9rem",
              }}
            />
            <Box
              onClick={hasSearchValue ? handleClear : handleSearch}
              sx={{
                px: 1,
                mr: 0.5,
                display: "flex",
                alignItems: "center",
                cursor: "pointer",
                borderRadius: "50px",
                "&:hover": { backgroundColor: "rgba(195, 15, 26, 0.08)" },
              }}
            >
              <Typography
                sx={{
                  color: "#C30F1A",
                  fontSize: "0.85rem",
                  fontWeight: 600,
                }}
              >
                {hasSearchValue ? "Clear" : "Search"}
              </Typography>
            </Box>
          </Box>

          {/* Data type items */}
          <Box sx={{ flex: 1 }}>
            {dataTypes.map((type) => (
              <Box key={type.id}>
                {/* Main item */}
                <Box
                  onClick={() => handleTypeClick(type.id, type.expandable)}
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    py: 1.5,
                    cursor: "pointer",
                  }}
                >
                  <Box
                    component="img"
                    src={type.image}
                    alt={type.label}
                    sx={{ width: 32, height: 32, mr: 1.5, objectFit: "contain" }}
                  />
                  <Typography
                    sx={{
                      flex: 1,
                      fontSize: "0.9rem",
                      fontWeight: 400,
                      color: selectedType === type.id ? "#C30F1A" : "#000000",
                    }}
                  >
                    {type.label}
                  </Typography>
                  {type.expandable && (
                    scrnaExpanded ? (
                      <ExpandLessIcon sx={{ fontSize: 20, color: "#666" }} />
                    ) : (
                      <ExpandMoreIcon sx={{ fontSize: 20, color: "#666" }} />
                    )
                  )}
                </Box>

                {/* scRNA-seq expandable subtypes */}
                {type.expandable && (
                  <Collapse in={scrnaExpanded}>
                    <Box sx={{ pl: 2 }}>
                      {scrnaSubtypes.map((subtype) => (
                        <Box
                          key={subtype.id}
                          onClick={() => handleSubtypeClick(subtype.id)}
                          sx={{
                            py: 1,
                            pl: 2,
                            cursor: "pointer",
                          }}
                        >
                          <Typography
                            sx={{
                              fontSize: "0.85rem",
                              fontWeight: 400,
                              color: selectedType === subtype.id ? "#C30F1A" : "#000000",
                            }}
                          >
                            {subtype.label}
                          </Typography>
                        </Box>
                      ))}
                    </Box>
                  </Collapse>
                )}
              </Box>
            ))}
          </Box>

          {/* Instructions at bottom */}
          <Box sx={{ mt: "auto", pt: 4 }}>
            <Typography sx={{ fontSize: "0.8rem", color: "#000000", lineHeight: 1.6 }}>
              Instructions:
            </Typography>
            <Typography sx={{ fontSize: "0.8rem", color: "#000000", lineHeight: 1.6 }}>
              Drag to pan.
            </Typography>
            <Typography sx={{ fontSize: "0.8rem", color: "#000000", lineHeight: 1.6 }}>
              Scroll to zoom.
            </Typography>
            <Typography sx={{ fontSize: "0.8rem", color: "#000000", lineHeight: 1.6 }}>
              Reset returns full view.
            </Typography>
          </Box>
        </Box>

        {/* Right side - 50% of screen width, fixed height */}
        <Box
          sx={{
            width: "50vw",
            height: "70vh",
            flexShrink: 0,
            borderRadius: "12px",
            overflow: "hidden",
            position: "relative",
            // Static color when not loading, animated gradient when loading
            ...(isLoading
              ? {
                  background: "linear-gradient(-45deg, #FFD7CE 0%, #FFF2F3 25%, #FFD7CE 50%, #FFF2F3 75%, #FFD7CE 100%)",
                  backgroundSize: "400% 400%",
                  animation: `${gradientAnimation} 4s ease-in-out infinite`,
                }
              : {
                  backgroundColor: "#FFF2F3",
                }),
          }}
        >
          {isLoading ? (
            // Loading state with animated dots
            <Box sx={{ display: "flex", alignItems: "center", justifyContent: "center", width: "100%", height: "100%" }}>
              <Typography
                sx={{
                  color: "#C30F1A",
                  fontSize: "1.5rem",
                  fontWeight: 500,
                  "&::after": {
                    content: '"."',
                    animation: `${dotsAnimation} 1.5s steps(1, end) infinite`,
                  },
                }}
              >
                Loading
              </Typography>
            </Box>
          ) : hasError ? (
            // Error state
            <Box sx={{ display: "flex", alignItems: "center", justifyContent: "center", width: "100%", height: "100%" }}>
              <Typography
                sx={{
                  color: "#C30F1A",
                  fontSize: "1.5rem",
                  fontWeight: 500,
                }}
              >
                Error
              </Typography>
            </Box>
          ) : imageUrl ? (
            // Show image in zoom/pan viewer
            <TransformWrapper
              key={imageUrl}
              initialScale={0.6}
              minScale={0.3}
              maxScale={12}
              centerOnInit
              wheel={{ step: 0.12 }}
              doubleClick={{ disabled: true }}
              panning={{ velocityDisabled: true }}
            >
              {({ resetTransform, zoomIn, zoomOut }) => (
                <>
                  <TransformComponent
                    wrapperStyle={{ width: "100%", height: "100%" }}
                    contentStyle={{ 
                      width: "100%", 
                      height: "100%",
                      display: "flex",
                      alignItems: "center",
                      justifyContent: "center",
                    }}
                  >
                    <img
                      src={imageUrl}
                      alt={`${selectedType} visualization`}
                      draggable={false}
                      style={{
                        maxWidth: "100%",
                        maxHeight: "100%",
                        objectFit: "contain",
                        display: "block",
                        userSelect: "none",
                      }}
                      onLoad={() => setIsLoading(false)}
                      onError={() => {
                        setIsLoading(false);
                        setHasError(true);
                      }}
                    />
                  </TransformComponent>

                  {/* Zoom controls - bottom right, red outline buttons in white circles */}
                  <Box
                    sx={{
                      position: "absolute",
                      bottom: 16,
                      right: 16,
                      zIndex: 10,
                      display: "flex",
                      gap: 1,
                      pointerEvents: "auto",
                    }}
                  >
                    <IconButton
                      onClick={() => zoomOut()}
                      sx={{
                        width: 36,
                        height: 36,
                        backgroundColor: "#ffffff",
                        border: "2px solid #C30F1A",
                        color: "#C30F1A",
                        "&:hover": {
                          backgroundColor: "#fff5f5",
                        },
                      }}
                    >
                      <RemoveIcon fontSize="small" />
                    </IconButton>
                    <IconButton
                      onClick={() => zoomIn()}
                      sx={{
                        width: 36,
                        height: 36,
                        backgroundColor: "#ffffff",
                        border: "2px solid #C30F1A",
                        color: "#C30F1A",
                        "&:hover": {
                          backgroundColor: "#fff5f5",
                        },
                      }}
                    >
                      <AddIcon fontSize="small" />
                    </IconButton>
                    <Box
                      onClick={() => resetTransform()}
                      sx={{
                        height: 36,
                        px: 1.5,
                        display: "flex",
                        alignItems: "center",
                        justifyContent: "center",
                        backgroundColor: "#ffffff",
                        border: "2px solid #C30F1A",
                        borderRadius: "50px",
                        color: "#C30F1A",
                        fontSize: "0.85rem",
                        fontWeight: 600,
                        cursor: "pointer",
                        "&:hover": {
                          backgroundColor: "#fff5f5",
                        },
                      }}
                    >
                      Reset
                    </Box>
                  </Box>
                </>
              )}
            </TransformWrapper>
          ) : null}
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
