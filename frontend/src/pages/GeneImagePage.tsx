import { useMemo, useState } from "react";
import {
  Alert,
  Box,
  Button,
  Card,
  CardContent,
  CircularProgress,
  Stack,
  TextField,
  Typography,
} from "@mui/material";
import { TransformComponent, TransformWrapper } from "react-zoom-pan-pinch";

type Props = {
  title: string;
  baseUrl: string;
  pathPrefix?: string;
  defaultImageSrc: string;
};

export default function GeneImagePage({
  title,
  baseUrl,
  pathPrefix = "/genes",
  defaultImageSrc,
}: Props) {
  const [gene, setGene] = useState("");
  const [submittedGene, setSubmittedGene] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState("");

  const remoteImgSrc = useMemo(() => {
    const g = submittedGene.trim();
    if (!g) return "";
    return `${baseUrl}${pathPrefix}/${encodeURIComponent(g)}`;
  }, [baseUrl, pathPrefix, submittedGene]);

  const currentSrc = remoteImgSrc || defaultImageSrc;
  const viewerKey = currentSrc;

  function submit() {
    const g = gene.trim();
    if (!g) {
      setError("Enter a gene symbol.");
      return;
    }
    setError("");
    setLoading(true);
    setSubmittedGene(g);
  }

  function clear() {
    setGene("");
    setSubmittedGene("");
    setLoading(false);
    setError("");
  }

  return (
    <Stack spacing={2}>
      <Typography variant="h4" fontWeight={700}>
        {title}
      </Typography>

      {/* Input */}
      <Card sx={{ borderRadius: 0 }}>
        <CardContent>
          <Stack direction={{ xs: "column", sm: "row" }} spacing={1.5}>
            <TextField
              fullWidth
              label="Gene"
              placeholder="Example: INS, PTPRN, GAD1"
              value={gene}
              onChange={(e) => setGene(e.target.value)}
              onKeyDown={(e) => e.key === "Enter" && submit()}
            />
            <Button
              variant="contained"
              onClick={submit}
              sx={{ minWidth: 120, borderRadius: 0 }}
            >
              Load
            </Button>
            <Button
              variant="outlined"
              onClick={clear}
              sx={{ minWidth: 120, borderRadius: 0 }}
            >
              Clear
            </Button>
          </Stack>

          {error && (
            <Box sx={{ mt: 1.5 }}>
              <Alert severity="warning">{error}</Alert>
            </Box>
          )}
        </CardContent>
      </Card>

      {/* Viewer */}
      <Card sx={{ borderRadius: 0 }}>
        <CardContent>
          <Box
            sx={{
              position: "relative",
              height: { xs: 420, md: 600 },
              backgroundColor: "#f8f9fa",
              border: "1px solid rgba(0,0,0,0.08)",
              overflow: "hidden",
            }}
          >
            {/* Loading overlay */}
            {loading && remoteImgSrc && (
              <Box
                sx={{
                  position: "absolute",
                  inset: 0,
                  display: "flex",
                  alignItems: "center",
                  justifyContent: "center",
                  zIndex: 2,
                  backgroundColor: "rgba(248,249,250,0.7)",
                }}
              >
                <Stack spacing={1} alignItems="center">
                  <CircularProgress />
                  <Typography variant="body2" color="text.secondary">
                    Loading {submittedGene}…
                  </Typography>
                </Stack>
              </Box>
            )}

            <TransformWrapper
              key={viewerKey}
              initialScale={1}
              minScale={0.5}
              maxScale={12}
              centerOnInit
              limitToBounds
              wheel={{ step: 0.12 }}
              doubleClick={{ disabled: true }}
              panning={{ velocityDisabled: true }}
            >
              {({ resetTransform, zoomIn, zoomOut }) => (
                <>
                  <TransformComponent
                    wrapperStyle={{ width: "100%", height: "100%" }}
                    contentStyle={{ width: "100%", height: "100%" }}
                  >
                    <img
                      src={currentSrc}
                      alt="viewer"
                      draggable={false}
                      style={{
                        width: "100%",
                        height: "100%",
                        objectFit: "contain",
                        display: "block",
                        userSelect: "none",
                      }}
                      onLoad={() => remoteImgSrc && setLoading(false)}
                      onError={() => {
                        if (remoteImgSrc) {
                          setLoading(false);
                          setError(`Failed to load image for gene: ${submittedGene}`);
                          setSubmittedGene("");
                        } else {
                          setError("Failed to load default image.");
                        }
                      }}
                    />
                  </TransformComponent>

                  {/* Zoom controls: rendered after TransformComponent so they sit on top and receive clicks */}
                  <Box
                    sx={{
                      position: "absolute",
                      top: 8,
                      right: 8,
                      zIndex: 10,
                      display: "flex",
                      gap: 0.5,
                      backgroundColor: "rgba(255,255,255,0.9)",
                      border: "1px solid rgba(0,0,0,0.15)",
                      pointerEvents: "auto",
                    }}
                  >
                    <Button
                      type="button"
                      size="small"
                      onClick={() => zoomOut()}
                      sx={{ borderRadius: 0, minWidth: 36 }}
                    >
                      -
                    </Button>
                    <Button
                      type="button"
                      size="small"
                      onClick={() => zoomIn()}
                      sx={{ borderRadius: 0, minWidth: 36 }}
                    >
                      +
                    </Button>
                    <Button
                      type="button"
                      size="small"
                      onClick={() => resetTransform()}
                      sx={{ borderRadius: 0 }}
                    >
                      Reset
                    </Button>
                  </Box>
                </>
              )}
            </TransformWrapper>
          </Box>

          <Box sx={{ mt: 1.5 }}>
            <Typography variant="body2" color="text.secondary">
              Drag to pan. Scroll to zoom. Reset returns full view.
            </Typography>
            {remoteImgSrc ? (
              <Typography variant="body2" color="text.secondary">
                Source gene: {submittedGene}
              </Typography>
            ) : (
              <Typography variant="body2" color="text.secondary">
                Showing default image.
              </Typography>
            )}
          </Box>
        </CardContent>
      </Card>
    </Stack>
  );
}
