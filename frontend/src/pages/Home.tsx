import { useMemo, useState } from "react";
import { useNavigate } from "react-router-dom";
import {
  Box,
  Button,
  Card,
  CardContent,
  Chip,
  InputAdornment,
  Stack,
  TextField,
  Typography,
} from "@mui/material";
import SendIcon from "@mui/icons-material/Send";
import heartLogo from "../assets/heart_logo_1.png";

export default function Home() {
  const navigate = useNavigate();
  const [chatInput, setChatInput] = useState("");

  const suggested = useMemo(
    () => [
      "What is the function of Sinoatrial node in human heart?",
      "Where is the location of Sinoatrial node?",
      "Can you tell me the marker genes for the neuronal cells in the fetal heart?",
    ],
    []
  );

  const goChat = (text: string) => {
    const t = text.trim();
    if (!t) {
      navigate("/chat");
      return;
    }
    navigate("/chat", { state: { chatInput: t } });
  };

  return (
    <Stack spacing={3}>
      {/* Title + description */}
      <Stack spacing={1}>
        <Typography variant="h4" fontWeight={800}>
          Welcome to HeartOmicsAtlas
        </Typography>

        <Typography variant="body1" color="text.secondary" sx={{ lineHeight: 1.7 }}>
          HeartOmicsAtlas is an AI-powered, user-friendly, open-access platform for analyzing
          spatial and single-nucleus (sn-)multiomics data. It features spatial transcriptomics
          and sn-multiomics datasets from fetal heart samples, as well as scRNA-seq data split
          into three subchannels: ACM_VCM_SAN, SAN-PCO and mini-heart. This platform enables
          systematic analysis of transcriptomic signatures across the major cell populations
          within the sinoatrial node niche.
        </Typography>
      </Stack>

      {/* Suggested questions */}
      <Card sx={{ borderRadius: 0 }}>
        <CardContent>
          <Typography variant="h6" fontWeight={700} gutterBottom>
            Suggested questions
          </Typography>

          <Box sx={{ display: "flex", flexWrap: "wrap", gap: 1 }}>
            {suggested.map((q) => (
              <Chip
                key={q}
                label={q}
                onClick={() => goChat(q)}
                clickable
                sx={{
                  borderRadius: 0,
                  bgcolor: "#e9ecef",
                  "&:hover": { bgcolor: "#dde2e6" },
                }}
              />
            ))}
          </Box>
        </CardContent>
      </Card>

      {/* Home chat input */}
      <Card sx={{ borderRadius: 0 }}>
        <CardContent>
          <Typography variant="h6" fontWeight={700} gutterBottom>
            Chat with the AI assistant
          </Typography>

          <Stack direction={{ xs: "column", sm: "row" }} spacing={1.5}>
            <TextField
              fullWidth
              placeholder="Enter text to chat with AI assistant"
              value={chatInput}
              onChange={(e) => setChatInput(e.target.value)}
              onKeyDown={(e) => {
                if (e.key === "Enter") goChat(chatInput);
              }}
              InputProps={{
                startAdornment: (
                  <InputAdornment position="start">
                    <Box
                      component="img"
                      src={heartLogo}
                      alt=""
                      sx={{ width: 28, height: 28, display: "block" }}
                    />
                  </InputAdornment>
                ),
              }}
            />
            <Button
              variant="contained"
              endIcon={<SendIcon />}
              onClick={() => goChat(chatInput)}
              sx={{ borderRadius: 0, minWidth: 140 }}
            >
              Send
            </Button>
          </Stack>

          <Typography variant="body2" color="text.secondary" sx={{ mt: 1 }}>
            Tip: click a suggested question to jump straight into the chat with it pre-filled.
          </Typography>
        </CardContent>
      </Card>

    </Stack>
  );
}
