import { useEffect, useMemo, useRef, useState } from "react";
import { useLocation } from "react-router-dom";
import {
  Alert,
  Box,
  Button,
  Card,
  CardContent,
  Container,
  Grid,
  IconButton,
  List,
  ListItem,
  Modal,
  Paper,
  Stack,
  TextField,
  Typography,
} from "@mui/material";
import CloseIcon from "@mui/icons-material/Close";
import SendIcon from "@mui/icons-material/Send";
import { useTheme } from "@mui/material/styles";

// ----------------------
// Types
// ----------------------
type DisplayMessage =
  | { type: "user"; content: string }
  | { type: "text"; content: string }
  | { type: "image"; content: string };

type HistoryItem = { role: "user" | "assistant"; content: string };

type BackendResponse = {
  messages?: Array<{ type: "text" | "image"; content: string }>;
  history?: HistoryItem[];
};

// ----------------------
// Config
// ----------------------
const TEST_MODE = false;

// In your old code, BASEURL was set only in development and blank in production.
// This preserves that behavior but uses Vite env flags.
const BASEURL = import.meta.env.DEV ? "http://128.84.40.121" : "";
const AI_CHAT_URL = `${BASEURL}/chat`;

// LocalStorage keys (keep same so you do not lose existing history)
const LS_OPENAI = "openai-history";
const LS_DISPLAY = "display-history";

export default function AIChat() {
  const theme = useTheme();

  // Router location for initial input (optional)
  const location = useLocation();
  const initialInput = useMemo(() => {
    const s = (location.state as { chatInput?: unknown } | null)?.chatInput;
    return typeof s === "string" ? s : "";
  }, [location.state]);

  // Messages (rendered in UI)
  const [messages, setMessages] = useState<DisplayMessage[]>(() => {
    const stored = localStorage.getItem(LS_DISPLAY);
    if (stored) {
      try {
        return JSON.parse(stored) as DisplayMessage[];
      } catch {
        return [];
      }
    }
    // If user navigated here with pre-filled input, show it as first user message
    return initialInput ? [{ type: "user", content: initialInput }] : [];
  });

  // Input + UI state
  const [input, setInput] = useState<string>("");
  const [waiting, setWaiting] = useState<boolean>(false);

  const [lightboxOpen, setLightboxOpen] = useState<boolean>(false);
  const [lightboxImage, setLightboxImage] = useState<string>("");

  const [error, setError] = useState<string>("");

  // Scroll container ref
  const scrollRef = useRef<HTMLDivElement | null>(null);

  // Auto-scroll when messages update
  useEffect(() => {
    const el = scrollRef.current;
    if (!el) return;
    el.scrollTop = el.scrollHeight;
  }, [messages]);

  // Quick prompts
  const prompts = useMemo(
    () => [
      "What is the function of Sinoatrial node in human heart?",
      "How to identify Sinoatrial node?",
      "Where is the location of Sinoatrial node?",
      "Can you tell me the marker genes for the neuronal cells in the fetal heart?"
    ],
    []
  );

  // Lightbox
  const handleImageClick = (src: string) => {
    setLightboxImage(src);
    setLightboxOpen(true);
  };

  const handleCloseLightbox = () => {
    setLightboxOpen(false);
    setLightboxImage("");
  };

  // Clear history
  const clearHistory = () => {
    if (waiting) return;
    localStorage.setItem(LS_OPENAI, JSON.stringify([]));
    localStorage.setItem(LS_DISPLAY, JSON.stringify([]));
    setMessages([]);
    setError("");
  };

  // Optional: test simulation
  const simulateBackendResponse = async (content: string): Promise<BackendResponse> => {
    await new Promise((r) => setTimeout(r, 700));
    if (content.toLowerCase().includes("image") || content.toLowerCase().includes("png")) {
      // Use any reachable image URL for testing
      return {
        messages: [
          { type: "text", content: "Here is the image you requested:" },
          { type: "image", content: "https://via.placeholder.com/512.png?text=Example" },
        ],
        history: JSON.parse(localStorage.getItem(LS_OPENAI) || "[]"),
      };
    }
    return {
      messages: [{ type: "text", content: `Test response: "${content}"` }],
      history: JSON.parse(localStorage.getItem(LS_OPENAI) || "[]"),
    };
  };

  // Send message
  const sendMessage = async (content: string) => {
    if (waiting) return;

    const trimmed = content.trim();
    if (!trimmed) return;

    setWaiting(true);
    setError("");

    const userMsg: DisplayMessage = { type: "user", content: trimmed };
    const loadingMsg: DisplayMessage = { type: "text", content: "Loading..." };

    // Build openai history
    const openaiHistory: HistoryItem[] = [
      ...(JSON.parse(localStorage.getItem(LS_OPENAI) || "[]") as HistoryItem[]),
      { role: "user", content: trimmed },
    ];

    // Optimistically show user + loading
    const nextDisplay = [...messages, userMsg, loadingMsg];
    setMessages(nextDisplay);
    localStorage.setItem(LS_OPENAI, JSON.stringify(openaiHistory));
    localStorage.setItem(LS_DISPLAY, JSON.stringify(nextDisplay));

    try {
      let data: BackendResponse;

      if (TEST_MODE) {
        data = await simulateBackendResponse(trimmed);
      } else {
        const res = await fetch(AI_CHAT_URL, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify(openaiHistory),
        });

        if (!res.ok) {
          throw new Error(`HTTP ${res.status}`);
        }

        data = (await res.json()) as BackendResponse;
      }

      const processed: DisplayMessage[] = (data.messages || []).map((m) => {
        if (m.type === "image") return { type: "image", content: m.content };
        return { type: "text", content: m.content };
      });

      // Replace the "Loading..." message by rebuilding final list from prior stable state.
      // Use messages (state) + userMsg + processed (not nextDisplay, because nextDisplay included Loading).
      const finalMessages = [...messages, userMsg, ...processed];
      setMessages(finalMessages);
      localStorage.setItem(LS_DISPLAY, JSON.stringify(finalMessages));

      if (data.history) {
        localStorage.setItem(LS_OPENAI, JSON.stringify(data.history));
      }
    } catch (e) {
      const msg = e instanceof Error ? e.message : "Unknown error";
      setError(msg);
      // Remove the loading message, keep the user message
      const finalMessages = [...messages, userMsg, { type: "text", content: `Error: ${msg}` }];
      setMessages(finalMessages);
      localStorage.setItem(LS_DISPLAY, JSON.stringify(finalMessages));
    } finally {
      setWaiting(false);
    }
  };

  // Handle prompt click
  const handlePromptClick = (prompt: string) => {
    setInput("");
    void sendMessage(prompt);
  };

  // Handle send
  const handleSend = () => {
    if (!input.trim()) return;
    const current = input;
    setInput("");
    void sendMessage(current);
  };

  // Send initial input on mount, if provided
  useEffect(() => {
    if (!initialInput) return;
    void sendMessage(initialInput);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  return (
    <Container sx={{ mb: 2 }}>
      <Stack spacing={2}>
        <Typography variant="h4" sx={{ fontWeight: 700, textAlign: "center", mt: 1 }}>
          AI Chat
        </Typography>

        <Box sx={{ display: "flex", justifyContent: "flex-end" }}>
          <Button
            variant="outlined"
            color="error"
            onClick={clearHistory}
            disabled={waiting}
            sx={{ borderRadius: 999, px: 3 }}
          >
            Clear History
          </Button>
        </Box>

        {error ? <Alert severity="error">{error}</Alert> : null}

        <Paper
          elevation={0}
          sx={{
            height: 620,
            borderRadius: 0,
            overflow: "hidden",
            border: "1px solid rgba(0,0,0,0.08)",
            backgroundColor: "#f8f9fa",
            display: "flex",
            flexDirection: "column",
          }}
        >
          {/* Scrollable messages */}
          <Box
            ref={scrollRef}
            sx={{
              flex: 1,
              overflowY: "auto",
              p: 2,
            }}
          >
            {messages.length === 0 ? (
              <Box
                sx={{
                  height: "100%",
                  display: "flex",
                  flexDirection: "column",
                  alignItems: "center",
                  justifyContent: "center",
                  gap: 2,
                  px: 2,
                }}
              >
                <Typography variant="h6" textAlign="center">
                  Try one of these prompts:
                </Typography>

                <Grid container spacing={2} justifyContent="center">
                  {prompts.map((p, idx) => (
                    <Grid item xs={12} sm={6} key={idx}>
                      <Card
                        onClick={() => handlePromptClick(p)}
                        sx={{
                          cursor: "pointer",
                          borderRadius: 0,
                          border: "1px solid rgba(0,0,0,0.08)",
                          "&:hover": { borderColor: "rgba(0,0,0,0.25)" },
                        }}
                      >
                        <CardContent>
                          <Typography>{p}</Typography>
                        </CardContent>
                      </Card>
                    </Grid>
                  ))}
                </Grid>
              </Box>
            ) : (
              <List sx={{ p: 0 }}>
                {messages.map((msg, idx) => {
                  const isUser = msg.type === "user";
                  const isImage = msg.type === "image";

                  return (
                    <ListItem
                      key={idx}
                      sx={{
                        display: "flex",
                        justifyContent: isUser ? "flex-end" : "flex-start",
                        alignItems: "flex-start",
                        px: 0,
                      }}
                    >
                      <Box
                        sx={{
                          maxWidth: "72%",
                          px: 1.5,
                          py: 1.25,
                          borderRadius: 0,
                          border: "1px solid rgba(0,0,0,0.08)",
                          backgroundColor: isUser ? theme.palette.primary.main : "#e9ecef",
                          color: isUser ? "#fff" : theme.palette.text.primary,
                          whiteSpace: "pre-wrap",
                          wordBreak: "break-word",
                        }}
                      >
                        {isImage ? (
                          <Box sx={{ display: "flex", justifyContent: "center" }}>
                            <img
                              src={msg.content}
                              alt="response"
                              style={{
                                maxWidth: 260,
                                maxHeight: 260,
                                cursor: "pointer",
                                display: "block",
                              }}
                              onClick={() => handleImageClick(msg.content)}
                            />
                          </Box>
                        ) : (
                          <Typography variant="body1">{msg.content}</Typography>
                        )}
                      </Box>
                    </ListItem>
                  );
                })}
              </List>
            )}
          </Box>

          {/* Input bar */}
          <Box
            sx={{
              borderTop: "1px solid rgba(0,0,0,0.08)",
              backgroundColor: "#ffffff",
              p: 1.5,
            }}
          >
            <Box
              sx={{
                display: "flex",
                gap: 1,
                alignItems: "center",
              }}
            >
              <TextField
                fullWidth
                placeholder="Ask AI anything..."
                value={input}
                onChange={(e) => setInput(e.target.value)}
                onKeyDown={(e) => {
                  if (e.key === "Enter") handleSend();
                }}
                disabled={waiting}
              />
              <Button
                variant="contained"
                onClick={handleSend}
                disabled={waiting}
                sx={{
                  minWidth: 52,
                  height: 56,
                  borderRadius: 0,
                }}
              >
                <SendIcon />
              </Button>
            </Box>
          </Box>
        </Paper>
      </Stack>

      {/* Lightbox */}
      <Modal
        open={lightboxOpen}
        onClose={handleCloseLightbox}
        sx={{
          display: "flex",
          alignItems: "center",
          justifyContent: "center",
          p: 2,
        }}
      >
        <Box
          sx={{
            position: "relative",
            maxWidth: "92vw",
            maxHeight: "92vh",
            bgcolor: "background.paper",
            borderRadius: 0,
            border: "1px solid rgba(0,0,0,0.15)",
            p: 1,
          }}
        >
          <IconButton
            onClick={handleCloseLightbox}
            sx={{
              position: "absolute",
              right: 8,
              top: 8,
              bgcolor: "rgba(0,0,0,0.55)",
              color: "white",
              borderRadius: 0,
              "&:hover": { bgcolor: "rgba(0,0,0,0.75)" },
            }}
          >
            <CloseIcon />
          </IconButton>

          <img
            src={lightboxImage}
            alt="Full size"
            style={{
              maxWidth: "90vw",
              maxHeight: "90vh",
              display: "block",
            }}
          />
        </Box>
      </Modal>
    </Container>
  );
}
