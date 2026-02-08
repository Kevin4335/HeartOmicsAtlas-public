import { useEffect, useMemo, useRef, useState } from "react";
import { useLocation } from "react-router-dom";
import {
  Box,
  Button,
  IconButton,
  InputBase,
  Modal,
  Typography,
} from "@mui/material";
import CloseIcon from "@mui/icons-material/Close";
import FavoriteIcon from "@mui/icons-material/Favorite";
import ArrowOutwardOutlinedIcon from "@mui/icons-material/ArrowOutwardOutlined";

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

const BASEURL = import.meta.env.DEV ? "http://128.84.40.121" : "";
const AI_CHAT_URL = `${BASEURL}/chat`;

const LS_OPENAI = "openai-history";
const LS_DISPLAY = "display-history";

export default function AIChat() {
  const location = useLocation();
  const initialInput = useMemo(() => {
    const s = (location.state as { chatInput?: unknown } | null)?.chatInput;
    return typeof s === "string" ? s : "";
  }, [location.state]);

  const [messages, setMessages] = useState<DisplayMessage[]>(() => {
    const stored = localStorage.getItem(LS_DISPLAY);
    if (stored) {
      try {
        return JSON.parse(stored) as DisplayMessage[];
      } catch {
        return [];
      }
    }
    return initialInput ? [{ type: "user", content: initialInput }] : [];
  });

  const [input, setInput] = useState<string>("");
  const [waiting, setWaiting] = useState<boolean>(false);

  const [lightboxOpen, setLightboxOpen] = useState<boolean>(false);
  const [lightboxImage, setLightboxImage] = useState<string>("");

  const scrollRef = useRef<HTMLDivElement | null>(null);

  useEffect(() => {
    const el = scrollRef.current;
    if (!el) return;
    el.scrollTop = el.scrollHeight;
  }, [messages]);

  const handleImageClick = (src: string) => {
    setLightboxImage(src);
    setLightboxOpen(true);
  };

  const handleCloseLightbox = () => {
    setLightboxOpen(false);
    setLightboxImage("");
  };

  const simulateBackendResponse = async (content: string): Promise<BackendResponse> => {
    await new Promise((r) => setTimeout(r, 700));
    if (content.toLowerCase().includes("image") || content.toLowerCase().includes("png")) {
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

  const sendMessage = async (content: string) => {
    if (waiting) return;

    const trimmed = content.trim();
    if (!trimmed) return;

    setWaiting(true);

    const userMsg: DisplayMessage = { type: "user", content: trimmed };

    const openaiHistory: HistoryItem[] = [
      ...(JSON.parse(localStorage.getItem(LS_OPENAI) || "[]") as HistoryItem[]),
      { role: "user", content: trimmed },
    ];

    const nextDisplay: DisplayMessage[] = [...messages, userMsg];
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

      const processed: DisplayMessage[] = (data.messages || []).map((m): DisplayMessage => {
        if (m.type === "image") return { type: "image", content: m.content };
        return { type: "text", content: m.content };
      });

      const finalMessages: DisplayMessage[] = [...messages, userMsg, ...processed];
      setMessages(finalMessages);
      localStorage.setItem(LS_DISPLAY, JSON.stringify(finalMessages));

      if (data.history) {
        localStorage.setItem(LS_OPENAI, JSON.stringify(data.history));
      }
    } catch (e) {
      const msg = e instanceof Error ? e.message : "Unknown error";
      const errMsg: DisplayMessage = { type: "text", content: `Error: ${msg}` };
      const finalMessages: DisplayMessage[] = [...messages, userMsg, errMsg];
      setMessages(finalMessages);
      localStorage.setItem(LS_DISPLAY, JSON.stringify(finalMessages));
    } finally {
      setWaiting(false);
    }
  };

  const handleSend = () => {
    if (!input.trim()) return;
    const current = input;
    setInput("");
    void sendMessage(current);
  };

  useEffect(() => {
    if (!initialInput) return;
    void sendMessage(initialInput);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  return (
    <Box
      sx={{
        height: "calc(100vh - 8.82vh)",
        display: "flex",
        flexDirection: "column",
        backgroundColor: "#ffffff",
        overflow: "hidden",
      }}
    >
      {/* Main content area - scrollable messages */}
      <Box
        ref={scrollRef}
        sx={{
          flex: 1,
          overflowY: "auto",
          pt: "8.82vh", // Distance from navbar = navbar height
          display: "flex",
          flexDirection: "column",
          alignItems: "center",
        }}
      >
        {/* Messages container - same width as input (42.75%) */}
        <Box
          sx={{
            width: "42.75%",
            display: "flex",
            flexDirection: "column",
            gap: 3,
            pb: 3,
          }}
        >
          {messages.map((msg, idx) => {
            const isUser = msg.type === "user";
            const isImage = msg.type === "image";

            if (isUser) {
              // User message - right aligned, pink bubble
              return (
                <Box
                  key={idx}
                  sx={{
                    display: "flex",
                    justifyContent: "flex-end",
                  }}
                >
                  <Box
                    sx={{
                      maxWidth: "70%",
                      px: 1.5,
                      py: 1,
                      borderRadius: "12px",
                      backgroundColor: "#FFF2F3",
                      color: "#000000",
                    }}
                  >
                    <Typography
                      sx={{
                        fontSize: "0.85rem",
                        lineHeight: 1.5,
                        whiteSpace: "pre-wrap",
                        wordBreak: "break-word",
                      }}
                    >
                      {msg.content}
                    </Typography>
                  </Box>
                </Box>
              );
            } else {
              // AI response - left aligned, no bubble, heart icon
              return (
                <Box
                  key={idx}
                  sx={{
                    display: "flex",
                    flexDirection: "column",
                    alignItems: "flex-start",
                  }}
                >
                  {/* Heart icon */}
                  <FavoriteIcon
                    sx={{
                      color: "#BE1B23",
                      fontSize: 24,
                      mb: 1,
                    }}
                  />
                  {/* Response content below the icon */}
                  <Box sx={{ maxWidth: "85%" }}>
                    {isImage ? (
                      <img
                        src={msg.content}
                        alt="response"
                        style={{
                          maxWidth: "100%",
                          maxHeight: 400,
                          cursor: "pointer",
                          display: "block",
                          borderRadius: "8px",
                        }}
                        onClick={() => handleImageClick(msg.content)}
                      />
                    ) : (
                      <Typography
                        sx={{
                          fontSize: "0.85rem",
                          lineHeight: 1.6,
                          color: "#000000",
                          whiteSpace: "pre-wrap",
                          wordBreak: "break-word",
                        }}
                      >
                        {msg.content}
                      </Typography>
                    )}
                  </Box>
                </Box>
              );
            }
          })}

          {/* Loading indicator */}
          {waiting && (
            <Box
              sx={{
                display: "flex",
                flexDirection: "column",
                alignItems: "flex-start",
              }}
            >
              <FavoriteIcon
                sx={{
                  color: "#BE1B23",
                  fontSize: 24,
                  mb: 1,
                }}
              />
              <Typography
                sx={{
                  fontSize: "0.85rem",
                  color: "#888888",
                }}
              >
                Thinking...
              </Typography>
            </Box>
          )}
        </Box>
      </Box>

      {/* Chat input - fixed at bottom, centered */}
      <Box
        sx={{
          position: "relative",
          display: "flex",
          justifyContent: "center",
          pb: 0,
          pt: 2,
          backgroundColor: "transparent",
        }}
      >
        <Box
          sx={{
            width: "42.75%",
            height: "16.8vh",
            position: "relative",
            backgroundColor: "transparent",
            borderRadius: "12px",
            boxShadow: "0 2px 12px rgba(0, 0, 0, 0.1)",
          }}
        >
          <InputBase
            placeholder="Ask me anything about your data ..."
            value={input}
            onChange={(e) => setInput(e.target.value)}
            onKeyDown={(e) => {
              if (e.key === "Enter" && !e.shiftKey) {
                e.preventDefault();
                handleSend();
              }
            }}
            disabled={waiting}
            multiline
            sx={{
              width: "100%",
              height: "100%",
              px: 2.5,
              py: 2,
              fontSize: "1rem",
              alignItems: "flex-start",
              "& .MuiInputBase-input": {
                height: "100% !important",
                overflow: "auto !important",
              },
            }}
          />
          {/* Send button */}
          <IconButton
            onClick={handleSend}
            disabled={waiting}
            sx={{
              position: "absolute",
              bottom: 12,
              right: 12,
              backgroundColor: "#C30F1A",
              width: 36,
              height: 36,
              "&:hover": {
                backgroundColor: "#a00d16",
              },
              "&:disabled": {
                backgroundColor: "#cccccc",
              },
            }}
          >
            <ArrowOutwardOutlinedIcon sx={{ color: "#ffffff", fontSize: 20 }} />
          </IconButton>
        </Box>

        {/* "What to ask" button - pill shaped, in right margin, aligned with send button */}
        <Button
          variant="outlined"
          sx={{
            position: "absolute",
            left: "calc(50% + 42.75%/2 + 16px)", // Just to the right of input box
            bottom: 12, // Aligned with send button (12px from input box bottom)
            height: 36,
            color: "#C30F1A",
            borderColor: "#C30F1A",
            backgroundColor: "#ffffff",
            borderRadius: "50px",
            textTransform: "none",
            fontWeight: 500,
            px: 2,
            whiteSpace: "nowrap",
            "&:hover": {
              borderColor: "#C30F1A",
              backgroundColor: "#fff5f5",
            },
            "&:focus": {
              outline: "none",
            },
            "&:focus-visible": {
              outline: "none",
            },
          }}
        >
          What to ask
        </Button>
      </Box>

      {/* Simple footer */}
      <Box
        component="footer"
        sx={{
          py: 2,
          textAlign: "center",
          backgroundColor: "#ffffff",
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

      {/* Lightbox for images */}
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
            borderRadius: "8px",
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
    </Box>
  );
}
