import { createTheme } from "@mui/material/styles";

export function createAppTheme(mode: "light" | "dark") {
  // We ignore dark mode for now because your palette is explicitly light.
  // If you later want real dark mode, we define a second palette.
  // Right now consistency > pretending dark mode exists.

  return createTheme({
    palette: {
      mode: "light",

      primary: {
        main: "#8B0000",
        dark: "#a52a2a",
        contrastText: "#ffffff",
      },

      background: {
        default: "#f8f9fa",
        paper: "#f8f9fa", 
      },

      text: {
        primary: "#263238",
        secondary: "#495057",
      },
    },

    shape: {
      borderRadius: 12,
    },

    typography: {
      fontFamily: [
        "Inter",
        "system-ui",
        "-apple-system",
        "Segoe UI",
        "Roboto",
        "Helvetica",
        "Arial",
        "sans-serif",
      ].join(","),
      h4: {
        fontWeight: 700,
      },
      h6: {
        fontWeight: 600,
      },
    },

    components: {
      MuiCard: {
        styleOverrides: {
          root: {
            backgroundColor: "#e9ecef",
            borderRadius: 16,
          },
        },
      },

      MuiAppBar: {
        styleOverrides: {
          root: {
            backgroundColor: "#8B0000",
          },
        },
      },

      MuiButton: {
        styleOverrides: {
          root: {
            textTransform: "none",
            fontWeight: 600,
          },
        },
      },
    },
  });
}
