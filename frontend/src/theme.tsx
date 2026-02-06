import { createTheme } from "@mui/material/styles";

export function createAppTheme(_mode: "light" | "dark") {
  // We ignore dark mode for now because your palette is explicitly light.
  // If you later want real dark mode, we define a second palette.
  // Right now consistency > pretending dark mode exists.

  return createTheme({
    palette: {
      mode: "light",

      primary: {
        main: "#C30F1A",
        dark: "#9a0c15",
        contrastText: "#ffffff",
      },

      background: {
        default: "#ffffff",
        paper: "#ffffff", 
      },

      text: {
        primary: "#000000",
        secondary: "#333333",
      },
    },

    shape: {
      borderRadius: 12,
    },

    typography: {
      fontFamily: [
        "Effra Trial",
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
            backgroundColor: "#f5f5f5",
            borderRadius: 16,
          },
        },
      },

      MuiAppBar: {
        styleOverrides: {
          root: {
            backgroundColor: "#ffffff",
            color: "#000000",
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
