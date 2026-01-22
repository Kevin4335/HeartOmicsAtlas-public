import { CssBaseline, ThemeProvider } from "@mui/material";
import { createAppTheme } from "./theme";
import AppRoutes from "./routes/AppRoutes";
import AppShell from "./components/AppShell";

export default function App() {
  const theme = createAppTheme("light");

  return (
    <ThemeProvider theme={theme}>
      <CssBaseline />
      <AppShell>
        <AppRoutes />
      </AppShell>
    </ThemeProvider>
  );
}
