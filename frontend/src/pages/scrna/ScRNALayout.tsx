import { Outlet, useLocation, useNavigate } from "react-router-dom";
import { Box, Card, CardContent, Stack, Tabs, Tab, Typography } from "@mui/material";

const tabs = [
  { label: "ACM_VCM_SAN", path: "/scrna/acm-vcm-san" },
  { label: "SAN-PCO", path: "/scrna/san-pco" },
  { label: "Mini-heart", path: "/scrna/mini-heart" },
];

export default function ScRNALayout() {
  const location = useLocation();
  const navigate = useNavigate();

  const current = tabs.findIndex((t) => location.pathname.startsWith(t.path));
  const value = current >= 0 ? current : 0;

  return (
    <Stack spacing={2}>
      <Typography variant="h4" fontWeight={700}>
        scRNA
      </Typography>

      <Card>
        <CardContent>
          <Typography variant="h6" gutterBottom>
            Select Cell Type scRNA
          </Typography>
          <Typography color="text.secondary">
            The scRNA channel is split into three subchannels: ACM_VCM_SAN, SAN-PCO, and Mini-heart.
          </Typography>

          <Box sx={{ mt: 2 }}>
            <Tabs
              value={value}
              onChange={(_, newValue: number) => navigate(tabs[newValue].path)}
              variant="scrollable"
              scrollButtons="auto"
              sx={{
                "& .MuiTab-root": { textTransform: "none", fontWeight: 600 },
              }}
            >
              {tabs.map((t) => (
                <Tab key={t.path} label={t.label} />
              ))}
            </Tabs>
          </Box>
        </CardContent>
      </Card>

      <Outlet />
    </Stack>
  );
}
