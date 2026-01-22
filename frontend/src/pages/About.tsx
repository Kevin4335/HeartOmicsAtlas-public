import { Card, CardContent, Stack, Typography } from "@mui/material";

export default function About() {
  return (
    <Stack spacing={3}>
      <Typography variant="h4" fontWeight={700}>
        About
      </Typography>

      <Card>
        <CardContent>
          <Typography variant="body1">
            Put app info here. You can add routes under src/routes/AppRoutes.tsx.
          </Typography>
        </CardContent>
      </Card>
    </Stack>
  );
}
