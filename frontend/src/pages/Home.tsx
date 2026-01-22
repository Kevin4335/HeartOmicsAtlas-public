import { Card, CardContent, Stack, Typography } from "@mui/material";

export default function Home() {
  return (
    <Stack spacing={3}>
      <Typography variant="h4" fontWeight={700}>
        Home
      </Typography>

      <Card>
        <CardContent>
          <Typography variant="h6" gutterBottom>
            Welcome
          </Typography>
          <Typography variant="body1">
            This is a lightweight React + MUI + Router starter for HeartOmics.
          </Typography>
        </CardContent>
      </Card>
    </Stack>
  );
}
