import { Button, Stack, Typography } from "@mui/material";
import { Link as RouterLink } from "react-router-dom";

export default function NotFound() {
  return (
    <Stack spacing={2}>
      <Typography variant="h4" fontWeight={700}>
        404
      </Typography>
      <Typography variant="body1">That page does not exist.</Typography>
      <Button variant="contained" component={RouterLink} to="/">
        Go home
      </Button>
    </Stack>
  );
}
