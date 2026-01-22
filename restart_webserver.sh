#!/bin/bash
# restart_webserver.sh
# Shut down and restart the Python webserver in a screen session

BASE_DIR="/home/ubuntu/HeartOmics"

echo "=========================================="
echo "Restarting Web Server"
echo "=========================================="

# Step 1: Kill existing webserver screen session
echo "Step 1: Stopping existing webserver screen session..."
screen -S webserver -X quit 2>/dev/null || true
# Also kill any screen sessions matching the pattern
screen -ls | grep -E "(webserver|server|python.*server)" | awk '{print $1}' | while read session; do
    screen -S "$session" -X quit 2>/dev/null || true
done
sleep 2

# Step 2: Kill any processes using port 80 (backup method)
echo "Step 2: Killing processes on port 80..."
pid=$(lsof -ti:80 2>/dev/null)
if [ ! -z "$pid" ]; then
    echo "  Killing process on port 80 (PID: $pid)"
    kill -9 $pid 2>/dev/null
    sleep 1
fi

# Step 3: Start webserver in screen session
echo "Step 3: Starting webserver in screen session..."
cd "$BASE_DIR"
screen -dmS webserver sudo python3 server.py
sleep 2

# Step 4: Verify server is running
echo ""
echo "Step 4: Verifying server is running..."
sleep 1

screen -ls | grep webserver

echo ""
echo "Checking port 80..."
if lsof -ti:80 > /dev/null 2>&1; then
    echo "  ✓ Port 80 is active"
else
    echo "  ✗ Port 80 is NOT active"
fi

echo ""
echo "=========================================="
echo "Web Server restart complete!"
echo "=========================================="
echo ""
echo "To view the server's output, use:"
echo "  screen -r webserver"
echo ""
echo "To detach from the screen session: Ctrl+A, then D"
