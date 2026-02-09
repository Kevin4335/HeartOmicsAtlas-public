#!/bin/bash
# restart_r_servers.sh
# Shut down and restart all R servers in screen sessions

BASE_DIR="/home/ubuntu/HeartOmics"

echo "=========================================="
echo "Restarting R Servers"
echo "=========================================="

# Step 1: Kill existing screen sessions for R servers
echo "Step 1: Stopping existing R server screen sessions..."
for session in spatial multiomics acm_vcm san_pco mini_heart; do
    screen -S "$session" -X quit 2>/dev/null || true
done
# Also kill any screen sessions matching the pattern
screen -ls | grep -E "(spatial|multiomics|acm_vcm|san_pco|mini_heart)" | awk '{print $1}' | while read session; do
    screen -S "$session" -X quit 2>/dev/null || true
done
sleep 2

# Step 2: Kill any processes using the R server ports (backup method)
echo "Step 2: Killing processes on R server ports (9025-9029)..."
for port in 9025 9026 9027 9028 9029; do
    pid=$(lsof -ti:$port 2>/dev/null)
    if [ ! -z "$pid" ]; then
        echo "  Killing process on port $port (PID: $pid)"
        kill -9 $pid 2>/dev/null
    fi
done
sleep 1

# Step 3: Start all R servers in screen sessions
echo "Step 3: Starting R servers in screen sessions..."

cd "$BASE_DIR"

# Spatial server (port 9025)
echo "  Starting Spatial server (port 9025)..."
cd "$BASE_DIR/resources-NEW/spatial_data"
screen -dmS spatial Rscript Spatial_plot_function_new.R
sleep 1

# Multiomics server (port 9026)
echo "  Starting Multiomics server (port 9026)..."
cd "$BASE_DIR/resources-NEW/multi_omics_data"
screen -dmS multiomics Rscript Multiomics_plot_function_new.R
sleep 1

# ACM_VCM_SAN server (port 9027)
echo "  Starting ACM_VCM_SAN server (port 9027)..."
cd "$BASE_DIR/resources-NEW/SAN_ACM_VCM"
screen -dmS acm_vcm Rscript ACM_VCM_SAN_plot_function_new.R
sleep 1

# SAN-PCO server (port 9028)
echo "  Starting SAN-PCO server (port 9028)..."
cd "$BASE_DIR/resources-NEW/SAN-PCO"
screen -dmS san_pco Rscript SAN_PCO_plot_function_new.R
sleep 1

# Mini-heart server (port 9029)
echo "  Starting Mini-heart server (port 9029)..."
cd "$BASE_DIR/resources-NEW/Mini-heart"
screen -dmS mini_heart Rscript mini_heart_plot_function_new.R
sleep 1

# Step 4: Verify servers are running
echo ""
echo "Step 4: Verifying servers are running..."
sleep 2

screen -ls

echo ""
echo "Checking ports..."
for port in 9025 9026 9027 9028 9029; do
    if lsof -ti:$port > /dev/null 2>&1; then
        echo "  ✓ Port $port is active"
    else
        echo "  ✗ Port $port is NOT active"
    fi
done

echo ""
echo "=========================================="
echo "R Server restart complete!"
echo "=========================================="
echo ""
echo "To view a server's output, use:"
echo "  screen -r <session_name>"
echo ""
echo "Available sessions:"
echo "  - spatial (port 9025)"
echo "  - multiomics (port 9026)"
echo "  - acm_vcm (port 9027)"
echo "  - san_pco (port 9028)"
echo "  - mini_heart (port 9029)"
echo ""
echo "To detach from a screen session: Ctrl+A, then D"
