#!/bin/bash
# Quick setup script for ATGCFLOW Backend Service

set -e

echo "🚀 ATGCFLOW Backend Service Setup"
echo "=================================="
echo ""

# Get current user
CURRENT_USER=$(whoami)
BACKEND_DIR="/media/drprabudh/m3/Nextflow-Script/WholeExome/backend"

echo "📋 Current user: $CURRENT_USER"
echo "📁 Backend directory: $BACKEND_DIR"
echo ""

# Step 1: Create logs directory
echo "📝 Step 1: Creating logs directory..."
mkdir -p "$BACKEND_DIR/logs"
echo "✅ Logs directory created"
echo ""

# Step 2: Update service file with current user
echo "📝 Step 2: Updating service file..."
sed "s/%USER%/$CURRENT_USER/g" "$BACKEND_DIR/atgcflow-backend.service" > /tmp/atgcflow-backend.service
echo "✅ Service file updated"
echo ""

# Step 3: Install systemd service
echo "📝 Step 3: Installing systemd service..."
echo "   (This requires sudo password)"
sudo cp /tmp/atgcflow-backend.service /etc/systemd/system/
sudo systemctl daemon-reload
echo "✅ Service installed"
echo ""

# Step 4: Enable and start service
echo "📝 Step 4: Enabling and starting service..."
sudo systemctl enable atgcflow-backend.service
sudo systemctl start atgcflow-backend.service
echo "✅ Service enabled and started"
echo ""

# Step 5: Check status
echo "📝 Step 5: Checking service status..."
sleep 2
sudo systemctl status atgcflow-backend.service --no-pager
echo ""

# Step 6: Start ngrok in background
echo "📝 Step 6: Starting ngrok tunnel..."

# Check if ngrok is installed
if ! command -v ngrok &> /dev/null; then
    echo "⚠️  ngrok is not installed!"
    echo "   Install ngrok from: https://ngrok.com/download"
    echo "   Or skip ngrok and use a different tunneling solution"
    echo ""
else
    # Kill any existing ngrok processes
    pkill -f ngrok || true

    # Start ngrok in background
    nohup ngrok http 8000 > "$BACKEND_DIR/logs/ngrok.log" 2>&1 &
    NGROK_PID=$!
    echo $NGROK_PID > "$BACKEND_DIR/ngrok.pid"

    echo "✅ ngrok started (PID: $NGROK_PID)"
    echo "   Waiting for ngrok to initialize..."
    sleep 3

    # Get ngrok URL
    NGROK_URL=$(curl -s http://localhost:4040/api/tunnels 2>/dev/null | grep -o '"public_url":"https://[^"]*' | grep -o 'https://[^"]*' | head -1)

    if [ -n "$NGROK_URL" ]; then
        echo "✅ ngrok tunnel active!"
        echo "   Public URL: $NGROK_URL"
        echo ""
        echo "⚠️  IMPORTANT: Update your backend/.env file:"
        echo "   Add this URL to CORS_ORIGINS: $NGROK_URL"
        echo "   Then restart backend: sudo systemctl restart atgcflow-backend.service"
        echo ""
    else
        echo "⚠️  Could not get ngrok URL"
        echo "   Check manually: http://localhost:4040"
        echo "   Or view logs: tail -f $BACKEND_DIR/logs/ngrok.log"
        echo ""
    fi
fi

echo "🎉 Setup Complete!"
echo ""
echo "📚 Useful commands:"
echo "   View status:  sudo systemctl status atgcflow-backend.service"
echo "   View logs:    sudo journalctl -u atgcflow-backend.service -f"
echo "   Restart:      sudo systemctl restart atgcflow-backend.service"
echo "   Stop:         sudo systemctl stop atgcflow-backend.service"
echo "   ngrok URL:    curl -s http://localhost:4040/api/tunnels | grep public_url"
echo ""
echo "📊 Log files:"
echo "   Backend:      tail -f $BACKEND_DIR/logs/backend.log"
echo "   Errors:       tail -f $BACKEND_DIR/logs/backend-error.log"
echo "   ngrok:        tail -f $BACKEND_DIR/logs/ngrok.log"
echo ""
echo "✅ Your backend is now running 24/7!"
echo "   Backend service: Running and auto-restarts on failure"
echo "   ngrok tunnel: Running (PID saved in ngrok.pid)"
echo "   Auto-start on boot: Enabled for backend service"
echo ""
echo "🔍 To view ngrok dashboard: http://localhost:4040"
echo ""
