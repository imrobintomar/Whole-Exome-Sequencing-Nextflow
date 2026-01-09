#!/bin/bash
# Complete ATGCFLOW Setup Script - Backend + ngrok

set -e

echo "🚀 ATGCFLOW Complete Service Setup"
echo "===================================="
echo "This will setup:"
echo "  1. Backend API (FastAPI + uvicorn)"
echo "  2. ngrok tunnel (optional)"
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

# Step 2: Setup Backend Service
echo "📝 Step 2: Setting up Backend Service..."
sed "s/%USER%/$CURRENT_USER/g" "$BACKEND_DIR/atgcflow-backend.service" > /tmp/atgcflow-backend.service
sudo cp /tmp/atgcflow-backend.service /etc/systemd/system/
sudo systemctl daemon-reload
sudo systemctl enable atgcflow-backend.service
sudo systemctl restart atgcflow-backend.service
echo "✅ Backend service installed and started"
echo ""

# Step 3: Check if ngrok is installed
echo "📝 Step 3: Checking for ngrok..."
if ! command -v ngrok &> /dev/null; then
    echo "⚠️  ngrok is not installed"
    echo ""
    echo "Do you want to skip ngrok setup? (y/n)"
    read -p "Skip ngrok? [y/N]: " SKIP_NGROK

    if [[ ! "$SKIP_NGROK" =~ ^[Yy]$ ]]; then
        echo ""
        echo "Please install ngrok first:"
        echo "  1. Visit: https://ngrok.com/download"
        echo "  2. Download and install ngrok"
        echo "  3. Run: ngrok authtoken YOUR_TOKEN"
        echo "  4. Run this script again"
        echo ""
        exit 1
    fi
    SETUP_NGROK=false
else
    echo "✅ ngrok is installed"
    SETUP_NGROK=true
fi
echo ""

# Step 4: Setup ngrok service (if ngrok is available)
if [ "$SETUP_NGROK" = true ]; then
    echo "📝 Step 4: Setting up ngrok service..."

    # Find ngrok path
    NGROK_PATH=$(which ngrok)

    # Update ngrok service file with correct path and user
    sed -e "s|%USER%|$CURRENT_USER|g" \
        -e "s|/usr/local/bin/ngrok|$NGROK_PATH|g" \
        "$BACKEND_DIR/atgcflow-ngrok.service" > /tmp/atgcflow-ngrok.service

    sudo cp /tmp/atgcflow-ngrok.service /etc/systemd/system/
    sudo systemctl daemon-reload
    sudo systemctl enable atgcflow-ngrok.service
    sudo systemctl restart atgcflow-ngrok.service
    echo "✅ ngrok service installed and started"
    echo ""

    # Wait for ngrok to start
    echo "   Waiting for ngrok to initialize..."
    sleep 5

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
        echo "⚠️  Could not get ngrok URL yet"
        echo "   Check status: sudo systemctl status atgcflow-ngrok.service"
        echo "   Or view dashboard: http://localhost:4040"
        echo ""
    fi
else
    echo "📝 Step 4: Skipping ngrok setup"
    echo ""
fi

# Step 5: Show service status
echo "📝 Step 5: Checking services status..."
echo ""
echo "Backend Service:"
sudo systemctl status atgcflow-backend.service --no-pager | head -10
echo ""

if [ "$SETUP_NGROK" = true ]; then
    echo "ngrok Service:"
    sudo systemctl status atgcflow-ngrok.service --no-pager | head -10
    echo ""
fi

echo "🎉 Setup Complete!"
echo ""
echo "📚 Useful commands:"
echo ""
echo "Backend:"
echo "  Status:   sudo systemctl status atgcflow-backend.service"
echo "  Logs:     tail -f $BACKEND_DIR/logs/backend.log"
echo "  Restart:  sudo systemctl restart atgcflow-backend.service"
echo "  Stop:     sudo systemctl stop atgcflow-backend.service"
echo ""

if [ "$SETUP_NGROK" = true ]; then
    echo "ngrok:"
    echo "  Status:   sudo systemctl status atgcflow-ngrok.service"
    echo "  Get URL:  ./manage-ngrok.sh url"
    echo "  Restart:  sudo systemctl restart atgcflow-ngrok.service"
    echo "  Stop:     sudo systemctl stop atgcflow-ngrok.service"
    echo "  Manage:   ./manage-ngrok.sh {start|stop|restart|status|url}"
    echo ""
fi

echo "📊 Quick Status Check:"
echo "  Run: ./check-status.sh"
echo ""
echo "🔍 ngrok Dashboard: http://localhost:4040"
echo ""
echo "✅ Your services are now running 24/7!"
echo "   Both services will automatically start on system boot"
echo "   Both services will automatically restart on failure"
echo ""
