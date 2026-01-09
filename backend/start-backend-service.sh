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

echo "🎉 Setup Complete!"
echo ""
echo "📚 Useful commands:"
echo "   View status:  sudo systemctl status atgcflow-backend.service"
echo "   View logs:    sudo journalctl -u atgcflow-backend.service -f"
echo "   Restart:      sudo systemctl restart atgcflow-backend.service"
echo "   Stop:         sudo systemctl stop atgcflow-backend.service"
echo ""
echo "📊 Log files:"
echo "   Backend:      tail -f $BACKEND_DIR/logs/backend.log"
echo "   Errors:       tail -f $BACKEND_DIR/logs/backend-error.log"
echo ""
echo "✅ Your backend is now running 24/7!"
echo "   It will automatically start on system boot."
echo "   It will automatically restart if it crashes."
echo ""
