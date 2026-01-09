#!/bin/bash
# Fix script for conda Python path issue

echo "🔧 Fixing ATGCFLOW Backend Service for Conda Python"
echo "===================================================="
echo ""

CURRENT_USER=$(whoami)
BACKEND_DIR="/media/drprabudh/m3/Nextflow-Script/WholeExome/backend"

echo "📝 Step 1: Updating service file with conda Python..."
sed "s/%USER%/$CURRENT_USER/g" "$BACKEND_DIR/atgcflow-backend.service" > /tmp/atgcflow-backend.service
echo "✅ Service file updated"
echo ""

echo "📝 Step 2: Reinstalling systemd service..."
echo "   (This requires sudo password)"
sudo cp /tmp/atgcflow-backend.service /etc/systemd/system/
sudo systemctl daemon-reload
echo "✅ Service reinstalled"
echo ""

echo "📝 Step 3: Restarting service..."
sudo systemctl restart atgcflow-backend.service
echo "✅ Service restarted"
echo ""

echo "📝 Step 4: Checking service status..."
sleep 2
sudo systemctl status atgcflow-backend.service --no-pager
echo ""

echo "🎉 Fix Applied!"
echo ""
echo "If you see 'active (running)' above, the fix worked!"
echo ""
