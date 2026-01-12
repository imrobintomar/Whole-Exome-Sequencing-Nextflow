#!/bin/bash
# Improved ATGCFLOW Backend Service Setup with ngrok Integration

set -e

echo "🚀 ATGCFLOW Backend Service Setup (with ngrok)"
echo "================================================"
echo ""

# Configuration
CURRENT_USER=$(whoami)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BACKEND_DIR="$SCRIPT_DIR"
PROJECT_DIR="$(dirname "$BACKEND_DIR")"
NGROK_CONFIG="$PROJECT_DIR/ngrok.yml"

echo "📋 Configuration:"
echo "   User:              $CURRENT_USER"
echo "   Backend directory: $BACKEND_DIR"
echo "   Project directory: $PROJECT_DIR"
echo "   ngrok config:      $NGROK_CONFIG"
echo ""

# Step 1: Validate prerequisites
echo "🔍 Step 1: Validating prerequisites..."

# Check if ngrok is installed
if ! command -v ngrok &> /dev/null; then
    echo "❌ ngrok is not installed!"
    echo ""
    echo "   Install ngrok:"
    echo "   1. Visit: https://ngrok.com/download"
    echo "   2. Or run: curl -s https://ngrok-agent.s3.amazonaws.com/ngrok.asc | sudo tee /etc/apt/trusted.gpg.d/ngrok.asc >/dev/null"
    echo "              echo \"deb https://ngrok-agent.s3.amazonaws.com buster main\" | sudo tee /etc/apt/sources.list.d/ngrok.list"
    echo "              sudo apt update && sudo apt install ngrok"
    echo ""
    exit 1
fi
echo "   ✅ ngrok installed: $(which ngrok)"

# Check if ngrok config exists
if [ ! -f "$NGROK_CONFIG" ]; then
    echo "⚠️  ngrok config not found: $NGROK_CONFIG"
    echo "   Creating default config..."
    cat > "$NGROK_CONFIG" << 'EOF'
version: "2"
authtoken: YOUR_NGROK_AUTHTOKEN_HERE

tunnels:
  wes-backend:
    proto: http
    addr: 8000
    schemes:
      - https
    inspect: true
    response_headers:
      add:
        - "ngrok-skip-browser-warning: true"
EOF
    echo "   ✅ Created $NGROK_CONFIG"
    echo ""
    echo "⚠️  IMPORTANT: Update ngrok.yml with your authtoken!"
    echo "   Run: ngrok config add-authtoken YOUR_AUTH_TOKEN"
    echo "   Get token from: https://dashboard.ngrok.com/get-started/your-authtoken"
    echo ""
    exit 1
fi
echo "   ✅ ngrok config found"

# Validate ngrok authtoken in config
if grep -q "YOUR_NGROK_AUTHTOKEN_HERE" "$NGROK_CONFIG"; then
    echo "❌ ngrok authtoken not configured!"
    echo ""
    echo "   Update $NGROK_CONFIG with your authtoken"
    echo "   Or run: ngrok config add-authtoken YOUR_AUTH_TOKEN"
    echo "   Get token from: https://dashboard.ngrok.com/get-started/your-authtoken"
    echo ""
    exit 1
fi
echo "   ✅ ngrok authenticated"

echo "✅ Prerequisites validated"
echo ""

# Step 2: Create logs directory
echo "📝 Step 2: Creating logs directory..."
mkdir -p "$BACKEND_DIR/logs"
echo "✅ Logs directory created"
echo ""

# Step 3: Update service file with current user
echo "📝 Step 3: Preparing systemd service files..."

# Create backend service
sed "s|%USER%|$CURRENT_USER|g" "$BACKEND_DIR/atgcflow-backend.service" > /tmp/atgcflow-backend.service

# Create ngrok service
cat > /tmp/atgcflow-ngrok.service << EOF
[Unit]
Description=ngrok tunnel for ATGCFLOW Backend
After=network.target atgcflow-backend.service

[Service]
Type=simple
User=$CURRENT_USER
WorkingDirectory=$PROJECT_DIR
ExecStart=/usr/bin/env ngrok start --config=$NGROK_CONFIG wes-backend
Restart=always
RestartSec=10
StandardOutput=append:$BACKEND_DIR/logs/ngrok.log
StandardError=append:$BACKEND_DIR/logs/ngrok-error.log

[Install]
WantedBy=multi-user.target
EOF

echo "✅ Service files prepared"
echo ""

# Step 4: Install systemd services
echo "📝 Step 4: Installing systemd services..."
echo "   (This requires sudo password)"

sudo cp /tmp/atgcflow-backend.service /etc/systemd/system/
sudo cp /tmp/atgcflow-ngrok.service /etc/systemd/system/
sudo systemctl daemon-reload

echo "✅ Services installed"
echo ""

# Step 5: Enable and start services
echo "📝 Step 5: Enabling and starting services..."

sudo systemctl enable atgcflow-backend.service
sudo systemctl enable atgcflow-ngrok.service

sudo systemctl start atgcflow-backend.service
echo "   ✅ Backend service started"

# Wait for backend to be ready before starting ngrok
echo "   ⏳ Waiting for backend to be ready..."
MAX_ATTEMPTS=30
ATTEMPT=0
while [ $ATTEMPT -lt $MAX_ATTEMPTS ]; do
    if curl -s http://localhost:8000/ > /dev/null 2>&1; then
        echo "   ✅ Backend is ready"
        break
    fi
    ATTEMPT=$((ATTEMPT + 1))
    if [ $ATTEMPT -eq $MAX_ATTEMPTS ]; then
        echo "   ❌ Backend failed to start"
        echo "      Check logs: sudo journalctl -u atgcflow-backend.service -n 50"
        exit 1
    fi
    sleep 1
done

sudo systemctl start atgcflow-ngrok.service
echo "   ✅ ngrok service started"

echo "✅ Services enabled and started"
echo ""

# Step 6: Wait for ngrok and get URL
echo "📝 Step 6: Getting ngrok public URL..."
echo "   ⏳ Waiting for ngrok tunnel to establish..."

sleep 3
MAX_ATTEMPTS=20
ATTEMPT=0
NGROK_URL=""

while [ $ATTEMPT -lt $MAX_ATTEMPTS ]; do
    NGROK_URL=$(curl -s http://localhost:4040/api/tunnels 2>/dev/null | grep -o '"public_url":"https://[^"]*' | grep -o 'https://[^"]*' | head -1)

    if [ -n "$NGROK_URL" ]; then
        break
    fi

    ATTEMPT=$((ATTEMPT + 1))
    if [ $ATTEMPT -eq $MAX_ATTEMPTS ]; then
        echo "⚠️  Could not retrieve ngrok URL automatically"
        echo "   Check manually: http://localhost:4040"
        echo "   Or view logs: sudo journalctl -u atgcflow-ngrok.service -n 50"
        echo ""
        break
    fi
    sleep 1
done

if [ -n "$NGROK_URL" ]; then
    echo "✅ ngrok tunnel established!"
    echo ""
    echo "   Public URL: $NGROK_URL"
    echo ""

    # Save URL to file for easy access
    echo "$NGROK_URL" > "$BACKEND_DIR/ngrok-url.txt"

    # Try to automatically update CORS in .env
    ENV_FILE="$BACKEND_DIR/.env"
    if [ -f "$ENV_FILE" ]; then
        echo "📝 Updating CORS configuration..."

        # Check if CORS_ORIGINS exists in .env
        if grep -q "^CORS_ORIGINS=" "$ENV_FILE"; then
            # Get current CORS origins
            CURRENT_CORS=$(grep "^CORS_ORIGINS=" "$ENV_FILE" | cut -d'=' -f2-)

            # Check if ngrok URL already in CORS
            if ! echo "$CURRENT_CORS" | grep -q "$NGROK_URL"; then
                # Add ngrok URL to CORS_ORIGINS
                NEW_CORS="${CURRENT_CORS},${NGROK_URL}"
                sed -i "s|^CORS_ORIGINS=.*|CORS_ORIGINS=${NEW_CORS}|g" "$ENV_FILE"

                echo "   ✅ Added ngrok URL to CORS_ORIGINS"
                echo "   ⏳ Restarting backend to apply changes..."
                sudo systemctl restart atgcflow-backend.service
                sleep 2
                echo "   ✅ Backend restarted"
            else
                echo "   ℹ️  ngrok URL already in CORS_ORIGINS"
            fi
        else
            echo "   ⚠️  CORS_ORIGINS not found in .env"
            echo "      Please add manually: CORS_ORIGINS=http://localhost:3000,$NGROK_URL"
        fi
    else
        echo "⚠️  .env file not found: $ENV_FILE"
        echo "   Please create .env with: CORS_ORIGINS=$NGROK_URL"
    fi
    echo ""

    # Update frontend .env.production
    ENV_PROD_FILE="$PROJECT_DIR/frontend/.env.production"
    if [ -f "$ENV_PROD_FILE" ]; then
        echo "📝 Updating frontend production config..."
        sed -i "s|^NEXT_PUBLIC_API_URL=.*|NEXT_PUBLIC_API_URL=$NGROK_URL|g" "$ENV_PROD_FILE"
        echo "   ✅ Updated frontend/.env.production"
        echo "   ⚠️  Remember to update Vercel and redeploy!"
        echo "      Run: $BACKEND_DIR/update-vercel-backend-url.sh"
    fi
    echo ""
fi

# Step 7: Check status
echo "📝 Step 7: Service status check..."
echo ""
echo "Backend Service:"
sudo systemctl status atgcflow-backend.service --no-pager | head -15
echo ""
echo "ngrok Service:"
sudo systemctl status atgcflow-ngrok.service --no-pager | head -15
echo ""

# Cleanup temp files
rm -f /tmp/atgcflow-backend.service /tmp/atgcflow-ngrok.service

# Final summary
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "🎉 Setup Complete!"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "📊 Service Status:"
echo "   ✅ Backend service: Running on http://localhost:8000"
if [ -n "$NGROK_URL" ]; then
    echo "   ✅ ngrok tunnel:    $NGROK_URL"
else
    echo "   ⚠️  ngrok tunnel:    Check http://localhost:4040"
fi
echo "   ✅ Auto-start:      Enabled on boot"
echo ""
echo "📚 Useful commands:"
echo "   View backend status:  sudo systemctl status atgcflow-backend.service"
echo "   View ngrok status:    sudo systemctl status atgcflow-ngrok.service"
echo "   View backend logs:    sudo journalctl -u atgcflow-backend.service -f"
echo "   View ngrok logs:      sudo journalctl -u atgcflow-ngrok.service -f"
echo "   Restart backend:      sudo systemctl restart atgcflow-backend.service"
echo "   Restart ngrok:        sudo systemctl restart atgcflow-ngrok.service"
echo "   Stop all:             sudo systemctl stop atgcflow-backend.service atgcflow-ngrok.service"
echo "   Get ngrok URL:        curl -s http://localhost:4040/api/tunnels | jq -r '.tunnels[0].public_url'"
echo "   Or read from file:    cat $BACKEND_DIR/ngrok-url.txt"
echo "   Update Vercel:        $BACKEND_DIR/update-vercel-backend-url.sh"
echo ""
echo "🔍 Web interfaces:"
echo "   Backend API:      http://localhost:8000"
if [ -n "$NGROK_URL" ]; then
    echo "   Public API:       $NGROK_URL"
fi
echo "   ngrok Dashboard:  http://localhost:4040"
echo ""
echo "💡 Tips:"
echo "   - ngrok URL is saved in: $BACKEND_DIR/ngrok-url.txt"
echo "   - Free ngrok URLs change on restart (consider ngrok Pro for static domains)"
echo "   - Both services auto-restart on failure and boot"
echo "   - Check service health: curl http://localhost:8000/"
echo ""
echo "⚠️  IMPORTANT for https://atgcflow.com/:"
echo "   Your production site needs the ngrok URL updated in Vercel!"
echo "   Run: $BACKEND_DIR/update-vercel-backend-url.sh"
echo "   Or manually update in Vercel Dashboard and redeploy"
echo ""
echo "✅ Your backend is now running 24/7 with public ngrok access!"
echo ""
