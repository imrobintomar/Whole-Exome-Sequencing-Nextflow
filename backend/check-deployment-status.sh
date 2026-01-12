#!/bin/bash
# Quick deployment status checker for ATGCFLOW

set -e

BACKEND_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$BACKEND_DIR")"

echo "🔍 ATGCFLOW Deployment Status Check"
echo "===================================="
echo ""

# 1. Check Backend Service
echo "1️⃣  Backend Service Status:"
if systemctl is-active --quiet atgcflow-backend.service; then
    echo "   ✅ Backend service is running"

    # Test local endpoint
    if curl -s http://localhost:8000/ > /dev/null 2>&1; then
        echo "   ✅ Backend responds on http://localhost:8000"
    else
        echo "   ❌ Backend not responding on http://localhost:8000"
    fi
else
    echo "   ❌ Backend service is NOT running"
    echo "      Start with: sudo systemctl start atgcflow-backend.service"
fi
echo ""

# 2. Check ngrok Service
echo "2️⃣  ngrok Tunnel Status:"
if systemctl is-active --quiet atgcflow-ngrok.service; then
    echo "   ✅ ngrok service is running"

    # Get ngrok URL
    NGROK_URL=$(curl -s http://localhost:4040/api/tunnels 2>/dev/null | grep -o '"public_url":"https://[^"]*' | grep -o 'https://[^"]*' | head -1)

    if [ -n "$NGROK_URL" ]; then
        echo "   ✅ ngrok tunnel active: $NGROK_URL"

        # Test ngrok endpoint
        if curl -s -o /dev/null -w "%{http_code}" "$NGROK_URL/" 2>/dev/null | grep -q "200\|307\|301"; then
            echo "   ✅ Backend accessible via ngrok"
        else
            echo "   ⚠️  ngrok tunnel exists but backend not accessible"
        fi
    else
        echo "   ❌ Could not retrieve ngrok URL"
    fi
else
    echo "   ❌ ngrok service is NOT running"
    echo "      Start with: sudo systemctl start atgcflow-ngrok.service"
    NGROK_URL=""
fi
echo ""

# 3. Check Frontend Configuration
echo "3️⃣  Frontend Configuration:"
ENV_PROD_FILE="$PROJECT_DIR/frontend/.env.production"
if [ -f "$ENV_PROD_FILE" ]; then
    FRONTEND_API_URL=$(grep "^NEXT_PUBLIC_API_URL=" "$ENV_PROD_FILE" | cut -d'=' -f2-)
    echo "   Frontend expects: $FRONTEND_API_URL"

    if [ -n "$NGROK_URL" ]; then
        if [ "$FRONTEND_API_URL" = "$NGROK_URL" ]; then
            echo "   ✅ Frontend .env.production matches current ngrok URL"
        else
            echo "   ⚠️  Frontend .env.production does NOT match current ngrok URL"
            echo "      Expected: $NGROK_URL"
            echo "      Got:      $FRONTEND_API_URL"
        fi
    fi
else
    echo "   ❌ Frontend .env.production not found"
fi
echo ""

# 4. Check CORS Configuration
echo "4️⃣  CORS Configuration:"
ENV_FILE="$BACKEND_DIR/.env"
if [ -f "$ENV_FILE" ]; then
    CORS_ORIGINS=$(grep "^CORS_ORIGINS=" "$ENV_FILE" | cut -d'=' -f2-)

    if [ -n "$NGROK_URL" ]; then
        if echo "$CORS_ORIGINS" | grep -q "$NGROK_URL"; then
            echo "   ✅ Current ngrok URL is in CORS_ORIGINS"
        else
            echo "   ⚠️  Current ngrok URL NOT in CORS_ORIGINS"
            echo "      Add it with: Restart backend service or run setup script"
        fi
    fi

    if echo "$CORS_ORIGINS" | grep -q "https://atgcflow.com"; then
        echo "   ✅ Production domain (atgcflow.com) in CORS_ORIGINS"
    else
        echo "   ⚠️  Production domain (atgcflow.com) NOT in CORS_ORIGINS"
    fi
else
    echo "   ❌ Backend .env not found"
fi
echo ""

# 5. Vercel Deployment Check
echo "5️⃣  Vercel Deployment:"
if [ -n "$NGROK_URL" ] && [ -n "$FRONTEND_API_URL" ]; then
    if [ "$FRONTEND_API_URL" = "$NGROK_URL" ]; then
        echo "   ℹ️  Local config is correct"
        echo "   ⚠️  Check Vercel environment variables:"
        echo "      https://vercel.com/your-project/settings/environment-variables"
        echo ""
        echo "      NEXT_PUBLIC_API_URL should be: $NGROK_URL"
    else
        echo "   ⚠️  Mismatch detected!"
        echo "      Run: $BACKEND_DIR/update-vercel-backend-url.sh"
    fi
else
    echo "   ⚠️  Cannot verify - ngrok URL or frontend config missing"
fi
echo ""

# Summary
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "📊 Summary"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

if [ -n "$NGROK_URL" ]; then
    echo "🌐 Current ngrok URL:"
    echo "   $NGROK_URL"
    echo ""
fi

echo "📝 Next Steps:"
if ! systemctl is-active --quiet atgcflow-backend.service; then
    echo "   1. Start backend: sudo systemctl start atgcflow-backend.service"
fi

if ! systemctl is-active --quiet atgcflow-ngrok.service; then
    echo "   2. Start ngrok: sudo systemctl start atgcflow-ngrok.service"
fi

if [ -n "$NGROK_URL" ] && [ "$FRONTEND_API_URL" != "$NGROK_URL" ]; then
    echo "   3. Update frontend config: $BACKEND_DIR/update-vercel-backend-url.sh"
fi

echo "   4. Update Vercel environment variable:"
echo "      Variable: NEXT_PUBLIC_API_URL"
echo "      Value:    $NGROK_URL"
echo "   5. Redeploy frontend on Vercel"
echo ""

echo "💡 Quick Commands:"
echo "   Status check:     $BACKEND_DIR/check-deployment-status.sh"
echo "   Update Vercel:    $BACKEND_DIR/update-vercel-backend-url.sh"
echo "   View logs:        sudo journalctl -u atgcflow-backend.service -f"
echo "   ngrok dashboard:  http://localhost:4040"
echo ""
