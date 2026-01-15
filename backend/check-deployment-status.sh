#!/bin/bash
# Quick deployment status checker for ATGCFLOW

set -e

BACKEND_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$BACKEND_DIR")"
PRODUCTION_API="https://api.atgcflow.com"

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

# 2. Check Production API (Cloudflare Tunnel)
echo "2️⃣  Production API Status (Cloudflare):"
if curl -s -o /dev/null -w "%{http_code}" "$PRODUCTION_API/" 2>/dev/null | grep -q "200\|307\|301"; then
    echo "   ✅ Production API accessible: $PRODUCTION_API"
else
    echo "   ⚠️  Production API not responding: $PRODUCTION_API"
    echo "      Check Cloudflare Tunnel configuration"
fi
echo ""

# 3. Check Frontend Configuration
echo "3️⃣  Frontend Configuration:"
ENV_PROD_FILE="$PROJECT_DIR/frontend/.env.production"
if [ -f "$ENV_PROD_FILE" ]; then
    FRONTEND_API_URL=$(grep "^NEXT_PUBLIC_API_URL=" "$ENV_PROD_FILE" | cut -d'=' -f2-)
    echo "   Frontend expects: $FRONTEND_API_URL"

    if [ "$FRONTEND_API_URL" = "$PRODUCTION_API" ]; then
        echo "   ✅ Frontend .env.production configured correctly"
    else
        echo "   ⚠️  Frontend .env.production should be: $PRODUCTION_API"
        echo "      Got: $FRONTEND_API_URL"
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

    if echo "$CORS_ORIGINS" | grep -q "https://atgcflow.com\|https://api.atgcflow.com"; then
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
echo "   ℹ️  Production frontend: https://atgcflow.com"
echo "   ℹ️  Production API: $PRODUCTION_API"
echo ""
echo "   Check Vercel environment variables:"
echo "   https://vercel.com/your-project/settings/environment-variables"
echo ""
echo "   NEXT_PUBLIC_API_URL should be: $PRODUCTION_API"
echo ""

# Summary
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "📊 Summary"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "🌐 Production Configuration:"
echo "   API:      $PRODUCTION_API"
echo "   Frontend: https://atgcflow.com"
echo ""

echo "📝 Next Steps:"
if ! systemctl is-active --quiet atgcflow-backend.service; then
    echo "   1. Start backend: sudo systemctl start atgcflow-backend.service"
fi

echo ""
echo "💡 Quick Commands:"
echo "   Status check:     $BACKEND_DIR/check-deployment-status.sh"
echo "   View logs:        sudo journalctl -u atgcflow-backend.service -f"
echo "   Restart backend:  sudo systemctl restart atgcflow-backend.service"
echo ""
