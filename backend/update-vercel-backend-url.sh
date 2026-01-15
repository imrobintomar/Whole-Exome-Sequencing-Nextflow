#!/bin/bash
# Update Vercel environment variable with production API URL

set -e

BACKEND_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$BACKEND_DIR")"
PRODUCTION_API="https://api.atgcflow.com"

echo "🔄 Updating Vercel with Production API URL"
echo "==========================================="
echo ""
echo "✅ Production API URL: $PRODUCTION_API"
echo ""

# Update local .env.production file
ENV_PROD_FILE="$PROJECT_DIR/frontend/.env.production"
if [ -f "$ENV_PROD_FILE" ]; then
    sed -i "s|^NEXT_PUBLIC_API_URL=.*|NEXT_PUBLIC_API_URL=$PRODUCTION_API|g" "$ENV_PROD_FILE"
    echo "✅ Updated $ENV_PROD_FILE"
else
    echo "⚠️  Frontend .env.production not found"
fi

# Check if Vercel CLI is installed
if ! command -v vercel &> /dev/null; then
    echo ""
    echo "⚠️  Vercel CLI not installed!"
    echo ""
    echo "   Install with: npm i -g vercel"
    echo ""
    echo "   After installation, run this script again to auto-update Vercel"
    echo "   Or manually update in Vercel Dashboard:"
    echo "   https://vercel.com/your-project/settings/environment-variables"
    echo ""
    echo "   Set: NEXT_PUBLIC_API_URL = $NGROK_URL"
    echo ""
    exit 1
fi

# Update Vercel environment variable
echo "📝 Updating Vercel environment variable..."
echo ""

cd "$PROJECT_DIR/frontend" || exit 1

# Remove old variable (if exists)
vercel env rm NEXT_PUBLIC_API_URL production --yes 2>/dev/null || true

# Add new variable
echo "$PRODUCTION_API" | vercel env add NEXT_PUBLIC_API_URL production

echo ""
echo "✅ Vercel environment variable updated"
echo ""

# Ask to redeploy
read -p "Do you want to trigger a Vercel deployment now? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "🚀 Deploying to Vercel..."
    vercel --prod
    echo ""
    echo "✅ Deployment triggered!"
else
    echo ""
    echo "⚠️  Remember to redeploy manually:"
    echo "   cd $PROJECT_DIR/frontend && vercel --prod"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "🎉 Update Complete!"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "📊 Current Configuration:"
echo "   Production API:   $PRODUCTION_API"
echo "   Frontend (prod):  https://atgcflow.com/"
echo ""
echo "💡 Tips:"
echo "   - Production API is served via Cloudflare Tunnel"
echo "   - API URL is stable and doesn't change on restart"
echo "   - Cloudflare provides automatic HTTPS and DDoS protection"
echo ""
