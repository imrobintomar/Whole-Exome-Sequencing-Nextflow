#!/bin/bash
# Auto-update Vercel environment variable with current ngrok URL

set -e

BACKEND_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$BACKEND_DIR")"

echo "🔄 Updating Vercel with ngrok URL"
echo "=================================="
echo ""

# Get current ngrok URL
NGROK_URL=$(curl -s http://localhost:4040/api/tunnels 2>/dev/null | grep -o '"public_url":"https://[^"]*' | grep -o 'https://[^"]*' | head -1)

if [ -z "$NGROK_URL" ]; then
    echo "❌ Could not retrieve ngrok URL"
    echo "   Is ngrok running? Check: sudo systemctl status atgcflow-ngrok.service"
    exit 1
fi

echo "✅ Current ngrok URL: $NGROK_URL"
echo ""

# Update local .env.production file
ENV_PROD_FILE="$PROJECT_DIR/frontend/.env.production"
if [ -f "$ENV_PROD_FILE" ]; then
    sed -i "s|^NEXT_PUBLIC_API_URL=.*|NEXT_PUBLIC_API_URL=$NGROK_URL|g" "$ENV_PROD_FILE"
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
echo "$NGROK_URL" | vercel env add NEXT_PUBLIC_API_URL production

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
echo "   ngrok URL:        $NGROK_URL"
echo "   Frontend (prod):  https://atgcflow.com/"
echo ""
echo "💡 Tips:"
echo "   - Free ngrok URLs change on restart"
echo "   - Run this script after ngrok restarts"
echo "   - Consider ngrok Pro for static domains"
echo ""
