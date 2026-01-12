#!/bin/bash
# Test authentication flow

echo "🔍 Testing Authentication Flow"
echo "=============================="
echo ""

NGROK_URL="https://disbursable-ennoblingly-jerlene.ngrok-free.dev"

echo "1️⃣  Testing root endpoint (should work without auth)"
echo "-------------------------------------------------------"
RESPONSE=$(curl -s -w "\nHTTP_CODE:%{http_code}" -H "ngrok-skip-browser-warning: true" "$NGROK_URL/")
HTTP_CODE=$(echo "$RESPONSE" | grep "HTTP_CODE" | cut -d: -f2)
BODY=$(echo "$RESPONSE" | grep -v "HTTP_CODE")

echo "Response: $BODY"
echo "HTTP Code: $HTTP_CODE"
echo ""

if [ "$HTTP_CODE" = "200" ]; then
    echo "✅ Root endpoint working"
else
    echo "❌ Root endpoint failed"
fi
echo ""

echo "2️⃣  Testing /jobs endpoint WITHOUT auth (should return 401 or 403)"
echo "-------------------------------------------------------------------"
RESPONSE=$(curl -s -w "\nHTTP_CODE:%{http_code}" -H "ngrok-skip-browser-warning: true" "$NGROK_URL/jobs")
HTTP_CODE=$(echo "$RESPONSE" | grep "HTTP_CODE" | cut -d: -f2)
BODY=$(echo "$RESPONSE" | grep -v "HTTP_CODE")

echo "Response: $BODY"
echo "HTTP Code: $HTTP_CODE"
echo ""

if [ "$HTTP_CODE" = "401" ] || [ "$HTTP_CODE" = "403" ]; then
    echo "✅ /jobs correctly requires authentication"
else
    echo "❌ /jobs should require authentication but returned: $HTTP_CODE"
fi
echo ""

echo "3️⃣  Testing /jobs endpoint WITH invalid token"
echo "----------------------------------------------"
RESPONSE=$(curl -s -w "\nHTTP_CODE:%{http_code}" \
    -H "ngrok-skip-browser-warning: true" \
    -H "Authorization: Bearer invalid_token_12345" \
    "$NGROK_URL/jobs")
HTTP_CODE=$(echo "$RESPONSE" | grep "HTTP_CODE" | cut -d: -f2)
BODY=$(echo "$RESPONSE" | grep -v "HTTP_CODE")

echo "Response: $BODY"
echo "HTTP Code: $HTTP_CODE"
echo ""

if [ "$HTTP_CODE" = "401" ]; then
    echo "✅ Invalid token correctly rejected"
else
    echo "⚠️  Expected 401, got: $HTTP_CODE"
fi
echo ""

echo "4️⃣  Check CORS headers"
echo "----------------------"
curl -v -H "Origin: https://atgcflow.com" \
    -H "ngrok-skip-browser-warning: true" \
    "$NGROK_URL/" 2>&1 | grep -i "access-control\|cors" | head -5
echo ""

echo "5️⃣  Testing OPTIONS preflight request"
echo "--------------------------------------"
curl -X OPTIONS -v \
    -H "Origin: https://atgcflow.com" \
    -H "Access-Control-Request-Method: GET" \
    -H "Access-Control-Request-Headers: Authorization" \
    "$NGROK_URL/jobs" 2>&1 | grep -i "HTTP\|access-control" | head -10
echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Summary"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "If you're getting 401 errors after logging in, check:"
echo "  1. Browser Console → Network tab → /jobs request → Headers"
echo "  2. Look for 'Authorization: Bearer <token>'"
echo "  3. If missing, frontend auth interceptor isn't working"
echo "  4. If present but still 401, check backend logs:"
echo "     tail -f /media/drprabudh/m3/Nextflow-Script/WholeExome/backend/logs/backend-error.log"
echo ""
