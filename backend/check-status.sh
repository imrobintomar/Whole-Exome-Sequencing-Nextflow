#!/bin/bash
# Quick status check script for ATGCFLOW Backend

echo "🔍 ATGCFLOW Backend Status Check"
echo "================================="
echo ""

# Check systemd service status
echo "📊 Service Status:"
if systemctl is-active --quiet atgcflow-backend.service; then
    echo "   ✅ Service is RUNNING"
    echo ""
    sudo systemctl status atgcflow-backend.service --no-pager | head -15
else
    echo "   ❌ Service is NOT running"
    echo ""
    sudo systemctl status atgcflow-backend.service --no-pager
fi
echo ""

# Check port
echo "🔌 Port 8000 Status:"
if lsof -i :8000 > /dev/null 2>&1; then
    echo "   ✅ Port 8000 is in use (backend is listening)"
    lsof -i :8000 | grep LISTEN
else
    echo "   ❌ Port 8000 is not in use"
fi
echo ""

# Test API endpoint
echo "🌐 API Health Check:"
if curl -s http://localhost:8000/docs > /dev/null; then
    echo "   ✅ Backend API is responding"
    echo "   📄 API Docs: http://localhost:8000/docs"
else
    echo "   ❌ Backend API is not responding"
fi
echo ""

# Show recent logs
echo "📝 Recent Logs (last 10 lines):"
echo "   ------------------------------------------------------------"
tail -10 logs/backend.log 2>/dev/null || echo "   No logs found"
echo "   ------------------------------------------------------------"
echo ""

# Check ngrok status
echo "🌍 ngrok Tunnel Status:"
if systemctl is-active --quiet atgcflow-ngrok.service 2>/dev/null; then
    echo "   ✅ ngrok service is RUNNING"
    NGROK_URL=$(curl -s http://localhost:4040/api/tunnels 2>/dev/null | grep -o '"public_url":"https://[^"]*' | grep -o 'https://[^"]*' | head -1)
    if [ -n "$NGROK_URL" ]; then
        echo "   🔗 Public URL: $NGROK_URL"
        echo "   📊 Dashboard: http://localhost:4040"
    else
        echo "   ⚠️  Could not get ngrok URL"
        echo "   📊 Check dashboard: http://localhost:4040"
    fi
elif pgrep -f "ngrok http" > /dev/null 2>&1; then
    echo "   ✅ ngrok is RUNNING (not as service)"
    NGROK_URL=$(curl -s http://localhost:4040/api/tunnels 2>/dev/null | grep -o '"public_url":"https://[^"]*' | grep -o 'https://[^"]*' | head -1)
    if [ -n "$NGROK_URL" ]; then
        echo "   🔗 Public URL: $NGROK_URL"
        echo "   📊 Dashboard: http://localhost:4040"
    else
        echo "   ⚠️  Could not get ngrok URL"
    fi
else
    echo "   ❌ ngrok is NOT running"
    echo "   💡 Start ngrok: ./manage-ngrok.sh start"
    echo "   💡 Or setup as service: ./setup-all-services.sh"
fi
echo ""

echo "💡 Useful Commands:"
echo "   Backend:"
echo "     View live logs:  tail -f logs/backend.log"
echo "     View errors:     tail -f logs/backend-error.log"
echo "     Restart:         sudo systemctl restart atgcflow-backend.service"
echo "     Stop:            sudo systemctl stop atgcflow-backend.service"
echo ""
echo "   ngrok:"
echo "     Get URL:         ./manage-ngrok.sh url"
echo "     Manage:          ./manage-ngrok.sh {start|stop|restart|status}"
echo "     Restart:         sudo systemctl restart atgcflow-ngrok.service"
echo ""
