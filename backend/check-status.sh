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

echo "💡 Useful Commands:"
echo "   View live logs:  tail -f logs/backend.log"
echo "   View errors:     tail -f logs/backend-error.log"
echo "   Restart:         sudo systemctl restart atgcflow-backend.service"
echo "   Stop:            sudo systemctl stop atgcflow-backend.service"
echo ""
