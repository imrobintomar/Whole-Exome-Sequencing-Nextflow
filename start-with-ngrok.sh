#!/bin/bash

# Start Backend with ngrok Tunnel
# This script starts the backend server and creates a public ngrok tunnel

echo "========================================"
echo "  WES Pipeline - Backend + ngrok Setup"
echo "========================================"
echo ""

# Check if ngrok is installed
if ! command -v ngrok &> /dev/null; then
    echo "❌ ngrok is not installed!"
    echo ""
    echo "Please install ngrok:"
    echo "  Visit: https://ngrok.com/download"
    echo "  Or run: sudo apt install ngrok (on Ubuntu/Debian)"
    echo ""
    exit 1
fi

# Check if backend virtual environment exists
if [ ! -d "backend/venv" ]; then
    echo "❌ Backend virtual environment not found!"
    echo "Please run: cd backend && python -m venv venv && source venv/bin/activate && pip install -r requirements.txt"
    exit 1
fi

echo "0️⃣  Cleaning up any existing processes..."
echo ""

# Kill any existing backend processes
lsof -ti:8000 | xargs kill -9 2>/dev/null || true
pkill -f "python main.py" 2>/dev/null || true
pkill -f ngrok 2>/dev/null || true
sleep 2

echo "✅ Cleanup complete"
echo ""
echo "1️⃣  Starting Backend Server..."
echo ""

# Start backend in background
cd backend
source venv/bin/activate
python main.py > backend.log 2>&1 &
BACKEND_PID=$!
cd ..

echo "✅ Backend started (PID: $BACKEND_PID)"
echo "   Logs: backend/backend.log"
echo ""

# Wait for backend to be ready
echo "⏳ Waiting for backend to start..."
MAX_ATTEMPTS=15
ATTEMPT=0
while [ $ATTEMPT -lt $MAX_ATTEMPTS ]; do
    if curl -s http://localhost:8000/ > /dev/null 2>&1; then
        echo "✅ Backend is running at http://localhost:8000"
        break
    fi
    ATTEMPT=$((ATTEMPT + 1))
    if [ $ATTEMPT -eq $MAX_ATTEMPTS ]; then
        echo "❌ Backend failed to start after 15 seconds"
        echo "   Check logs: tail -f backend/backend.log"
        kill $BACKEND_PID 2>/dev/null
        exit 1
    fi
    sleep 1
done

echo ""
echo "2️⃣  Starting ngrok tunnel..."
echo ""

# Start ngrok in background
ngrok http 8000 --log=stdout > ngrok.log 2>&1 &
NGROK_PID=$!

echo "✅ ngrok started (PID: $NGROK_PID)"
echo ""

# Wait for ngrok to be ready
echo "⏳ Waiting for ngrok tunnel..."
MAX_ATTEMPTS=10
ATTEMPT=0
NGROK_URL=""

while [ $ATTEMPT -lt $MAX_ATTEMPTS ]; do
    NGROK_URL=$(curl -s http://localhost:4040/api/tunnels 2>/dev/null | grep -o '"public_url":"https://[^"]*' | grep -o 'https://[^"]*' | head -n 1)
    
    if [ ! -z "$NGROK_URL" ]; then
        break
    fi
    
    ATTEMPT=$((ATTEMPT + 1))
    if [ $ATTEMPT -eq $MAX_ATTEMPTS ]; then
        echo "❌ Failed to get ngrok URL. Check if ngrok is authenticated."
        echo "   Run: ngrok config add-authtoken YOUR_AUTH_TOKEN"
        echo "   Get token from: https://dashboard.ngrok.com/get-started/your-authtoken"
        echo ""
        echo "   Or check logs: cat ngrok.log"
        kill $BACKEND_PID 2>/dev/null
        kill $NGROK_PID 2>/dev/null
        exit 1
    fi
    sleep 1
done

echo ""
echo "========================================"
echo "  🎉 Setup Complete!"
echo "========================================"
echo ""
echo "📍 Backend (Local):  http://localhost:8000"
echo "🌐 Backend (Public): $NGROK_URL"
echo ""
echo "========================================"
echo "  Next Steps:"
echo "========================================"
echo ""
echo "1. Copy your ngrok URL: $NGROK_URL"
echo ""
echo "2. Set Vercel Environment Variable:"
echo "   NEXT_PUBLIC_API_URL=$NGROK_URL"
echo ""
echo "3. After getting your Vercel URL, update backend CORS:"
echo "   Edit backend/.env and add your Vercel URL to CORS_ORIGINS"
echo ""
echo "4. Monitor services:"
echo "   - Backend logs: tail -f backend/backend.log"
echo "   - ngrok dashboard: http://localhost:4040"
echo "   - Backend health: curl http://localhost:8000/"
echo ""
echo "5. Stop services:"
echo "   ./stop-services.sh"
echo "   OR"
echo "   kill $BACKEND_PID $NGROK_PID"
echo ""
echo "⚠️  NOTE: Free ngrok URL changes on restart!"
echo "   Consider upgrading to ngrok Pro for static domain."
echo ""

# Save PIDs to file for easy cleanup
echo "$BACKEND_PID" > .backend.pid
echo "$NGROK_PID" > .ngrok.pid

# Show services are running
echo "========================================"
echo "  Services Running"
echo "========================================"
echo ""
echo "Press Ctrl+C to stop all services..."
echo ""
echo "Monitoring services... (checking every 10 seconds)"
echo ""

# Trap Ctrl+C
trap "echo ''; echo '⏹️  Stopping services...'; kill $BACKEND_PID $NGROK_PID 2>/dev/null; rm -f .backend.pid .ngrok.pid; echo '✅ Services stopped'; exit 0" INT TERM

# Monitor services
CHECK_COUNT=0
while true; do
    sleep 10
    CHECK_COUNT=$((CHECK_COUNT + 1))
    
    # Check if backend is still running
    if ! kill -0 $BACKEND_PID 2>/dev/null; then
        echo "❌ Backend stopped unexpectedly!"
        echo "   Last 20 lines of log:"
        tail -n 20 backend/backend.log
        kill $NGROK_PID 2>/dev/null
        rm -f .backend.pid .ngrok.pid
        exit 1
    fi
    
    # Check if ngrok is still running
    if ! kill -0 $NGROK_PID 2>/dev/null; then
        echo "❌ ngrok stopped unexpectedly!"
        echo "   Last 20 lines of log:"
        tail -n 20 ngrok.log
        kill $BACKEND_PID 2>/dev/null
        rm -f .backend.pid .ngrok.pid
        exit 1
    fi
    
    # Every minute, show a heartbeat
    if [ $((CHECK_COUNT % 6)) -eq 0 ]; then
        echo "💚 Services running ($(date '+%H:%M:%S')) - Backend: ✓  ngrok: ✓"
    fi
done
