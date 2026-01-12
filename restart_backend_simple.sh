#!/bin/bash
# Simple script to restart the backend server

echo "🔄 Restarting backend server..."

# Kill any existing uvicorn processes
echo "Stopping existing backend processes..."
pkill -f "uvicorn main:app" 2>/dev/null
sleep 2

# Navigate to backend directory
cd backend || exit 1

# Start uvicorn in background
echo "Starting backend server..."
python3 -m uvicorn main:app --reload --host 0.0.0.0 --port 8000 > ../backend.log 2>&1 &

BACKEND_PID=$!
echo "✅ Backend started with PID: $BACKEND_PID"
echo "📝 Logs: tail -f backend.log"

# Wait for server to start
sleep 3

# Test if it's responding
if curl -s http://localhost:8000/ | grep -q "ok"; then
    echo "✅ Backend is responding!"

    # Check for new endpoints
    if curl -s http://localhost:8000/openapi.json | grep -q "/admin/users/{user_uid}/details"; then
        echo "✅ New user details endpoints are loaded!"
    else
        echo "⚠️  User details endpoints not found yet. May need a moment to load."
    fi
else
    echo "❌ Backend not responding. Check backend.log"
fi
