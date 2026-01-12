#!/bin/bash
# Script to restart the backend server

echo "🔄 Restarting backend server..."

# Kill any existing uvicorn processes
echo "Stopping existing backend processes..."
pkill -f "uvicorn main:app" 2>/dev/null

# Wait a moment
sleep 2

# Start the backend server
echo "Starting backend server..."
cd "$(dirname "$0")/backend" || exit 1

# Check if virtual environment exists
if [ -d "venv" ]; then
    echo "Activating virtual environment..."
    source venv/bin/activate
fi

# Start uvicorn with reload
echo "Starting uvicorn..."
nohup uvicorn main:app --reload --host 0.0.0.0 --port 8000 > ../backend.log 2>&1 &

# Get the PID
BACKEND_PID=$!
echo "✅ Backend started with PID: $BACKEND_PID"

# Wait a moment for it to start
sleep 3

# Check if it's running
if ps -p $BACKEND_PID > /dev/null; then
    echo "✅ Backend is running successfully!"
    echo "📝 Logs are being written to: backend.log"
    echo "🌐 API available at: http://localhost:8000"
    echo "📚 API docs at: http://localhost:8000/docs"

    # Verify the new endpoints are available
    echo ""
    echo "🔍 Verifying new endpoints..."
    sleep 2

    if curl -s http://localhost:8000/openapi.json | grep -q "users.*details"; then
        echo "✅ User details endpoints are loaded!"
    else
        echo "⚠️  Warning: User details endpoints not found. Check backend.log"
    fi
else
    echo "❌ Failed to start backend. Check backend.log for errors."
    exit 1
fi
