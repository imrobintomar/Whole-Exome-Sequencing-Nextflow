#!/bin/bash

# Stop Backend and ngrok Services

echo "Stopping WES Pipeline services..."

# Stop from PID files
if [ -f .backend.pid ]; then
    BACKEND_PID=$(cat .backend.pid)
    if kill -0 $BACKEND_PID 2>/dev/null; then
        kill $BACKEND_PID
        echo "✅ Stopped backend (PID: $BACKEND_PID)"
    fi
    rm .backend.pid
fi

if [ -f .ngrok.pid ]; then
    NGROK_PID=$(cat .ngrok.pid)
    if kill -0 $NGROK_PID 2>/dev/null; then
        kill $NGROK_PID
        echo "✅ Stopped ngrok (PID: $NGROK_PID)"
    fi
    rm .ngrok.pid
fi

# Fallback: kill by process name
pkill -f "python main.py" 2>/dev/null && echo "✅ Stopped backend process"
pkill -f "ngrok http" 2>/dev/null && echo "✅ Stopped ngrok process"

echo "Done!"
