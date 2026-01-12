#!/bin/bash
# Kill any running backend processes and restart via systemd

echo "🔄 Killing existing backend processes..."
pkill -f "uvicorn main:app" || echo "No running processes found"
sleep 2

echo "Starting backend via systemd..."
sudo systemctl restart atgcflow-backend.service
sudo systemctl restart atgcflow-ngrok.service

sleep 3

echo ""
echo "Backend Status:"
sudo systemctl status atgcflow-backend.service --no-pager | head -10

echo ""
echo "ngrok Status:"
sudo systemctl status atgcflow-ngrok.service --no-pager | head -10

echo ""
echo "Testing backend..."
curl -s http://localhost:8000/

echo ""
echo "✅ Services restarted"
