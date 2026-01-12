#!/bin/bash
# Quick backend restart

echo "🔄 Restarting backend service..."
sudo systemctl restart atgcflow-backend.service

echo "Waiting for backend to start..."
sleep 3

echo ""
echo "Status:"
sudo systemctl status atgcflow-backend.service --no-pager | head -10

echo ""
echo "Testing backend..."
curl -s http://localhost:8000/ || echo "Backend not responding"

echo ""
echo "✅ Backend restarted"
