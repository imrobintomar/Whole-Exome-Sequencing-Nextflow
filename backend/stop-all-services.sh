#!/bin/bash
# Stop all ATGCFLOW services

echo "🛑 Stopping ATGCFLOW Services"
echo "============================="
echo ""

echo "Stopping backend service..."
sudo systemctl stop atgcflow-backend.service

echo "Stopping ngrok service..."
sudo systemctl stop atgcflow-ngrok.service

echo ""
echo "Checking status..."
echo ""

echo "Backend Service:"
sudo systemctl status atgcflow-backend.service --no-pager | head -5

echo ""
echo "ngrok Service:"
sudo systemctl status atgcflow-ngrok.service --no-pager | head -5

echo ""
echo "✅ All services stopped"
