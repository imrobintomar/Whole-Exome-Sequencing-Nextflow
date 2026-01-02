#!/bin/bash

# WES Pipeline Backend Startup Script

echo "======================================"
echo "  WES Pipeline Backend Server"
echo "======================================"

# Check if virtual environment exists
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
fi

# Activate virtual environment
echo "Activating virtual environment..."
source venv/bin/activate

# Install/update dependencies
echo "Installing dependencies..."
pip install -q -r requirements.txt

# Check if .env exists
if [ ! -f ".env" ]; then
    echo "WARNING: .env file not found!"
    echo "Copying .env.example to .env..."
    cp .env.example .env
    echo "Please edit .env with your configuration before running again."
    exit 1
fi

# Create upload and results directories
echo "Creating required directories..."
source .env
mkdir -p "$UPLOAD_DIR" "$RESULTS_DIR"

# Start server
echo ""
echo "======================================"
echo "Starting FastAPI server..."
echo "API will be available at:"
echo "  - Local: http://localhost:8000"
echo "  - Docs: http://localhost:8000/docs"
echo "======================================"
echo ""

python main.py
