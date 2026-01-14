#!/bin/bash
# MinIO Stop Script
# Created: 2026-01-14

set -e

echo "Stopping MinIO server..."

# Find and kill MinIO process
if pgrep -f "minio server" > /dev/null; then
    pkill -f "minio server"
    echo "MinIO server stopped successfully"
else
    echo "MinIO server is not running"
fi
