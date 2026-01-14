#!/bin/bash
# MinIO Startup Script
# Created: 2026-01-14

set -e

# Load configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/minio-config.env"

# Ensure data directory exists
mkdir -p "${MINIO_DATA_DIR}"

echo "=========================================="
echo "Starting MinIO Server"
echo "=========================================="
echo "Data Directory: ${MINIO_DATA_DIR}"
echo "API Endpoint:   http://localhost:9000"
echo "Console:        http://localhost:9001"
echo "User:           ${MINIO_ROOT_USER}"
echo "=========================================="
echo ""

# Export environment variables
export MINIO_ROOT_USER
export MINIO_ROOT_PASSWORD
export MINIO_ADDRESS
export MINIO_CONSOLE_ADDRESS
export MINIO_REGION
export MINIO_BROWSER
export MINIO_PROMETHEUS_AUTH_TYPE
export MINIO_LICENSE

# Start MinIO
exec ~/bin/minio server "${MINIO_DATA_DIR}"
