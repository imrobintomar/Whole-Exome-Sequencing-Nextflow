#!/bin/bash
# Automated security patch for ATGCFLOW backend

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN_PY="$SCRIPT_DIR/main.py"
BACKUP_FILE="$MAIN_PY.backup.$(date +%Y%m%d_%H%M%S)"

echo "🔒 ATGCFLOW Security Patch"
echo "=========================="
echo ""

# Step 1: Check if security_middleware.py exists
if [ ! -f "$SCRIPT_DIR/security_middleware.py" ]; then
    echo "❌ security_middleware.py not found!"
    echo "   This should have been created. Please check the files."
    exit 1
fi
echo "✅ Found security_middleware.py"

# Step 2: Backup main.py
echo ""
echo "📋 Creating backup of main.py..."
cp "$MAIN_PY" "$BACKUP_FILE"
echo "✅ Backup created: $BACKUP_FILE"

# Step 3: Check if already patched
if grep -q "from security_middleware import SecurityMiddleware" "$MAIN_PY"; then
    echo ""
    echo "ℹ️  Security middleware is already imported in main.py"

    if grep -q "app.add_middleware(SecurityMiddleware" "$MAIN_PY"; then
        echo "✅ Security middleware is already enabled"
        echo ""
        echo "Your backend is already protected!"
        echo ""
        read -p "Do you want to restart the backend anyway? (y/n) " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            echo ""
            echo "🔄 Restarting backend..."
            sudo systemctl restart atgcflow-backend.service
            sleep 2
            echo "✅ Backend restarted"
        fi
        exit 0
    fi
fi

# Step 4: Apply patch
echo ""
echo "🔧 Applying security patch to main.py..."

# Create temporary file with patched content
TEMP_FILE=$(mktemp)

# Flag to track if we've added the middleware
MIDDLEWARE_ADDED=false

while IFS= read -r line; do
    echo "$line" >> "$TEMP_FILE"

    # Add import after other imports (after pandas import)
    if [[ "$line" == *"import pandas as pd"* ]] && ! $MIDDLEWARE_ADDED; then
        echo "from security_middleware import SecurityMiddleware" >> "$TEMP_FILE"
        echo "   ✅ Added security middleware import"
    fi

    # Add middleware before CORS middleware
    if [[ "$line" == *"app.add_middleware("* ]] && [[ "$line" == *"CORSMiddleware"* ]]; then
        # Insert security middleware BEFORE CORS
        cat >> "$TEMP_FILE" << 'EOF'

# Security Protection (MUST be before CORS)
app.add_middleware(
    SecurityMiddleware,
    max_requests_per_minute=60,   # Rate limit: 60 requests/min per IP
    ban_duration_minutes=15        # Ban duration for attackers
)

EOF
        MIDDLEWARE_ADDED=true
        echo "   ✅ Added security middleware initialization"
    fi
done < "$MAIN_PY"

# Replace original with patched version
mv "$TEMP_FILE" "$MAIN_PY"

echo "✅ Security patch applied successfully"
echo ""

# Step 5: Verify patch
echo "🔍 Verifying patch..."
if grep -q "from security_middleware import SecurityMiddleware" "$MAIN_PY" && \
   grep -q "app.add_middleware(SecurityMiddleware" "$MAIN_PY"; then
    echo "✅ Patch verification successful"
else
    echo "❌ Patch verification failed!"
    echo "   Restoring from backup..."
    cp "$BACKUP_FILE" "$MAIN_PY"
    echo "   Backup restored. Please apply patch manually."
    exit 1
fi

# Step 6: Block attacker IP
echo ""
echo "🚫 Blocking attacker IP (192.168.1.27)..."
if [ -f "$SCRIPT_DIR/block-attacker.sh" ]; then
    bash "$SCRIPT_DIR/block-attacker.sh"
else
    echo "⚠️  block-attacker.sh not found, skipping IP block"
    echo "   Manually block: sudo iptables -I INPUT -s 192.168.1.27 -j DROP"
fi

# Step 7: Restart backend
echo ""
echo "🔄 Restarting backend service..."
sudo systemctl restart atgcflow-backend.service

echo ""
echo "⏳ Waiting for backend to start..."
sleep 3

# Step 8: Verify backend is running
if curl -s http://localhost:8000/ > /dev/null 2>&1; then
    echo "✅ Backend is running and responding"
else
    echo "⚠️  Backend may not be running properly"
    echo "   Check logs: sudo journalctl -u atgcflow-backend.service -n 50"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "🎉 Security Patch Applied Successfully!"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "✅ Security Features Enabled:"
echo "   • Rate limiting (60 requests/min per IP)"
echo "   • Attack pattern detection"
echo "   • Automatic IP banning (15 min)"
echo "   • Request logging"
echo ""
echo "🔒 Protected Against:"
echo "   • IoT camera exploits"
echo "   • CGI injection attacks"
echo "   • Path traversal"
echo "   • Brute force attempts"
echo ""
echo "📊 Monitoring:"
echo "   Watch logs: sudo journalctl -u atgcflow-backend.service -f"
echo "   Look for: 🚨 ATTACK DETECTED"
echo ""
echo "🧪 Test Security:"
echo "   # Normal request (should work)"
echo "   curl http://localhost:8000/"
echo ""
echo "   # Attack path (should be blocked)"
echo "   curl http://localhost:8000/cgi-bin/test"
echo ""
echo "📁 Backup saved: $BACKUP_FILE"
echo ""
echo "✅ Your backend is now protected!"
echo ""
