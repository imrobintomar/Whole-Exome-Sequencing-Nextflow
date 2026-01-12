# Security Patch Instructions

## 🚨 Your Backend is Under Attack

You're experiencing an **IoT vulnerability scan** from a device on your local network (`192.168.1.27`).

## Immediate Actions (Do Now!)

### 1. Block the Attacker IP

```bash
cd /media/drprabudh/m3/Nextflow-Script/WholeExome/backend
./block-attacker.sh
```

This will immediately block IP `192.168.1.27` using iptables.

### 2. Identify the Attacking Device

```bash
# Find what device this is on your network
nmap -sn 192.168.1.0/24 | grep -B 2 "192.168.1.27"

# Or check your router's DHCP table
# Look for device with IP 192.168.1.27
```

**Possible culprits:**
- Compromised IoT camera
- Infected smart TV
- Compromised computer on your network
- Malware on a phone/tablet

### 3. Add Security Middleware to Backend

#### Step 1: Add security middleware to main.py

Edit `backend/main.py` and add these lines **after line 34** (after `app = FastAPI(...)`):

```python
# Add this import at the top with other imports
from security_middleware import SecurityMiddleware

# Add middleware BEFORE CORS (around line 51, before app.add_middleware(CORSMiddleware))
# This is critical - security middleware must run BEFORE CORS

# Security middleware - add this RIGHT AFTER app = FastAPI(...)
security_middleware = SecurityMiddleware(
    app,
    max_requests_per_minute=60,  # Max 60 requests per minute per IP
    ban_duration_minutes=15       # Ban for 15 minutes on attack
)
app.add_middleware(SecurityMiddleware,
    max_requests_per_minute=60,
    ban_duration_minutes=15
)
```

#### Step 2: Restart backend

```bash
sudo systemctl restart atgcflow-backend.service

# Or if using the old method
pkill -f "python main.py"
cd backend
source venv/bin/activate
python main.py &
```

---

## What the Security Middleware Does

1. **Rate Limiting**: Max 60 requests/minute per IP
2. **Attack Detection**: Blocks known IoT/camera exploit paths
3. **Auto-Ban**: Bans attacking IPs for 15 minutes
4. **Request Logging**: Logs all suspicious activity

### Protected Against:

- ✅ IoT camera exploits
- ✅ CGI injection attacks
- ✅ Path traversal attacks
- ✅ Brute force attempts
- ✅ Port scanning
- ✅ Credential harvesting

---

## Manual Patch Instructions

If you prefer to patch manually:

### Option 1: Quick Patch (Add to existing main.py)

1. **Add import** (line 22, after other imports):
```python
from security_middleware import SecurityMiddleware
```

2. **Add middleware** (line 52, BEFORE CORS middleware):
```python
# Security middleware (must be before CORS)
app.add_middleware(
    SecurityMiddleware,
    max_requests_per_minute=60,
    ban_duration_minutes=15
)
```

3. **Restart backend**:
```bash
sudo systemctl restart atgcflow-backend.service
```

### Option 2: Automatic Patch Script

Create `backend/apply-security-patch.sh`:

```bash
#!/bin/bash
set -e

echo "🔒 Applying security patch to backend..."

BACKEND_DIR="/media/drprabudh/m3/Nextflow-Script/WholeExome/backend"
MAIN_PY="$BACKEND_DIR/main.py"

# Backup original
cp "$MAIN_PY" "$MAIN_PY.backup.$(date +%Y%m%d_%H%M%S)"

# Check if already patched
if grep -q "SecurityMiddleware" "$MAIN_PY"; then
    echo "✅ Security middleware already installed"
else
    # Add import
    sed -i '/import pandas as pd/a from security_middleware import SecurityMiddleware' "$MAIN_PY"

    # Add middleware (before CORS)
    sed -i '/app.add_middleware($/i\# Security middleware\napp.add_middleware(\n    SecurityMiddleware,\n    max_requests_per_minute=60,\n    ban_duration_minutes=15\n)\n' "$MAIN_PY"

    echo "✅ Security patch applied"
fi

# Restart service
echo "🔄 Restarting backend..."
sudo systemctl restart atgcflow-backend.service

echo "✅ Backend secured and restarted"
echo ""
echo "Check logs: sudo journalctl -u atgcflow-backend.service -f"
```

---

## Verify Security is Working

### Test 1: Check backend is running
```bash
curl http://localhost:8000/
# Should return: {"status":"ok","message":"WES Pipeline API is running"}
```

### Test 2: Test attack detection (should be blocked)
```bash
curl http://localhost:8000/cgi-bin/test
# Should return: {"detail":"Not found"} with 404 status
```

### Test 3: Check logs for security events
```bash
sudo journalctl -u atgcflow-backend.service -f
# Look for 🚨 ATTACK DETECTED messages
```

---

## Long-Term Security Recommendations

### 1. Network Isolation

**Isolate IoT devices:**
- Put IP cameras and IoT devices on separate VLAN
- Don't allow IoT devices to access your work servers

### 2. Firewall Rules

**Only allow specific IPs to access backend:**
```bash
# Allow only your IP
sudo ufw allow from YOUR_IP to any port 8000

# Deny all others
sudo ufw deny 8000
```

### 3. Use HTTPS Only (ngrok already does this)

Your ngrok setup already provides HTTPS, which is good!

### 4. Enable IP Whitelisting (Optional)

If you want maximum security, edit `main.py` to only allow specific IPs:

```python
from security_middleware import IPWhitelistMiddleware

# Add after app creation
app.add_middleware(
    IPWhitelistMiddleware,
    allowed_ips=['192.168.1.100', '192.168.1.101']  # Your trusted IPs
)
```

### 5. Monitor Access Logs

```bash
# Watch for attacks
sudo journalctl -u atgcflow-backend.service -f | grep "ATTACK"

# See banned IPs
curl http://localhost:8000/admin/banned-ips  # (need to add endpoint)
```

---

## Investigating the Attack

### Find the attacking device:

```bash
# Method 1: Using nmap
sudo nmap -sn 192.168.1.0/24

# Method 2: Using arp
arp -a | grep 192.168.1.27

# Method 3: Check router DHCP table
# Log into router admin (usually 192.168.1.1)
# Check connected devices
```

### Common culprits:

1. **IP Camera** (Hikvision, Dahua, generic Chinese cameras)
   - Often come with backdoors or malware
   - Update firmware immediately
   - Change default passwords

2. **Smart TV**
   - May have malware from sketchy apps
   - Factory reset and update

3. **Compromised Computer**
   - Run antivirus scan
   - Check for botnet activity

4. **Router/Modem**
   - Someone may have gained access to your network
   - Change WiFi password
   - Update router firmware

---

## Emergency Network Lockdown

If attacks continue:

```bash
# 1. Block entire local subnet except your machine
sudo iptables -A INPUT -s 192.168.1.0/24 ! -s YOUR_IP -j DROP

# 2. Only allow backend access via localhost
# Edit atgcflow-backend.service:
ExecStart=/path/to/python -m uvicorn main:app --host 127.0.0.1 --port 8000

# 3. Use ngrok for external access only
# This way local network can't access backend directly
```

---

## Testing After Patch

```bash
# 1. Start backend with security
sudo systemctl restart atgcflow-backend.service

# 2. Monitor logs
sudo journalctl -u atgcflow-backend.service -f

# 3. Test normal access (should work)
curl http://localhost:8000/

# 4. Test attack path (should be blocked)
curl http://localhost:8000/cgi-bin/test

# 5. Test rate limiting
for i in {1..70}; do curl http://localhost:8000/ & done
# After 60 requests, should get: {"detail":"Too many requests"}
```

---

## Support

If you need help:

1. **Check logs**: `sudo journalctl -u atgcflow-backend.service -f`
2. **Verify middleware loaded**: Look for "Security middleware" in startup logs
3. **Test attack detection**: Try accessing `/cgi-bin/test` (should 404)

---

## Summary

1. ✅ Block attacker IP: `./block-attacker.sh`
2. ✅ Add security middleware to `main.py`
3. ✅ Restart backend: `sudo systemctl restart atgcflow-backend.service`
4. ✅ Find and isolate the attacking device
5. ✅ Monitor logs for future attacks

**Your backend will now automatically block similar attacks in the future!**
