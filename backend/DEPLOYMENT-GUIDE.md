# ATGCFLOW Backend Deployment Guide
## Keep Backend Running 24/7 Without VS Code

---

## Option 1: systemd Service (Recommended for Linux Server)

This is the best solution for production deployment on a Linux server. The backend will:
- Start automatically on system boot
- Restart automatically if it crashes
- Run in the background without any terminal/VS Code

### Step 1: Create Log Directory
```bash
cd /media/drprabudh/m3/Nextflow-Script/WholeExome/backend
mkdir -p logs
```

### Step 2: Update the Service File
Edit `atgcflow-backend.service` and replace `%USER%` with your actual username:
```bash
# Get your username
whoami

# Edit the service file
nano atgcflow-backend.service

# Replace %USER% with the output from whoami command
```

### Step 3: Install the Service
```bash
# Copy service file to systemd directory
sudo cp atgcflow-backend.service /etc/systemd/system/

# Reload systemd to recognize the new service
sudo systemctl daemon-reload

# Enable service to start on boot
sudo systemctl enable atgcflow-backend.service

# Start the service now
sudo systemctl start atgcflow-backend.service
```

### Step 4: Verify It's Running
```bash
# Check service status
sudo systemctl status atgcflow-backend.service

# Should show "active (running)" in green
```

### Step 5: View Logs
```bash
# View live logs
sudo journalctl -u atgcflow-backend.service -f

# Or check log files directly
tail -f logs/backend.log
tail -f logs/backend-error.log
```

### Managing the Service
```bash
# Stop the service
sudo systemctl stop atgcflow-backend.service

# Restart the service
sudo systemctl restart atgcflow-backend.service

# Check status
sudo systemctl status atgcflow-backend.service

# Disable auto-start on boot
sudo systemctl disable atgcflow-backend.service
```

---

## Option 2: Screen/tmux (Quick Alternative)

If you don't have sudo access or prefer a simpler solution:

### Using Screen:
```bash
# Install screen if not already installed
sudo apt install screen

# Start a new screen session
screen -S atgcflow-backend

# Navigate to backend directory
cd /media/drprabudh/m3/Nextflow-Script/WholeExome/backend

# Activate virtual environment if you have one
# source venv/bin/activate

# Start the backend
python3 -m uvicorn main:app --host 0.0.0.0 --port 8000

# Detach from screen: Press Ctrl+A, then D

# Reattach to screen later
screen -r atgcflow-backend

# List all screens
screen -ls

# Kill the screen session (to stop backend)
screen -S atgcflow-backend -X quit
```

### Using tmux:
```bash
# Install tmux if not already installed
sudo apt install tmux

# Start a new tmux session
tmux new -s atgcflow-backend

# Navigate to backend directory
cd /media/drprabudh/m3/Nextflow-Script/WholeExome/backend

# Start the backend
python3 -m uvicorn main:app --host 0.0.0.0 --port 8000

# Detach from tmux: Press Ctrl+B, then D

# Reattach to tmux later
tmux attach -t atgcflow-backend

# List all tmux sessions
tmux ls

# Kill the tmux session (to stop backend)
tmux kill-session -t atgcflow-backend
```

---

## Option 3: nohup (Simplest)

Simplest option, but no automatic restart on failure:

```bash
cd /media/drprabudh/m3/Nextflow-Script/WholeExome/backend

# Start backend in background
nohup python3 -m uvicorn main:app --host 0.0.0.0 --port 8000 > logs/backend.log 2>&1 &

# Get the process ID
echo $! > backend.pid

# View logs
tail -f logs/backend.log

# Stop the backend later
kill $(cat backend.pid)
```

---

## Option 4: Docker (Advanced)

For production deployment with isolation:

### Step 1: Create Dockerfile
```dockerfile
FROM python:3.10-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
```

### Step 2: Build and Run
```bash
# Build image
docker build -t atgcflow-backend .

# Run container
docker run -d \
  --name atgcflow-backend \
  --restart unless-stopped \
  -p 8000:8000 \
  -v /media/drprabudh/m3/Nextflow-Script/WholeExome/backend/.env:/app/.env \
  -v /media/drprabudh/m3/Nextflow-Script/WholeExome/backend/uploads:/app/uploads \
  -v /media/drprabudh/m3/Nextflow-Script/WholeExome/backend/results:/app/results \
  atgcflow-backend

# View logs
docker logs -f atgcflow-backend

# Stop container
docker stop atgcflow-backend

# Start container
docker start atgcflow-backend
```

---

## Recommended Setup for Your Use Case

Based on your setup, I recommend **Option 1 (systemd)** because:

1. ✅ Automatic restart on failure
2. ✅ Starts automatically on system boot
3. ✅ Easy log management
4. ✅ Standard Linux service management
5. ✅ No need to keep terminal open

---

## Exposing Backend to Internet

Your backend is currently running on `localhost:8000`. To make it accessible from the internet:

### Option A: ngrok (What you're currently using)
```bash
# Start ngrok in a separate terminal/screen
ngrok http 8000

# Update CORS_ORIGINS in .env with the new ngrok URL
# Note: ngrok URL changes on each restart (free plan)
```

### Option B: Cloudflare Tunnel (Recommended - Free & Persistent)
```bash
# Install cloudflared
wget https://github.com/cloudflare/cloudflared/releases/latest/download/cloudflared-linux-amd64.deb
sudo dpkg -i cloudflared-linux-amd64.deb

# Login to Cloudflare
cloudflared tunnel login

# Create a tunnel
cloudflared tunnel create atgcflow-backend

# Route the tunnel
cloudflared tunnel route dns atgcflow-backend api.atgcflow.com

# Run the tunnel
cloudflared tunnel run --url localhost:8000 atgcflow-backend
```

### Option C: Reverse Proxy with Nginx (Production)
If you have a VPS/server with public IP:

```nginx
# /etc/nginx/sites-available/atgcflow-backend
server {
    listen 80;
    server_name api.atgcflow.com;

    location / {
        proxy_pass http://localhost:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }
}
```

Then enable SSL with Let's Encrypt:
```bash
sudo apt install certbot python3-certbot-nginx
sudo certbot --nginx -d api.atgcflow.com
```

---

## Monitoring Backend Health

### Create a Health Check Script
```bash
#!/bin/bash
# health-check.sh

BACKEND_URL="http://localhost:8000"
LOG_FILE="/media/drprabudh/m3/Nextflow-Script/WholeExome/backend/logs/health-check.log"

response=$(curl -s -o /dev/null -w "%{http_code}" $BACKEND_URL)

if [ $response -eq 200 ]; then
    echo "$(date): Backend is healthy (HTTP $response)" >> $LOG_FILE
else
    echo "$(date): Backend is down (HTTP $response)" >> $LOG_FILE
    # Restart the service if using systemd
    sudo systemctl restart atgcflow-backend.service
fi
```

Add to crontab to run every 5 minutes:
```bash
chmod +x health-check.sh

# Edit crontab
crontab -e

# Add this line:
*/5 * * * * /media/drprabudh/m3/Nextflow-Script/WholeExome/backend/health-check.sh
```

---

## Summary

**For immediate deployment:**
1. Use systemd service (Option 1)
2. Set it up with the commands in "Step 1-5" above
3. Your backend will run 24/7, restart automatically, and survive reboots

**For internet access:**
- Continue using ngrok (update URL in .env after each restart)
- OR switch to Cloudflare Tunnel for persistent URL
- OR deploy to a VPS with Nginx reverse proxy

**Your frontend is already live at https://atgcflow.com** - no action needed!

Let me know which option you want to use and I can help you set it up!
