# Quick Start: Keep ATGCFLOW Running 24/7

## TL;DR - One Command Setup

To keep your backend running 24/7 without VS Code:

```bash
cd /media/drprabudh/m3/Nextflow-Script/WholeExome/backend
./start-backend-service.sh
```

That's it! Your backend will now run permanently.

---

## What This Does

✅ **Frontend (https://atgcflow.com)** - Already running 24/7 on Vercel
✅ **Backend** - Will run 24/7 as a system service

The script will:
1. Create log directories
2. Install a systemd service
3. Enable automatic startup on boot
4. Start the backend immediately
5. Configure automatic restart on failure

---

## After Setup

### Check if it's running:
```bash
sudo systemctl status atgcflow-backend.service
```

You should see `active (running)` in green.

### View live logs:
```bash
# System logs
sudo journalctl -u atgcflow-backend.service -f

# Or application logs
tail -f /media/drprabudh/m3/Nextflow-Script/WholeExome/backend/logs/backend.log
```

### Common commands:
```bash
# Restart backend
sudo systemctl restart atgcflow-backend.service

# Stop backend
sudo systemctl stop atgcflow-backend.service

# Start backend
sudo systemctl start atgcflow-backend.service
```

---

## What About ngrok?

Your backend is currently exposed via ngrok at:
```
https://disbursable-ennoblingly-jerlene.ngrok-free.dev
```

**Important:** ngrok URLs change each time you restart ngrok. You have two options:

### Option 1: Keep using ngrok (Current setup)
Each time you restart ngrok, you need to:
1. Start ngrok: `ngrok http 8000`
2. Copy the new URL
3. Update `.env` file CORS_ORIGINS with the new URL
4. Restart backend: `sudo systemctl restart atgcflow-backend.service`

### Option 2: Get a permanent URL (Recommended)

**Use Cloudflare Tunnel (Free & Better than ngrok):**

1. Install:
```bash
wget https://github.com/cloudflare/cloudflared/releases/latest/download/cloudflared-linux-amd64.deb
sudo dpkg -i cloudflared-linux-amd64.deb
```

2. Setup:
```bash
cloudflared tunnel login
cloudflared tunnel create atgcflow-backend
cloudflared tunnel route dns atgcflow-backend api.atgcflow.com
```

3. Run tunnel (in background):
```bash
cloudflared tunnel run --url localhost:8000 atgcflow-backend
```

4. Your backend will be available at: `https://api.atgcflow.com`

See [backend/DEPLOYMENT-GUIDE.md](backend/DEPLOYMENT-GUIDE.md) for more options.

---

## Verification Checklist

✅ Frontend at https://atgcflow.com - Working
✅ Backend service running - Check with: `sudo systemctl status atgcflow-backend.service`
✅ Backend accessible via ngrok - Test: Visit your ngrok URL/docs
✅ Email verification working - Verified in previous steps

---

## You Can Now Close VS Code!

Your application will continue running:
- Frontend: Hosted on Vercel (always online)
- Backend: Running as a system service (always online)
- Even if you restart your computer, the backend will start automatically

---

## Troubleshooting

**If backend is not running:**
```bash
# Check status
sudo systemctl status atgcflow-backend.service

# Check logs for errors
sudo journalctl -u atgcflow-backend.service -n 50

# Restart
sudo systemctl restart atgcflow-backend.service
```

**If you get "Failed to start":**
- Check that port 8000 is not already in use: `sudo lsof -i :8000`
- Check logs: `tail -f backend/logs/backend-error.log`
- Verify Python dependencies: `pip install -r backend/requirements.txt`

For detailed deployment options, see [backend/DEPLOYMENT-GUIDE.md](backend/DEPLOYMENT-GUIDE.md).
