# Start Services - Simple Manual Commands

The automated script had issues. Use these manual commands instead:

## Terminal 1: Start Backend

```bash
cd /media/drprabudh/m3/Nextflow-Script/WholeExome/backend
source venv/bin/activate
python main.py
```

Keep this terminal open. You should see:
```
INFO:     Uvicorn running on http://0.0.0.0:8000
✅ Firebase initialized...
```

## Terminal 2: Start ngrok

Open a NEW terminal window:

```bash
cd /media/drprabudh/m3/Nextflow-Script/WholeExome
ngrok http 8000
```

You'll see:
```
Forwarding   https://something.ngrok-free.dev -> http://localhost:8000
```

**Copy the HTTPS URL!**

## Terminal 3: Update CORS and Verify

Open a NEW terminal window:

```bash
cd /media/drprabudh/m3/Nextflow-Script/WholeExome

# Get your ngrok URL
curl -s http://localhost:4040/api/tunnels | grep -o '"public_url":"https://[^"]*' | grep -o 'https://[^"]*' | head -n 1

# Test backend
curl http://localhost:8000/

# Test ngrok
curl https://YOUR-NGROK-URL.ngrok-free.dev/
```

## Update Backend CORS

Edit `backend/.env`:
```bash
CORS_ORIGINS=http://localhost:3000,https://YOUR-NGROK-URL.ngrok-free.dev
```

Then restart the backend (in Terminal 1):
- Press Ctrl+C
- Run: `python main.py`

## Deploy to Vercel

Now you can deploy to Vercel with:
- Environment Variable: `NEXT_PUBLIC_API_URL` = `https://YOUR-NGROK-URL.ngrok-free.dev`

After deploying, add your Vercel URL to backend `.env`:
```bash
CORS_ORIGINS=http://localhost:3000,https://YOUR-NGROK-URL.ngrok-free.dev,https://your-app.vercel.app
```

## Stop Services

- Terminal 1 (Backend): Press Ctrl+C
- Terminal 2 (ngrok): Press Ctrl+C

## Quick Commands Reference

```bash
# Check if backend is running
curl http://localhost:8000/

# Get ngrok URL
curl -s http://localhost:4040/api/tunnels | grep public_url

# View backend logs
tail -f backend/backend.log

# Kill all (if needed)
pkill -f "python main.py"
pkill -f ngrok
```
