# Fix Vercel API URL - URGENT! 🚨

Your Vercel deployment is using a placeholder URL. Here's how to fix it:

## Your Actual URLs

- **ngrok URL**: `https://disbursable-ennoblingly-jerlene.ngrok-free.dev`
- **Vercel URL**: `https://exomeanalysis-q9d6ris9l-ripsscon-2025.vercel.app`

## Fix Steps

### 1. Update Vercel Environment Variable

Go to your Vercel project:
1. Visit: https://vercel.com/ripsscon-2025/exomeanalysis/settings/environment-variables
2. Find `NEXT_PUBLIC_API_URL`
3. Edit it to: `https://disbursable-ennoblingly-jerlene.ngrok-free.dev`
4. Click "Save"

### 2. Redeploy Vercel

After changing the environment variable:
1. Go to: https://vercel.com/ripsscon-2025/exomeanalysis/deployments
2. Click the "..." menu on the latest deployment
3. Click "Redeploy"
4. Wait for deployment to complete (~2 minutes)

### 3. Restart Backend

Your backend CORS is now updated. Restart it:

```bash
# Terminal 1 - If backend is running, stop it (Ctrl+C), then:
cd /media/drprabudh/m3/Nextflow-Script/WholeExome/backend
source venv/bin/activate
python main.py
```

### 4. Verify Everything Works

After redeployment, visit:
- https://exomeanalysis-q9d6ris9l-ripsscon-2025.vercel.app

And check:
- ✅ No CORS errors in console
- ✅ Can login
- ✅ Can see jobs list

## Quick Verification Commands

```bash
# Check backend is running
curl http://localhost:8000/

# Check ngrok tunnel
curl https://disbursable-ennoblingly-jerlene.ngrok-free.dev/

# Test CORS from command line
curl -H "Origin: https://exomeanalysis-q9d6ris9l-ripsscon-2025.vercel.app" \
     -H "Access-Control-Request-Method: GET" \
     -X OPTIONS \
     https://disbursable-ennoblingly-jerlene.ngrok-free.dev/jobs
```

If the CORS test shows `access-control-allow-origin` header, you're good! ✅

## Current Status

✅ Backend CORS updated to include:
- http://localhost:3000
- https://disbursable-ennoblingly-jerlene.ngrok-free.dev  
- https://exomeanalysis-q9d6ris9l-ripsscon-2025.vercel.app
- https://exomeanalysis.vercel.app

⏳ Need to:
1. Update Vercel environment variable
2. Redeploy on Vercel
3. Restart backend

## Troubleshooting

### Still getting CORS errors?
1. Make sure backend is restarted after `.env` change
2. Check backend logs: `tail -f backend/backend.log`
3. Verify Vercel env var is saved (not just edited)
4. Make sure you clicked "Redeploy" after changing env var

### ngrok URL changed?
If you restart ngrok, the URL will change. You'll need to:
1. Get new ngrok URL
2. Update Vercel env var
3. Update backend `.env`
4. Restart backend
5. Redeploy Vercel

Consider upgrading to ngrok Pro ($8/month) for a static URL to avoid this!
