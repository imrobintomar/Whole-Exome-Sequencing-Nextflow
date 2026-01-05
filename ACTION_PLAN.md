# 🎯 Action Plan - Fix Your Vercel Deployment

## Current Status

✅ **Backend**: Running with updated CORS  
✅ **ngrok**: Active at `https://disbursable-ennoblingly-jerlene.ngrok-free.dev`  
✅ **Vercel**: Deployed at `https://exomeanalysis-q9d6ris9l-ripsscon-2025.vercel.app`  
❌ **Problem**: Vercel is using wrong API URL (placeholder instead of real ngrok URL)

## The Issue

Your Vercel app is trying to connect to:
```
https://your-ngrok-url.ngrok-free.app  ❌ (This doesn't exist!)
```

But it should connect to:
```
https://disbursable-ennoblingly-jerlene.ngrok-free.dev  ✅ (Your real ngrok URL)
```

## Fix in 3 Steps (5 minutes)

### Step 1: Update Vercel Environment Variable

1. Go to: **https://vercel.com/ripsscon-2025/exomeanalysis/settings/environment-variables**

2. Find the environment variable `NEXT_PUBLIC_API_URL`

3. Click "Edit" or "..." menu

4. Change the value from:
   ```
   https://your-ngrok-url.ngrok-free.app
   ```
   
   To:
   ```
   https://disbursable-ennoblingly-jerlene.ngrok-free.dev
   ```

5. Click **Save**

### Step 2: Redeploy on Vercel

1. Go to: **https://vercel.com/ripsscon-2025/exomeanalysis/deployments**

2. Find the latest deployment (top of the list)

3. Click the "**...**" (three dots) menu

4. Click "**Redeploy**"

5. Confirm the redeploy

6. Wait 1-2 minutes for deployment to complete

### Step 3: Test Your App

Visit: **https://exomeanalysis-q9d6ris9l-ripsscon-2025.vercel.app**

Test these features:
- ✅ Login with Firebase (should work)
- ✅ No CORS errors in browser console (F12)
- ✅ Can see jobs list
- ✅ Can upload files
- ✅ Can view variant analysis

## Verification

Open browser console (F12) and check:

**Before fix**:
```
❌ CORS header 'Access-Control-Allow-Origin' missing
❌ https://your-ngrok-url.ngrok-free.app/jobs 404
```

**After fix**:
```
✅ No CORS errors
✅ Requests to https://disbursable-ennoblingly-jerlene.ngrok-free.dev/ succeed
```

## Services Status

Your backend and ngrok are running properly:

```bash
# Backend
curl http://localhost:8000/
# Should return: {"status":"ok","message":"WES Pipeline API is running"}

# ngrok  
curl https://disbursable-ennoblingly-jerlene.ngrok-free.dev/
# Should return the same message

# ngrok dashboard (optional)
http://localhost:4040
```

## CORS Configuration

Your backend is now configured to accept requests from:
- ✅ `http://localhost:3000` (local development)
- ✅ `https://disbursable-ennoblingly-jerlene.ngrok-free.dev` (ngrok)
- ✅ `https://exomeanalysis-q9d6ris9l-ripsscon-2025.vercel.app` (Vercel preview)
- ✅ `https://exomeanalysis.vercel.app` (Vercel production)

## Important Notes

### If ngrok URL Changes

⚠️ **ngrok free tier**: URL changes every restart!

When you restart ngrok:
1. Get new URL from ngrok
2. Update Vercel env var (Step 1)
3. Update `backend/.env` CORS_ORIGINS
4. Restart backend
5. Redeploy Vercel (Step 2)

**Recommendation**: Upgrade to ngrok Pro ($8/month) for static domain.

### Keep Services Running

**Backend**:
```bash
# Check if running
ps aux | grep "python main.py"

# Restart if needed
cd backend && source venv/bin/activate && python main.py
```

**ngrok**:
```bash
# Check if running
ps aux | grep ngrok

# Restart if needed
ngrok http 8000
```

## Troubleshooting

### Still Getting CORS Errors After Redeploy?

1. **Clear browser cache**: Ctrl+Shift+R or Cmd+Shift+R
2. **Check env var was saved**: Visit Vercel settings again
3. **Verify backend restarted**: `curl http://localhost:8000/`
4. **Check deployment used new var**: Look at deployment logs

### Can't Access Backend?

```bash
# Check backend is running
curl http://localhost:8000/

# If not, restart
lsof -ti:8000 | xargs kill -9
cd backend && source venv/bin/activate && python main.py &

# Check ngrok
curl http://localhost:4040/api/tunnels
```

## Next Steps After Fix

Once everything works:

1. ✅ Test all features thoroughly
2. 📝 Document your production URLs
3. 🔐 Set up Firebase Security Rules
4. 🔄 Consider PM2 for auto-restart: `npm install -g pm2 && pm2 start main.py --name wes-backend`
5. 💰 Consider ngrok Pro for static URL

---

**Need help?** Check other guides:
- [FIX_VERCEL_URL.md](FIX_VERCEL_URL.md) - Detailed fix instructions
- [START_SERVICES.md](START_SERVICES.md) - How to start/restart services
- [DEPLOYMENT_GUIDE.md](DEPLOYMENT_GUIDE.md) - Full deployment documentation
