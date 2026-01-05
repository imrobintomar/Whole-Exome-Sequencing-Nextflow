# Quick Deploy Guide

Deploy your WES Pipeline frontend to Vercel in 5 minutes!

## Prerequisites

- ngrok account (free): https://dashboard.ngrok.com/signup
- Vercel account (free): https://vercel.com/signup
- GitHub account (free): https://github.com/signup

## Step 1: Install & Setup ngrok (2 minutes)

```bash
# Install ngrok
sudo apt install ngrok

# Authenticate (get token from https://dashboard.ngrok.com/get-started/your-authtoken)
ngrok config add-authtoken YOUR_AUTH_TOKEN
```

## Step 2: Start Backend with ngrok (1 minute)

```bash
cd /media/drprabudh/m3/Nextflow-Script/WholeExome
./start-with-ngrok.sh
```

You'll see output like:
```
🌐 Backend (Public): https://abcd-1234.ngrok-free.app
```

**Copy this URL** - you'll need it for Vercel!

## Step 3: Deploy Frontend to Vercel (2 minutes)

### Option A: Deploy via Vercel Dashboard (Easiest)

1. **Push to GitHub**:
   ```bash
   cd frontend
   git init
   git add .
   git commit -m "Deploy to Vercel"
   
   # Create repo on GitHub, then:
   git remote add origin https://github.com/YOUR_USERNAME/wes-pipeline.git
   git push -u origin main
   ```

2. **Deploy on Vercel**:
   - Go to https://vercel.com/new
   - Import your GitHub repository
   - Add environment variable:
     - Name: `NEXT_PUBLIC_API_URL`
     - Value: Your ngrok URL (e.g., `https://abcd-1234.ngrok-free.app`)
   - Click "Deploy"

3. **Update Backend CORS**:
   After deployment, copy your Vercel URL (e.g., `https://wes-pipeline.vercel.app`)
   
   Edit `backend/.env`:
   ```bash
   CORS_ORIGINS=http://localhost:3000,https://abcd-1234.ngrok-free.app,https://wes-pipeline.vercel.app
   ```
   
   Restart backend:
   ```bash
   ./stop-services.sh
   ./start-with-ngrok.sh
   ```

### Option B: Deploy via Vercel CLI

```bash
# Install Vercel CLI
npm install -g vercel

# Login
vercel login

# Deploy
cd frontend
vercel

# Set environment variable
vercel env add NEXT_PUBLIC_API_URL
# Enter your ngrok URL when prompted

# Deploy to production
vercel --prod
```

## Step 4: Test Your Deployment! 🎉

Visit your Vercel URL and test:
- ✅ Login with Firebase
- ✅ Submit a job
- ✅ View results
- ✅ View ACMG classification
- ✅ View IGV browser

## Important Notes

⚠️ **ngrok Free Tier**: URL changes every time you restart. Solutions:
- Upgrade to ngrok Pro ($8/month) for static domain
- Update `NEXT_PUBLIC_API_URL` in Vercel when ngrok restarts
- Use your public IP with port forwarding instead

🔄 **Redeploying Frontend**:
```bash
cd frontend
git add .
git commit -m "Update"
git push origin main
# Vercel auto-deploys
```

🛑 **Stop Services**:
```bash
./stop-services.sh
```

## Troubleshooting

### CORS Error
- Check backend `.env` has your Vercel URL in `CORS_ORIGINS`
- Restart backend after changing `.env`

### 404 Errors
- Verify ngrok is running: `curl http://localhost:4040/api/tunnels`
- Check backend logs: `tail -f backend/backend.log`

### Authentication Error
- Verify Firebase config in `frontend/lib/firebase.ts`
- Check Firebase service account JSON in `backend/`

## Need Help?

See full guide: [DEPLOYMENT_GUIDE.md](DEPLOYMENT_GUIDE.md)
