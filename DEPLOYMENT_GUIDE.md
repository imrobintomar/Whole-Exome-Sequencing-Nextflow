# Deployment Guide: Vercel Frontend + Local Backend

This guide will help you deploy the frontend to Vercel while keeping the backend running on your local system.

## Overview

- **Frontend**: Deployed to Vercel (public access)
- **Backend**: Running on your local machine (exposed via ngrok)

## Prerequisites

1. Vercel account (free tier works)
2. ngrok account (free tier works)
3. Git repository (GitHub, GitLab, or Bitbucket)

---

## Part 1: Expose Backend via ngrok

### Step 1: Install ngrok

```bash
# Download and install ngrok
# Visit: https://ngrok.com/download

# Or on Linux:
curl -s https://ngrok-agent.s3.amazonaws.com/ngrok.asc | sudo tee /etc/apt/trusted.gpg.d/ngrok.asc >/dev/null
echo "deb https://ngrok-agent.s3.amazonaws.com buster main" | sudo tee /etc/apt/sources.list.d/ngrok.list
sudo apt update
sudo apt install ngrok
```

### Step 2: Authenticate ngrok

```bash
# Sign up at https://dashboard.ngrok.com/signup
# Get your auth token from https://dashboard.ngrok.com/get-started/your-authtoken

ngrok config add-authtoken YOUR_AUTH_TOKEN
```

### Step 3: Start Backend Server

Make sure your backend is running:

```bash
cd /media/drprabudh/m3/Nextflow-Script/WholeExome/backend
source venv/bin/activate
python main.py
```

### Step 4: Expose Backend with ngrok

In a new terminal:

```bash
ngrok http 8000
```

You'll see output like:
```
Forwarding   https://abcd-1234-5678.ngrok-free.app -> http://localhost:8000
```

**Copy the HTTPS URL** (e.g., `https://abcd-1234-5678.ngrok-free.app`)

### Step 5: Update Backend CORS Settings

Edit `/media/drprabudh/m3/Nextflow-Script/WholeExome/backend/.env`:

```bash
# Add your Vercel domain (you'll get this after deploying)
CORS_ORIGINS=http://localhost:3000,https://your-app.vercel.app
```

Or update `config.py`:

```python
CORS_ORIGINS: str = "http://localhost:3000,https://your-app.vercel.app,https://*.vercel.app"
```

Restart the backend after updating.

---

## Part 2: Deploy Frontend to Vercel

### Option A: Deploy via Vercel Dashboard (Easiest)

1. **Push code to GitHub**:
   ```bash
   cd /media/drprabudh/m3/Nextflow-Script/WholeExome/frontend
   git init
   git add .
   git commit -m "Initial commit for Vercel deployment"

   # Create a new repository on GitHub, then:
   git remote add origin https://github.com/YOUR_USERNAME/wes-pipeline-frontend.git
   git branch -M main
   git push -u origin main
   ```

2. **Go to Vercel Dashboard**:
   - Visit https://vercel.com/
   - Click "Add New Project"
   - Import your GitHub repository

3. **Configure Environment Variables**:
   - In Vercel dashboard, go to Settings → Environment Variables
   - Add: `NEXT_PUBLIC_API_URL` = `https://your-ngrok-url.ngrok-free.app`
   - Make sure to include the ngrok URL you got from Step 4 above

4. **Deploy**:
   - Click "Deploy"
   - Wait for build to complete
   - Copy your Vercel URL (e.g., `https://wes-pipeline.vercel.app`)

5. **Update CORS in Backend**:
   - Update backend `.env` or `config.py` with your Vercel URL
   - Restart backend server

### Option B: Deploy via Vercel CLI

1. **Install Vercel CLI**:
   ```bash
   npm install -g vercel
   ```

2. **Login to Vercel**:
   ```bash
   vercel login
   ```

3. **Deploy Frontend**:
   ```bash
   cd /media/drprabudh/m3/Nextflow-Script/WholeExome/frontend
   vercel
   ```

4. **Follow prompts**:
   - Set up and deploy? Yes
   - Scope: Your account
   - Link to existing project? No
   - Project name: wes-pipeline-frontend
   - Directory: ./
   - Override settings? No

5. **Set Environment Variable**:
   ```bash
   vercel env add NEXT_PUBLIC_API_URL
   # When prompted, enter your ngrok URL: https://your-ngrok-url.ngrok-free.app
   # Select: Production, Preview, Development
   ```

6. **Deploy to Production**:
   ```bash
   vercel --prod
   ```

---

## Part 3: Testing the Deployment

1. **Access your Vercel URL** (e.g., `https://wes-pipeline.vercel.app`)
2. **Login with your Firebase account**
3. **Test features**:
   - Submit a job (file upload)
   - View job list
   - View variant analysis
   - View ACMG classification
   - View IGV browser

---

## Important Notes

### ngrok Free Tier Limitations

- **URL changes on restart**: Every time you restart ngrok, you get a new URL
- **Solution**:
  - Upgrade to ngrok paid plan for a static domain
  - OR: Update `NEXT_PUBLIC_API_URL` in Vercel every time ngrok restarts
  - OR: Use a static domain setup (see below)

### Static Domain Setup (Recommended for Production)

If you have a static IP or domain:

1. **Use your public IP**:
   ```bash
   # Get your public IP
   curl ifconfig.me

   # Use: http://YOUR_PUBLIC_IP:8000
   ```

2. **Configure port forwarding** on your router:
   - Forward port 8000 to your machine's local IP
   - Use: `http://YOUR_PUBLIC_IP:8000`

3. **Use a domain** (if you have one):
   - Point domain A record to your public IP
   - Use: `https://api.yourdomain.com`

### Keeping Backend Running 24/7

Use a process manager like `systemd` or `pm2`:

```bash
# Install PM2
npm install -g pm2

# Start backend with PM2
cd /media/drprabudh/m3/Nextflow-Script/WholeExome/backend
source venv/bin/activate
pm2 start main.py --name wes-backend --interpreter python

# Save PM2 process list
pm2 save

# Setup PM2 to start on boot
pm2 startup
```

### Firebase Configuration

Your Firebase config is already public in the code (this is normal for Firebase):
- API keys in frontend code are NOT secrets
- Firebase security is handled by Security Rules
- Make sure your Firebase Security Rules are properly configured

---

## Updating the Deployment

### Update Frontend
```bash
cd /media/drprabudh/m3/Nextflow-Script/WholeExome/frontend
git add .
git commit -m "Update features"
git push origin main
# Vercel will auto-deploy
```

### Update Backend
```bash
cd /media/drprabudh/m3/Nextflow-Script/WholeExome/backend
git pull  # or make your changes
# Restart backend server
```

---

## Troubleshooting

### CORS Errors
- Make sure Vercel URL is in backend `CORS_ORIGINS`
- Check browser console for exact error
- Verify ngrok is running and URL is correct

### Authentication Errors
- Verify Firebase config is correct
- Check Firebase console for authentication logs
- Ensure Firebase service account JSON is in backend

### File Upload Issues
- Check backend logs: `tail -f backend.log`
- Verify upload directory permissions
- Check available disk space

### IGV Browser Not Loading
- Ensure BAM/VCF files are accessible
- Check backend download endpoints
- Verify Firebase auth token is being sent

---

## Security Considerations

1. **Never commit secrets**:
   - `.env` files should be in `.gitignore`
   - Firebase service account JSON should NOT be in git
   - Use Vercel environment variables for secrets

2. **Enable HTTPS**:
   - Vercel provides HTTPS automatically
   - ngrok provides HTTPS in free tier
   - Use HTTPS for production

3. **Limit CORS origins**:
   - Don't use `*` in CORS_ORIGINS
   - Only allow your Vercel domain

4. **Secure Firebase**:
   - Configure Firebase Security Rules
   - Enable Firebase Authentication
   - Review Firebase Console logs

---

## Cost Estimates

- **Vercel**: Free tier (sufficient for this app)
- **ngrok**: Free tier (URL changes on restart)
- **ngrok Pro**: $8/month (static domain)
- **Firebase**: Free tier (generous limits)
- **Backend hosting**: Your local machine (free)

---

## Alternative: Deploy Backend to Cloud

If you want to deploy backend to cloud instead of running locally:

### Options:
1. **Railway.app** (easy, good free tier)
2. **Render.com** (easy, good free tier)
3. **DigitalOcean** (VPS, $4-6/month)
4. **AWS EC2** (complex, free tier available)
5. **Google Cloud Run** (serverless, pay-per-use)

Would you like a guide for any of these backend hosting options?
