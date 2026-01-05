# Vercel Deployment - Ready to Deploy! ✅

All TypeScript errors have been fixed. Your app is ready for Vercel deployment.

## Your ngrok URL (Active Now)
```
https://disbursable-ennoblingly-jerlene.ngrok-free.dev
```

## Backend Status
✅ Running on: http://localhost:8000  
✅ Public URL: https://disbursable-ennoblingly-jerlene.ngrok-free.dev  
✅ CORS configured for ngrok URL

## Quick Deploy to Vercel

### Method 1: Vercel Dashboard (Recommended)

1. **Commit and Push to GitHub**:
   ```bash
   cd /media/drprabudh/m3/Nextflow-Script/WholeExome/frontend
   git init
   git add .
   git commit -m "Initial deployment"
   
   # Create repository on GitHub, then:
   git remote add origin https://github.com/YOUR_USERNAME/wes-pipeline-frontend.git
   git branch -M main
   git push -u origin main
   ```

2. **Deploy on Vercel**:
   - Go to: https://vercel.com/new
   - Click "Import Git Repository"
   - Select your repository
   - Click "Import"
   
3. **Add Environment Variable**:
   - In deployment settings, add:
     - **Name**: `NEXT_PUBLIC_API_URL`
     - **Value**: `https://disbursable-ennoblingly-jerlene.ngrok-free.dev`
   - Click "Deploy"

4. **After Deployment**:
   - Copy your Vercel URL (e.g., `https://wes-pipeline.vercel.app`)
   - Update backend CORS:
     ```bash
     # Edit backend/.env
     CORS_ORIGINS=http://localhost:3000,https://disbursable-ennoblingly-jerlene.ngrok-free.dev,https://wes-pipeline.vercel.app
     ```
   - Restart backend:
     ```bash
     cd /media/drprabudh/m3/Nextflow-Script/WholeExome
     pkill -f "python main.py"
     cd backend && source venv/bin/activate && nohup python main.py > backend.log 2>&1 &
     ```

### Method 2: Vercel CLI

```bash
# Install CLI
npm install -g vercel

# Login
vercel login

# Deploy
cd /media/drprabudh/m3/Nextflow-Script/WholeExome/frontend
vercel

# Add environment variable
vercel env add NEXT_PUBLIC_API_URL
# Enter: https://disbursable-ennoblingly-jerlene.ngrok-free.dev

# Deploy to production
vercel --prod
```

## Testing After Deployment

Visit your Vercel URL and test:
- ✅ Login with Firebase
- ✅ Upload FASTQ files
- ✅ View jobs
- ✅ View variant analysis charts
- ✅ View ACMG classification
- ✅ Open IGV browser

## Important Notes

### ngrok URL Changes
⚠️ Your ngrok URL will change every time you restart ngrok.

**When ngrok restarts**:
1. Get new URL from ngrok
2. Update in Vercel: Settings → Environment Variables → `NEXT_PUBLIC_API_URL`
3. Redeploy: Deployments → Click "..." → Redeploy
4. Update backend CORS in `.env`

**To avoid this**: Upgrade to ngrok Pro ($8/month) for a static domain.

### Keep Backend Running 24/7

Use PM2 for automatic restart:
```bash
npm install -g pm2
cd /media/drprabudh/m3/Nextflow-Script/WholeExome/backend
source venv/bin/activate
pm2 start main.py --name wes-backend --interpreter python
pm2 save
pm2 startup  # Follow instructions
```

### Monitoring

```bash
# Backend logs
tail -f /media/drprabudh/m3/Nextflow-Script/WholeExome/backend/backend.log

# ngrok status
curl http://localhost:4040/api/tunnels | jq

# Check if backend is up
curl http://localhost:8000/
```

## Troubleshooting

### CORS Error in Browser
- Check backend `.env` has correct Vercel URL
- Restart backend after changing `.env`
- Check browser console for exact origin

### 502 Bad Gateway
- Verify backend is running: `curl http://localhost:8000/`
- Check ngrok is running: `curl http://localhost:4040/api/tunnels`
- Restart both if needed

### Build Failed on Vercel
- All TypeScript errors are now fixed ✅
- Build should succeed
- If not, check Vercel logs for details

## Files Changed for Deployment

✅ Fixed files:
- `frontend/components/VariantVisualization.tsx` - Fixed Legend and Pie chart TypeScript errors
- `frontend/components/theme-provider.tsx` - Fixed ThemeProviderProps import
- `backend/.env` - Added ngrok URL to CORS

✅ Created files:
- `frontend/vercel.json` - Vercel configuration
- `frontend/.env.production` - Production environment template
- `frontend/.gitignore` - Prevents committing secrets
- Various deployment guides

## Cost Summary

- **Vercel**: Free tier (sufficient for this app)
- **Firebase**: Free tier (generous limits)
- **ngrok Free**: $0/month (URL changes)
- **ngrok Pro**: $8/month (static domain)
- **Backend**: Your local machine (free)

## Next Steps After Deployment

1. Test all features on your Vercel URL
2. Consider ngrok Pro for static URL
3. Set up PM2 for backend auto-restart
4. Configure Firebase Security Rules
5. Monitor usage in Vercel and Firebase dashboards

---

**Need Help?**
- Full guide: [DEPLOYMENT_GUIDE.md](DEPLOYMENT_GUIDE.md)
- Quick start: [QUICK_DEPLOY.md](QUICK_DEPLOY.md)
- Checklist: [DEPLOYMENT_CHECKLIST.md](DEPLOYMENT_CHECKLIST.md)
