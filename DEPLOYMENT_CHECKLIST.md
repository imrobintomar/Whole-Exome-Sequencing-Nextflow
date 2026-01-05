# Deployment Checklist ✅

Use this checklist to ensure everything is configured correctly.

## Before Deployment

### Backend Setup
- [ ] Backend dependencies installed (`pip install -r requirements.txt`)
- [ ] Firebase service account JSON downloaded to `backend/firebase-service-account.json`
- [ ] gnomAD constraint data downloaded (`bash backend/download_gnomad_constraint.sh`)
- [ ] Backend starts successfully (`python backend/main.py`)
- [ ] Backend health check works (`curl http://localhost:8000/`)

### Frontend Setup
- [ ] Frontend dependencies installed (`npm install`)
- [ ] Firebase config updated in `frontend/lib/firebase.ts`
- [ ] Frontend builds successfully (`npm run build`)
- [ ] Frontend runs locally (`npm run dev`)

## Deployment Steps

### 1. ngrok Setup
- [ ] ngrok installed (`sudo apt install ngrok` or from https://ngrok.com/download)
- [ ] ngrok authenticated (`ngrok config add-authtoken YOUR_TOKEN`)
- [ ] ngrok tunnel created (`./start-with-ngrok.sh`)
- [ ] ngrok URL copied (e.g., `https://abcd-1234.ngrok-free.app`)

### 2. Backend CORS Configuration
- [ ] Edit `backend/.env`
- [ ] Add ngrok URL to CORS_ORIGINS: `CORS_ORIGINS=http://localhost:3000,https://your-ngrok-url.ngrok-free.app`
- [ ] Backend restarted after CORS update

### 3. Frontend Git Setup
- [ ] Git repository initialized in `frontend/` directory
- [ ] `.gitignore` configured (prevents committing secrets)
- [ ] Code committed (`git add . && git commit -m "Initial commit"`)
- [ ] GitHub repository created
- [ ] Code pushed to GitHub (`git push -u origin main`)

### 4. Vercel Deployment
- [ ] Vercel account created (https://vercel.com/signup)
- [ ] New project created in Vercel dashboard
- [ ] GitHub repository imported
- [ ] Environment variable added:
  - Name: `NEXT_PUBLIC_API_URL`
  - Value: Your ngrok URL
- [ ] Deployment successful
- [ ] Vercel URL copied (e.g., `https://wes-pipeline.vercel.app`)

### 5. Final Backend CORS Update
- [ ] Edit `backend/.env` again
- [ ] Add Vercel URL to CORS_ORIGINS: `CORS_ORIGINS=http://localhost:3000,https://your-ngrok-url.ngrok-free.app,https://wes-pipeline.vercel.app`
- [ ] Backend restarted with new CORS settings

## Testing

### Local Testing (Before Deployment)
- [ ] Login works (`http://localhost:3000`)
- [ ] Can upload FASTQ files
- [ ] Can view job list
- [ ] Can view variant analysis
- [ ] Can view ACMG classification
- [ ] IGV browser loads

### Production Testing (After Deployment)
- [ ] Login works on Vercel URL
- [ ] Can upload FASTQ files
- [ ] Can view job list
- [ ] Can view variant analysis
- [ ] Can view ACMG classification
- [ ] IGV browser loads
- [ ] No CORS errors in browser console
- [ ] No authentication errors

## Security Checklist

- [ ] Firebase security rules configured
- [ ] Firebase service account JSON NOT committed to git
- [ ] `.env` files added to `.gitignore`
- [ ] Backend uses HTTPS (ngrok provides this)
- [ ] Frontend uses HTTPS (Vercel provides this)
- [ ] CORS origins restricted (no `*` wildcard)
- [ ] Strong SECRET_KEY set in backend `.env`

## Monitoring

- [ ] Backend logs accessible (`tail -f backend/backend.log`)
- [ ] ngrok dashboard accessible (`http://localhost:4040`)
- [ ] Vercel deployment logs accessible (Vercel dashboard)
- [ ] Firebase authentication logs reviewed

## Troubleshooting Commands

```bash
# Check if backend is running
curl http://localhost:8000/

# Check ngrok tunnel
curl http://localhost:4040/api/tunnels

# View backend logs
tail -f backend/backend.log

# Restart backend
./stop-services.sh
./start-with-ngrok.sh

# Test CORS
curl -H "Origin: https://wes-pipeline.vercel.app" \
  -H "Access-Control-Request-Method: POST" \
  -X OPTIONS \
  http://localhost:8000/jobs

# Check backend processes
ps aux | grep python
ps aux | grep ngrok
```

## Common Issues

### Issue: CORS Error
**Solution**: 
- Verify Vercel URL is in `backend/.env` CORS_ORIGINS
- Restart backend after changing `.env`
- Check browser console for exact origin

### Issue: ngrok URL Changed
**Solution**:
- Update `NEXT_PUBLIC_API_URL` in Vercel environment variables
- Redeploy frontend on Vercel
- Update `backend/.env` CORS_ORIGINS

### Issue: Authentication Failed
**Solution**:
- Check Firebase service account JSON exists in `backend/`
- Verify Firebase config in `frontend/lib/firebase.ts`
- Check Firebase console for authentication logs

### Issue: File Upload Fails
**Solution**:
- Check backend logs: `tail -f backend/backend.log`
- Verify upload directory exists and is writable
- Check available disk space: `df -h`

### Issue: IGV Browser Not Loading
**Solution**:
- Check if BAM/VCF files exist in results directory
- Verify download endpoints return 200 status
- Check browser console for errors
- Verify Firebase auth token is being sent

## Maintenance

### Daily
- [ ] Check backend logs for errors
- [ ] Monitor disk space usage
- [ ] Verify ngrok tunnel is active

### Weekly
- [ ] Review Firebase usage
- [ ] Review Vercel usage
- [ ] Check for completed jobs to archive

### Monthly
- [ ] Update dependencies (`npm update`, `pip install -U`)
- [ ] Review security logs
- [ ] Backup database

## Cost Monitoring

- **Vercel**: Check usage at https://vercel.com/dashboard/usage
- **Firebase**: Check usage at https://console.firebase.google.com/
- **ngrok**: Free tier limits at https://dashboard.ngrok.com/

## Need Help?

- Full deployment guide: [DEPLOYMENT_GUIDE.md](DEPLOYMENT_GUIDE.md)
- Quick start: [QUICK_DEPLOY.md](QUICK_DEPLOY.md)
- GitHub Issues: (if repository is public)
