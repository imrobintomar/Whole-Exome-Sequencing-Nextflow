# Firebase Authentication - Quick Start

## ✅ What's Done

Firebase Authentication has been fully integrated into your WES Pipeline!

### Frontend Changes
- ✅ Firebase SDK installed
- ✅ Login/Register now use Firebase Auth
- ✅ Automatic token management
- ✅ User state persists across page refreshes

### Backend Changes  
- ✅ Firebase Admin SDK for token verification
- ✅ Database updated with `firebase_uid` field
- ✅ JWT auth removed
- ✅ Automatic user creation on first login

## 🚀 Next Steps to Get Running

### 1. Get Firebase Service Account Key (5 minutes)

1. Visit: https://console.firebase.google.com/project/variant-ac1c6/settings/serviceaccounts/adminsdk
2. Click **"Generate new private key"**
3. Download the JSON file
4. Save it as: `/media/drprabudh/m3/Nextflow-Script/WholeExome/backend/firebase-service-account.json`

### 2. Update Backend Config

Edit `backend/firebase_auth.py` line 13:

```python
SERVICE_ACCOUNT_PATH = os.getenv(
    "GOOGLE_APPLICATION_CREDENTIALS",
    "/media/drprabudh/m3/Nextflow-Script/WholeExome/backend/firebase-service-account.json"  # ← Update this path
)
```

### 3. Enable Email/Password Auth in Firebase

1. Visit: https://console.firebase.google.com/project/variant-ac1c6/authentication/providers
2. Click on **"Email/Password"**
3. Toggle **Enable** and **Save**

### 4. Install Dependencies & Run

**Backend:**
```bash
cd /media/drprabudh/m3/Nextflow-Script/WholeExome/backend

# Remove old virtual env
rm -rf venv

# Run setup script
./run.sh
```

**Frontend:**
```bash
cd /media/drprabudh/m3/Nextflow-Script/WholeExome/frontend

# Install dependencies (includes firebase)
npm install

# Start dev server
npm run dev
```

### 5. Test It!

1. Open: http://localhost:3000
2. Click **"Register here"**
3. Create account with email/password
4. You should be automatically logged in
5. Upload FASTQ files - it works!

## 🔥 What You Get

1. **Secure Authentication**: Industry-standard Firebase Auth
2. **Easy Password Reset**: Built-in email password reset
3. **No Token Management**: Firebase handles everything
4. **Email Verification**: Can be enabled in Firebase Console
5. **Social Login**: Easy to add Google, GitHub, etc. later

## ⚙️ How It Works

```
User enters email/password
  ↓
Firebase authenticates
  ↓
Frontend gets ID token
  ↓
All API calls include: Authorization: Bearer {firebase-token}
  ↓
Backend verifies token with Firebase
  ↓
Backend finds/creates user in database
  ↓
Request processed
```

## 🐛 Common Issues

### Issue: "Invalid authentication credentials"
**Solution**: Make sure you downloaded the service account JSON and updated the path in `firebase_auth.py`

### Issue: "Email already in use"
**Solution**: User already registered. Use different email or reset password.

### Issue: Backend fails to start
```bash
# Check if firebase-admin is installed
pip list | grep firebase

# Should show: firebase-admin==6.6.0
# If not, install:
pip install firebase-admin==6.6.0
```

### Issue: Frontend shows Firebase error
**Solution**: Check browser console for specific error. Make sure Email/Password is enabled in Firebase Console.

## 📚 More Info

For detailed setup and troubleshooting: see [FIREBASE_SETUP.md](./FIREBASE_SETUP.md)

## 🎯 Summary

You've successfully upgraded from basic JWT to enterprise-grade Firebase Authentication! Users can now securely register and login to use your WES Pipeline.

**Next**: Deploy frontend to Vercel and backend to your server for production use!
