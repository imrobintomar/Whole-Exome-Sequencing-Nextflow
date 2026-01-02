# 🚀 START HERE - Firebase Setup

Your WES Pipeline is **almost ready**! Just need to set up Firebase credentials.

## Current Status

✅ **Frontend**: Firebase SDK installed and configured
✅ **Backend**: Running on http://localhost:8000 (but auth won't work yet)
⚠️  **Firebase Auth**: Needs service account JSON file

## Quick Setup (5 minutes)

### Step 1: Download Firebase Service Account

1. **Visit**: https://console.firebase.google.com/project/variant-ac1c6/settings/serviceaccounts/adminsdk

2. **Click**: "Generate new private key" button

3. **Download**: The JSON file will download

4. **Save it as**: 
   ```
   /media/drprabudh/m3/Nextflow-Script/WholeExome/backend/firebase-service-account.json
   ```

### Step 2: Enable Email/Password Authentication

1. **Visit**: https://console.firebase.google.com/project/variant-ac1c6/authentication/providers

2. **Click** on "Email/Password"

3. **Toggle** "Enable" to ON

4. **Click** "Save"

### Step 3: Restart Backend

```bash
cd /media/drprabudh/m3/Nextflow-Script/WholeExome/backend

# Stop current backend (Ctrl+C if running)

# Start again
./run.sh
```

You should see:
```
✅ Firebase initialized with service account: ./firebase-service-account.json
```

### Step 4: Start Frontend

```bash
cd /media/drprabudh/m3/Nextflow-Script/WholeExome/frontend

npm run dev
```

### Step 5: Test It!

1. Open: http://localhost:3000

2. Click "Register here"

3. Create account:
   - Email: your@email.com
   - Display Name: Your Name (optional)
   - Password: (minimum 6 characters)

4. You should be logged in automatically!

5. Try uploading FASTQ files - it works!

## ✅ That's It!

Once you complete these steps, your WES Pipeline will have full Firebase Authentication!

## 🐛 Troubleshooting

### Backend shows: "Firebase service account not found"
**Solution**: You need to download the JSON file from Step 1 above

### Frontend shows: "Module not found: firebase/auth"
**Solution**: 
```bash
cd frontend
rm -rf node_modules
npm install
```

### Firebase Console: "Email/Password not enabled"
**Solution**: Follow Step 2 above to enable it

## 📚 More Help

- Quick Guide: [FIREBASE_QUICK_START.md](./FIREBASE_QUICK_START.md)
- Detailed Setup: [FIREBASE_SETUP.md](./FIREBASE_SETUP.md)
- General Setup: [QUICKSTART.md](./QUICKSTART.md)

## 🎯 Summary

```
Backend:  ✅ Running (needs Firebase JSON)
Frontend: ✅ Ready (Firebase SDK installed)
Firebase: ⏳ Needs service account JSON + enable Email/Password

Total time: ~5 minutes
```

**Next**: Complete Steps 1-2 above, restart backend, and you're done!
