# Firebase Authentication Setup Guide

## Overview

This WES Pipeline now uses Firebase Authentication for user sign-in and sign-up instead of JWT tokens.

## What Changed

### Frontend
- ✅ Firebase SDK installed (`firebase` package)
- ✅ Firebase configuration added ([lib/firebase.ts](frontend/lib/firebase.ts))
- ✅ Login/Register forms updated to use Firebase Auth
- ✅ API client updated to use Firebase ID tokens
- ✅ User state managed with Firebase `onAuthStateChanged`

### Backend
- ✅ Firebase Admin SDK added for token verification
- ✅ Database schema updated (added `firebase_uid` field)
- ✅ JWT auth removed, replaced with Firebase token verification
- ✅ User creation handled automatically on first login

## Backend Setup Required

### Step 1: Get Firebase Service Account Key

1. Go to [Firebase Console](https://console.firebase.google.com/)
2. Select your project: `variant-ac1c6`
3. Click the gear icon → **Project settings**
4. Go to **Service accounts** tab
5. Click **Generate new private key**
6. Download the JSON file

### Step 2: Configure Backend

Create `.env` file in `backend/` directory (if not exists):

```bash
cd backend
cp .env.example .env
```

Add Firebase credentials to `.env`:

```bash
# Existing settings...
SECRET_KEY=your-secret-key-here
UPLOAD_DIR=/media/drprabudh/m3/Nextflow-Script/WholeExome/uploads
RESULTS_DIR=/media/drprabudh/m3/Nextflow-Script/WholeExome/results

# Firebase Admin SDK
FIREBASE_PROJECT_ID=variant-ac1c6
FIREBASE_PRIVATE_KEY_ID=your-private-key-id-from-json
FIREBASE_PRIVATE_KEY="-----BEGIN PRIVATE KEY-----\nYour private key here\n-----END PRIVATE KEY-----\n"
FIREBASE_CLIENT_EMAIL=firebase-adminsdk-xxxxx@variant-ac1c6.iam.gserviceaccount.com
FIREBASE_CLIENT_ID=your-client-id

# Or simpler: just point to the JSON file
GOOGLE_APPLICATION_CREDENTIALS=/path/to/your/firebase-service-account.json
```

### Step 3: Update `firebase_auth.py`

Option A: Use environment variables (recommended for production):

```python
# In firebase_auth.py, replace the cred initialization:
import os

if os.getenv("GOOGLE_APPLICATION_CREDENTIALS"):
    cred = credentials.Certificate(os.getenv("GOOGLE_APPLICATION_CREDENTIALS"))
else:
    cred = credentials.Certificate({
        "type": "service_account",
        "project_id": os.getenv("FIREBASE_PROJECT_ID"),
        "private_key_id": os.getenv("FIREBASE_PRIVATE_KEY_ID"),
        "private_key": os.getenv("FIREBASE_PRIVATE_KEY").replace('\\n', '\n'),
        "client_email": os.getenv("FIREBASE_CLIENT_EMAIL"),
        "client_id": os.getenv("FIREBASE_CLIENT_ID"),
        # ... rest of config
    })
```

Option B: Use JSON file directly (easier for testing):

```python
cred = credentials.Certificate("/path/to/firebase-service-account.json")
```

### Step 4: Install Dependencies

```bash
cd backend
source venv/bin/activate
pip install -r requirements.txt
```

### Step 5: Recreate Database (if needed)

Since we added `firebase_uid` field, you may need to recreate the database:

```bash
# Backup existing database (if any)
cp wes_pipeline.db wes_pipeline.db.backup

# Remove old database
rm wes_pipeline.db

# Database will be recreated automatically on next run
```

### Step 6: Start Backend

```bash
./run.sh
```

## Frontend Setup

### Step 1: Install Dependencies

```bash
cd frontend
npm install
```

### Step 2: Configure API URL

Create `.env.local`:

```bash
NEXT_PUBLIC_API_URL=http://localhost:8000
```

### Step 3: Start Frontend

```bash
npm run dev
```

Visit `http://localhost:3000`

## Firebase Console Setup

### Enable Authentication Methods

1. Go to [Firebase Console](https://console.firebase.google.com/)
2. Select `variant-ac1c6` project
3. Click **Authentication** in left sidebar
4. Go to **Sign-in method** tab
5. Enable **Email/Password**

### Optional: Add Authorized Domains

1. In Authentication → Settings
2. Add authorized domains:
   - `localhost` (for development)
   - Your Vercel domain (for production)
   - Your custom domain (if any)

## How It Works

### User Registration Flow

1. User fills registration form on frontend
2. Frontend calls Firebase `createUserWithEmailAndPassword()`
3. Firebase creates user account
4. User is automatically logged in
5. Frontend gets Firebase ID token
6. Backend receives ID token, verifies it, creates User record

### User Login Flow

1. User fills login form
2. Frontend calls Firebase `signInWithEmailAndPassword()`
3. Firebase authenticates user
4. Frontend gets Firebase ID token
5. All API requests include: `Authorization: Bearer {firebase-id-token}`
6. Backend verifies token with Firebase Admin SDK
7. Backend retrieves/creates user from database

### API Request Flow

```
Frontend → API Request with Firebase ID Token
  ↓
Backend → Verify token with Firebase Admin SDK
  ↓
Backend → Get user from database (create if not exists)
  ↓
Backend → Process request
  ↓
Backend → Return response
```

## Benefits of Firebase Auth

1. **Security**: Industry-standard authentication
2. **Easy Setup**: No custom JWT implementation
3. **Features**: Built-in password reset, email verification
4. **Scalability**: Handles millions of users
5. **Multi-platform**: Same auth works for web, iOS, Android
6. **Free Tier**: 50,000 MAU (Monthly Active Users) free

## Testing

1. **Register**: Create a new account with email/password
2. **Login**: Sign in with created credentials
3. **Upload**: Submit FASTQ files - should work
4. **Logout**: Sign out
5. **Login Again**: Sign back in - jobs should persist

## Troubleshooting

### "Invalid authentication credentials"
- Check Firebase service account JSON is correct
- Verify `GOOGLE_APPLICATION_CREDENTIALS` path
- Check backend logs for detailed error

### "Email already in use"
- User already registered in Firebase
- Use password reset or different email

### "Token expired"
- Firebase tokens expire after 1 hour
- Frontend automatically refreshes tokens
- Re-login if issues persist

### Backend won't start
```bash
# Check dependencies
pip list | grep firebase

# Should show:
# firebase-admin==6.6.0

# Reinstall if missing
pip install firebase-admin==6.6.0
```

## Migration from JWT to Firebase

If you had users with the old JWT system:

1. Users must re-register with Firebase
2. Old database will be recreated with new schema
3. Previous user data will be lost (unless manually migrated)

For production, consider writing a migration script to:
- Export old user emails
- Create Firebase accounts
- Map old user IDs to new Firebase UIDs

## Production Deployment

### Security Best Practices

1. **Never commit** Firebase service account JSON to git
2. Use environment variables for credentials
3. Enable Firebase App Check for additional security
4. Set up Firebase Security Rules
5. Enable Firebase Audit Logging

### Vercel Deployment

Add environment variables in Vercel:
```
NEXT_PUBLIC_API_URL=https://your-backend-url
```

Firebase config is already in the code (it's okay to commit apiKey, it's public).

### Backend Deployment

Set environment variables on your server:
```bash
export GOOGLE_APPLICATION_CREDENTIALS=/path/to/service-account.json
# OR
export FIREBASE_PROJECT_ID=variant-ac1c6
export FIREBASE_PRIVATE_KEY="..."
# etc.
```

## Support

For Firebase issues:
- [Firebase Auth Documentation](https://firebase.google.com/docs/auth)
- [Firebase Admin SDK for Python](https://firebase.google.com/docs/admin/setup)

For Pipeline issues:
- Check `backend/` logs
- Review `QUICKSTART.md` for general setup
