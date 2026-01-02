import os
import firebase_admin
from firebase_admin import credentials, auth as firebase_auth
from fastapi import Depends, HTTPException, status
from fastapi.security import HTTPBearer, HTTPAuthorizationCredentials
from sqlalchemy.orm import Session
from database import get_db, User

# Initialize Firebase Admin SDK
SERVICE_ACCOUNT_PATH = os.getenv(
    "GOOGLE_APPLICATION_CREDENTIALS",
    "./firebase-service-account.json"
)

FIREBASE_INITIALIZED = False

try:
    if os.path.exists(SERVICE_ACCOUNT_PATH):
        cred = credentials.Certificate(SERVICE_ACCOUNT_PATH)
        firebase_admin.initialize_app(cred)
        FIREBASE_INITIALIZED = True
        print(f"✅ Firebase initialized with service account: {SERVICE_ACCOUNT_PATH}")
    else:
        print(f"⚠️  Firebase service account not found at: {SERVICE_ACCOUNT_PATH}")
        print(f"⚠️  Download from: https://console.firebase.google.com/project/variant-ac1c6/settings/serviceaccounts/adminsdk")
        print(f"⚠️  Save as: {SERVICE_ACCOUNT_PATH}")
        print(f"⚠️  Authentication will not work until configured!")
except ValueError as e:
    print(f"ℹ️  Firebase app already initialized: {e}")
    FIREBASE_INITIALIZED = True
except Exception as e:
    print(f"❌ Firebase initialization failed: {e}")
    print(f"⚠️  Authentication will not work!")

security = HTTPBearer()

async def get_current_user(
    credentials: HTTPAuthorizationCredentials = Depends(security),
    db: Session = Depends(get_db)
):
    """
    Verify Firebase token and return user info
    """
    if not FIREBASE_INITIALIZED:
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Firebase Authentication not configured. Please download service account JSON.",
        )
    
    token = credentials.credentials
    
    try:
        # Verify the Firebase token
        decoded_token = firebase_auth.verify_id_token(token)
        uid = decoded_token['uid']
        email = decoded_token.get('email', '')
        display_name = decoded_token.get('name', email.split('@')[0] if email else uid[:8])
        
        # Check if user exists in database, if not create them
        user = db.query(User).filter(User.firebase_uid == uid).first()
        
        if not user:
            print(f"Creating new user: {email} (UID: {uid})")
            user = User(
                firebase_uid=uid,
                email=email,
                username=display_name
            )
            db.add(user)
            db.commit()
            db.refresh(user)
        
        return user
        
    except firebase_auth.InvalidIdTokenError:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid authentication token",
            headers={"WWW-Authenticate": "Bearer"},
        )
    except firebase_auth.ExpiredIdTokenError:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Authentication token expired. Please log in again.",
            headers={"WWW-Authenticate": "Bearer"},
        )
    except Exception as e:
        print(f"Auth error: {e}")
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Authentication failed",
            headers={"WWW-Authenticate": "Bearer"},
        )
