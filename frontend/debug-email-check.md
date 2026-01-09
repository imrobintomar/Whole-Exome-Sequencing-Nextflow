# Email Verification Debugging Checklist

## Browser Console Check

When you register a new account, open browser Developer Tools (F12) and check for:

1. **Console Tab** - Look for errors like:
   - `auth/operation-not-allowed`
   - `auth/unauthorized-continue-uri`
   - Network errors

2. **Network Tab** - Look for requests to:
   - `identitytoolkit.googleapis.com`
   - Check if the request succeeds (Status 200)

## Firebase Console Checklist

### ✅ Authentication Setup
- [ ] Go to Firebase Console → Authentication → Sign-in method
- [ ] Email/Password provider is **ENABLED**
- [ ] Status shows "Enabled" with green indicator

### ✅ Email Templates
- [ ] Go to Authentication → Templates
- [ ] Email address verification template exists
- [ ] Click "Edit template" and verify:
  - [ ] Sender name: (e.g., "ATGCFLOW" or your app name)
  - [ ] Sender email: `noreply@variant-ac1c6.firebaseapp.com`
  - [ ] Subject: Has a subject line
  - [ ] Body: Contains %LINK% placeholder
  - [ ] Action URL: `https://atgcflow.com/__/auth/action`

### ✅ Authorized Domains
- [ ] Go to Authentication → Settings → Authorized domains
- [ ] These domains are listed:
  - [ ] `atgcflow.com`
  - [ ] `localhost`
  - [ ] `variant-ac1c6.firebaseapp.com` (should be there by default)

## Testing Steps

### Step 1: Test with Console Output
Add this to your RegisterForm to see detailed logs:

```typescript
try {
  console.log('🔵 Starting registration...');
  const userCredential = await createUserWithEmailAndPassword(auth, email, password);
  console.log('✅ User created:', userCredential.user.uid);

  console.log('📧 Sending verification email...');
  await sendEmailVerification(userCredential.user, {
    url: continueUrl,
    handleCodeInApp: false
  });
  console.log('✅ Verification email sent!');

} catch (err) {
  console.error('❌ Error:', err.code, err.message);
}
```

### Step 2: Check Firebase Logs
- Go to Firebase Console → Project Overview → Usage
- Check if email sending quota is exhausted (unlikely for new projects)

### Step 3: Test with Different Email
- Try with Gmail, Outlook, Yahoo
- Some email providers have aggressive spam filters

## Quick Fix Attempts

### Option 1: Reset Firebase Auth Settings
1. Disable Email/Password provider
2. Re-enable it
3. Try again

### Option 2: Use Default Email Template
1. Go to Templates → Email address verification
2. Click "Reset to default template"
3. Save and try again

### Option 3: Check Email Sending Status
1. Firebase Console → Project Settings
2. Look for any warnings about email sending

## Expected Behavior

When working correctly:
1. User registers
2. Console shows: "Verification email sent!"
3. Email arrives within 1-2 minutes
4. Email comes from: `noreply@variant-ac1c6.firebaseapp.com`
5. Subject: "Verify your email for variant-ac1c6"

## If Still Not Working

Contact Firebase Support or check:
- Firebase Status Page: https://status.firebase.google.com/
- Make sure your Firebase project isn't in sandbox mode
- Verify billing is enabled (free tier should work, but check)
