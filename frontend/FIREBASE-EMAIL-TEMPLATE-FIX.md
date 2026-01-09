# Firebase Email Template Fix - Remove Space Character

## Problem
The verification email link contains an encoded space (`%20`) causing a 404 error:
```
https://atgcflow.com/__/auth/action%20?mode=verifyEmail...
                                    ^^^ This space causes 404
```

## Root Cause
There is a **space character** in your Firebase email template after the `%LINK%` placeholder.

## Fix Instructions (URGENT - Do This Now)

### Step 1: Open Firebase Console
1. Go to: https://console.firebase.google.com/
2. Select your project: `variant-ac1c6`
3. Navigate to: **Authentication** → **Templates** (in left sidebar)

### Step 2: Edit Email Verification Template
1. Find **Email address verification** template
2. Click the **pencil icon** (Edit) on the right side

### Step 3: Check the Template Body
Look for the verification link in your template. It should look like this:

**❌ INCORRECT (has space - causes 404):**
```html
<a href="%LINK% ">Verify Email</a>
              ^^^^ Space here!
```
OR
```html
Click here: %LINK% ?continue=...
                  ^^^^ Space here!
```

**✅ CORRECT (no space):**
```html
<a href="%LINK%">Verify Email</a>
             ^^^ No space!
```

### Step 4: Remove ALL Spaces Around %LINK%
1. Find `%LINK%` in your template
2. Remove any spaces **immediately after** `%LINK%`
3. Remove any spaces **immediately before** the next character
4. Make sure there's no space between `%LINK%` and the closing quote `"`

### Step 5: Save Template
1. Click **Save** button
2. Wait for confirmation message

### Step 6: Test Verification Flow
1. Register a **new test account** with a fresh email address
2. Check the verification email
3. Click the verification link
4. Confirm it now works without 404 error

---

## Complete Working Email Template

Here's the complete, tested template without any spacing issues:

### Template Configuration:
- **Sender name:** ATGCFLOW
- **Sender email:** noreply@atgcflow.com
- **Reply-to:** aiimsgenomics@gmail.com
- **Subject:** Verify your email for ATGCFLOW

### Template Body (Copy this exactly):

```html
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
</head>
<body style="margin: 0; padding: 0; background-color: #f3f4f6;">
  <table role="presentation" style="width: 100%; border-collapse: collapse;">
    <tr>
      <td align="center" style="padding: 40px 0;">
        <table role="presentation" style="width: 600px; max-width: 100%; border-collapse: collapse; background: white; border-radius: 12px; overflow: hidden; box-shadow: 0 4px 6px rgba(0,0,0,0.1);">

          <!-- Header -->
          <tr>
            <td style="background: linear-gradient(135deg, #2563eb 0%, #06b6d4 100%); padding: 40px; text-align: center;">
              <h1 style="color: white; margin: 0; font-size: 32px; font-family: Arial, sans-serif;">ATGCFLOW</h1>
              <p style="color: rgba(255,255,255,0.95); margin: 8px 0 0 0; font-size: 14px; font-family: Arial, sans-serif;">Whole Exome Sequencing Analysis Platform</p>
            </td>
          </tr>

          <!-- Content -->
          <tr>
            <td style="padding: 40px 30px;">
              <h2 style="color: #1f2937; margin: 0 0 20px 0; font-size: 24px; font-family: Arial, sans-serif;">Hello %DISPLAY_NAME%,</h2>

              <p style="color: #4b5563; line-height: 1.6; font-size: 16px; margin: 0 0 16px 0; font-family: Arial, sans-serif;">
                Thank you for registering with <strong>ATGCFLOW</strong>! We're excited to have you join our platform.
              </p>

              <p style="color: #4b5563; line-height: 1.6; font-size: 16px; margin: 0 0 30px 0; font-family: Arial, sans-serif;">
                To complete your account setup and start analyzing your genomic data, please verify your email address by clicking the button below:
              </p>

              <!-- Button -->
              <table role="presentation" style="width: 100%; border-collapse: collapse; margin: 30px 0;">
                <tr>
                  <td align="center">
                    <a href="%LINK%" style="background: #2563eb; color: white; padding: 16px 48px; text-decoration: none; border-radius: 8px; font-weight: 600; display: inline-block; font-size: 16px; font-family: Arial, sans-serif;">Verify Email Address</a>
                  </td>
                </tr>
              </table>

              <p style="color: #6b7280; font-size: 14px; line-height: 1.6; margin: 25px 0 0 0; font-family: Arial, sans-serif;">
                Or copy and paste this link into your browser:
              </p>
              <p style="color: #2563eb; font-size: 14px; word-break: break-all; margin: 8px 0 0 0; font-family: Arial, sans-serif;">
                %LINK%
              </p>

              <!-- Warning Box -->
              <table role="presentation" style="width: 100%; border-collapse: collapse; margin: 30px 0;">
                <tr>
                  <td style="background: #fef3c7; border-left: 4px solid #f59e0b; padding: 16px; border-radius: 6px;">
                    <p style="color: #92400e; margin: 0; font-size: 14px; font-family: Arial, sans-serif;">
                      <strong>⚠️ Important:</strong> This verification link will expire in 24 hours.
                    </p>
                  </td>
                </tr>
              </table>

              <p style="color: #6b7280; font-size: 14px; line-height: 1.6; margin: 25px 0 0 0; font-family: Arial, sans-serif;">
                If you didn't create an account with ATGCFLOW, you can safely ignore this email.
              </p>
            </td>
          </tr>

          <!-- Footer -->
          <tr>
            <td style="background: #f9fafb; padding: 30px; border-top: 1px solid #e5e7eb;">
              <p style="color: #1f2937; margin: 0 0 20px 0; font-weight: 600; font-family: Arial, sans-serif; font-size: 14px;">
                Best regards,<br/>
                The ATGCFLOW Team
              </p>

              <table role="presentation" style="width: 100%; border-collapse: collapse; border-top: 1px solid #e5e7eb; padding-top: 20px; margin-top: 20px;">
                <tr>
                  <td>
                    <p style="color: #6b7280; font-size: 13px; margin: 5px 0; font-family: Arial, sans-serif;">
                      📧 <a href="mailto:support@atgcflow.com" style="color: #2563eb; text-decoration: none;">support@atgcflow.com</a>
                    </p>
                    <p style="color: #6b7280; font-size: 13px; margin: 5px 0; font-family: Arial, sans-serif;">
                      🌐 <a href="https://atgcflow.com" style="color: #2563eb; text-decoration: none;">https://atgcflow.com</a>
                    </p>
                    <p style="color: #9ca3af; font-size: 12px; margin: 15px 0 0 0; font-style: italic; font-family: Arial, sans-serif;">
                      Industry-grade whole exome sequencing analysis | Research Use Only
                    </p>
                  </td>
                </tr>
              </table>
            </td>
          </tr>

        </table>
      </td>
    </tr>
  </table>
</body>
</html>
```

---

## Critical Points to Check

### ✅ In the template, verify:
1. `<a href="%LINK%"` has **NO space** between `%LINK%` and the closing `"`
2. The plain text link `%LINK%` stands alone with no spaces around it
3. There are no extra spaces or line breaks that could be encoded

### ✅ Action URL Configuration:
- Go to **Authentication** → **Templates** → **Action URL**
- Set to: `https://atgcflow.com/__/auth/action`
- **IMPORTANT:** No spaces, no trailing slash

---

## After Fixing

### Test the Complete Flow:
1. **Register** a new test account (use a fresh email)
2. **Check email** - should arrive from `noreply@atgcflow.com`
3. **Click link** - should go to `https://atgcflow.com/__/auth/action?mode=verifyEmail...`
4. **Verify success** - should see Firebase verification success page
5. **Login** - should be able to access dashboard

### Expected URL (No space!):
```
✅ CORRECT:
https://atgcflow.com/__/auth/action?mode=verifyEmail&oobCode=...
                                   ^^^ No space here!

❌ INCORRECT (what you have now):
https://atgcflow.com/__/auth/action%20?mode=verifyEmail&oobCode=...
                                   ^^^^^ Space encoded as %20
```

---

## Emergency Workaround (Temporary)

If you need to verify an existing email immediately while fixing the template:

1. Copy the verification link from your email
2. Manually remove the `%20` from the URL:
   ```
   Change: https://atgcflow.com/__/auth/action%20?mode=verifyEmail...
   To:     https://atgcflow.com/__/auth/action?mode=verifyEmail...
   ```
3. Paste the corrected URL in your browser
4. This will work, but you must fix the template for all future users

---

## Need Help?

If the issue persists after removing all spaces:

1. **Reset template to default:**
   - Firebase Console → Templates → Email address verification
   - Click "Reset to default template"
   - Save and test again

2. **Check browser network tab:**
   - Open Developer Tools (F12)
   - Click verification link
   - Check the actual URL being requested
   - Look for `%20` or other encoded characters

3. **Verify Next.js rewrites are working:**
   - Visit: `https://atgcflow.com/__/auth/handler`
   - Should get HTTP 200 response (we already verified this works)

---

## Summary

**The fix is simple:** Remove the space character in your Firebase email template after `%LINK%`.

This is a Firebase Console configuration issue, not a code issue. All the code is working correctly. The space in the email template is causing the URL to be malformed, resulting in a 404 error.

**Fix it now in Firebase Console and test with a new registration.**
