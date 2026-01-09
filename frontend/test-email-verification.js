// Test script to verify Firebase email verification is working
// Run with: node test-email-verification.js

const { initializeApp } = require('firebase/app');
const { getAuth, createUserWithEmailAndPassword, sendEmailVerification } = require('firebase/auth');

const firebaseConfig = {
  apiKey: "AIzaSyAcMhwJRa0diLwHw2bZaCLIr_akFZPAsWc",
  authDomain: "variant-ac1c6.firebaseapp.com",
  projectId: "variant-ac1c6",
  storageBucket: "variant-ac1c6.firebasestorage.app",
  messagingSenderId: "290831279354",
  appId: "1:290831279354:web:d1a3bb9be9a8345a455543"
};

// Initialize Firebase
const app = initializeApp(firebaseConfig);
const auth = getAuth(app);

async function testEmailVerification() {
  const testEmail = `test+${Date.now()}@example.com`;
  const testPassword = 'TestPassword123!';

  console.log('🧪 Testing Firebase Email Verification...\n');
  console.log(`📧 Test Email: ${testEmail}`);
  console.log(`🔑 Test Password: ${testPassword}\n`);

  try {
    // Create test user
    console.log('1️⃣ Creating test user...');
    const userCredential = await createUserWithEmailAndPassword(auth, testEmail, testPassword);
    console.log('✅ User created successfully');
    console.log(`   UID: ${userCredential.user.uid}`);
    console.log(`   Email: ${userCredential.user.email}`);
    console.log(`   Email Verified: ${userCredential.user.emailVerified}\n`);

    // Send verification email
    console.log('2️⃣ Sending verification email...');
    await sendEmailVerification(userCredential.user, {
      url: 'https://atgcflow.com/?verified=true',
      handleCodeInApp: false
    });
    console.log('✅ Verification email sent successfully!\n');

    console.log('📬 Check the inbox for:', testEmail);
    console.log('   (Note: This is a test email, check Firebase Console for actual email settings)\n');

    console.log('🎉 Test completed successfully!');
    console.log('\nIf you don\'t receive emails in production:');
    console.log('1. Check Firebase Console → Authentication → Templates');
    console.log('2. Verify sender email and template are properly configured');
    console.log('3. Check spam folder');
    console.log('4. Ensure domain is authorized in Firebase');

    process.exit(0);
  } catch (error) {
    console.error('❌ Test failed!');
    console.error('Error code:', error.code);
    console.error('Error message:', error.message);
    console.error('\nCommon issues:');
    console.error('- Email/Password sign-in not enabled in Firebase Console');
    console.error('- Authorized domains not configured');
    console.error('- Email template not set up');
    process.exit(1);
  }
}

testEmailVerification();
