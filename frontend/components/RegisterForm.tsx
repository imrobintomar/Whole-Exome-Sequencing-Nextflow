'use client';

import { useState } from 'react';
import { createUserWithEmailAndPassword, updateProfile, sendEmailVerification } from 'firebase/auth';
import { auth } from '@/lib/firebase';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from './ui/card';
import { Input } from './ui/input';
import { Label } from './ui/label';
import { Button } from './ui/button';
import { Dna, Loader, Mail, CheckCircle } from 'lucide-react';

interface RegisterFormProps {
  onRegister: () => void;
  onToggle: () => void;
}

export default function RegisterForm({ onRegister, onToggle }: RegisterFormProps) {
  const [email, setEmail] = useState('');
  const [displayName, setDisplayName] = useState('');
  const [password, setPassword] = useState('');
  const [confirmPassword, setConfirmPassword] = useState('');
  const [error, setError] = useState('');
  const [loading, setLoading] = useState(false);
  const [verificationSent, setVerificationSent] = useState(false);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setError('');

    if (password !== confirmPassword) {
      setError('Passwords do not match');
      return;
    }

    if (password.length < 6) {
      setError('Password must be at least 6 characters');
      return;
    }

    setLoading(true);

    try {
      console.log('🔵 [Registration] Starting user registration...');
      console.log('📧 [Registration] Email:', email);

      const userCredential = await createUserWithEmailAndPassword(auth, email, password);
      console.log('✅ [Registration] User created successfully!');
      console.log('👤 [Registration] UID:', userCredential.user.uid);
      console.log('📧 [Registration] Email:', userCredential.user.email);
      console.log('✉️ [Registration] Email Verified:', userCredential.user.emailVerified);

      // Update display name
      if (displayName) {
        console.log('📝 [Registration] Updating display name...');
        await updateProfile(userCredential.user, {
          displayName: displayName
        });
        console.log('✅ [Registration] Display name updated');
      }

      // Send email verification with continue URL
      const continueUrl = `${window.location.origin}/?verified=true`;
      console.log('📬 [Registration] Sending verification email...');
      console.log('🔗 [Registration] Continue URL:', continueUrl);

      await sendEmailVerification(userCredential.user, {
        url: continueUrl,
        handleCodeInApp: false
      });

      console.log('✅ [Registration] Verification email sent successfully!');
      console.log('📮 [Registration] Check inbox for:', email);
      console.log('💡 [Registration] If not received, check spam folder');

      setVerificationSent(true);

      // Sign out the user immediately after registration so they can't access the dashboard without verification
      console.log('🔐 [Registration] Signing out user...');
      await auth.signOut();
      console.log('✅ [Registration] User signed out');
    } catch (err: any) {
      console.error('❌ [Registration] Error:', err);
      console.error('❌ [Registration] Error code:', err.code);
      console.error('❌ [Registration] Error message:', err.message);
      const errorMessage = err.code === 'auth/email-already-in-use'
        ? 'Email already registered'
        : err.code === 'auth/weak-password'
        ? 'Password is too weak'
        : err.code === 'auth/invalid-email'
        ? 'Invalid email address'
        : 'Registration failed. Please try again.';
      setError(errorMessage);
    } finally {
      setLoading(false);
    }
  };

  // If verification email was sent, show success message
  if (verificationSent) {
    return (
      <div className="flex min-h-screen items-center justify-center bg-gradient-to-br from-blue-50 via-indigo-50 to-purple-50 p-4">
        <Card className="w-full max-w-md">
          <CardHeader className="space-y-1 text-center">
            <div className="mx-auto mb-2 flex h-16 w-16 items-center justify-center rounded-full bg-green-100">
              <CheckCircle className="h-10 w-10 text-green-600" />
            </div>
            <CardTitle className="text-2xl font-bold">Verify Your Email</CardTitle>
            <CardDescription>
              We've sent a verification link to your email
            </CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            <div className="rounded-lg bg-blue-50 border border-blue-200 p-4">
              <div className="flex items-start gap-3">
                <Mail className="h-5 w-5 text-blue-600 mt-0.5 flex-shrink-0" />
                <div className="flex-1">
                  <p className="text-sm font-medium text-blue-900 mb-2">Check your inbox</p>
                  <p className="text-sm text-blue-800 mb-3">
                    We sent a verification email to <strong>{email}</strong>
                  </p>
                  <ul className="text-sm text-blue-700 space-y-1 list-disc list-inside">
                    <li>Click the verification link in the email</li>
                    <li>Return to this page and sign in</li>
                    <li>Check your spam folder if you don't see it</li>
                  </ul>
                </div>
              </div>
            </div>

            <div className="bg-amber-50 border border-amber-200 rounded-lg p-4">
              <p className="text-sm text-amber-900">
                <strong>Note:</strong> You must verify your email before you can access the dashboard. The verification link expires in 24 hours.
              </p>
            </div>

            <Button
              onClick={onToggle}
              className="w-full"
              variant="default"
            >
              Go to Sign In
            </Button>
          </CardContent>
        </Card>
      </div>
    );
  }

  return (
    <div className="flex min-h-screen items-center justify-center bg-gradient-to-br from-blue-50 via-indigo-50 to-purple-50 p-4">
      <Card className="w-full max-w-md">
        <CardHeader className="space-y-1 text-center">
          <div className="mx-auto mb-2 flex h-12 w-12 items-center justify-center rounded-full bg-primary/10">
            <Dna className="h-6 w-6 text-primary" />
          </div>
          <CardTitle className="text-2xl font-bold">Create Account</CardTitle>
          <CardDescription>
            Join ATGCFLOW to start analyzing your Whole Exome data
          </CardDescription>
        </CardHeader>
        <CardContent>
          <form onSubmit={handleSubmit} className="space-y-4">
            <div className="space-y-2">
              <Label htmlFor="email">Email</Label>
              <Input
                id="email"
                type="email"
                value={email}
                onChange={(e) => setEmail(e.target.value)}
                placeholder="name@example.com"
                required
                disabled={loading}
              />
            </div>

            <div className="space-y-2">
              <Label htmlFor="displayName">Display Name</Label>
              <Input
                id="displayName"
                type="text"
                value={displayName}
                onChange={(e) => setDisplayName(e.target.value)}
                placeholder="Optional"
                disabled={loading}
              />
            </div>

            <div className="space-y-2">
              <Label htmlFor="password">Password</Label>
              <Input
                id="password"
                type="password"
                value={password}
                onChange={(e) => setPassword(e.target.value)}
                placeholder="At least 6 characters"
                required
                disabled={loading}
              />
            </div>

            <div className="space-y-2">
              <Label htmlFor="confirmPassword">Confirm Password</Label>
              <Input
                id="confirmPassword"
                type="password"
                value={confirmPassword}
                onChange={(e) => setConfirmPassword(e.target.value)}
                placeholder="Re-enter your password"
                required
                disabled={loading}
              />
            </div>

            {error && (
              <Card className="border-destructive/50 bg-destructive/5">
                <CardContent className="pt-6 pb-4">
                  <p className="text-sm text-destructive">{error}</p>
                </CardContent>
              </Card>
            )}

            <Button
              type="submit"
              disabled={loading}
              className="w-full"
            >
              {loading ? (
                <>
                  <Loader className="mr-2 h-4 w-4 animate-spin" />
                  Creating account...
                </>
              ) : (
                'Create Account'
              )}
            </Button>
          </form>

          <div className="mt-6 text-center text-sm">
            <span className="text-muted-foreground">Already have an account? </span>
            <button
              type="button"
              onClick={onToggle}
              className="font-semibold text-primary hover:underline"
            >
              Sign in
            </button>
          </div>
        </CardContent>
      </Card>
    </div>
  );
}
