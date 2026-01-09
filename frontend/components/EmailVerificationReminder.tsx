'use client';

import { useState } from 'react';
import { sendEmailVerification, signOut } from 'firebase/auth';
import { auth } from '@/lib/firebase';
import { Card, CardContent, CardHeader, CardTitle } from './ui/card';
import { Button } from './ui/button';
import { Mail, AlertCircle, CheckCircle, Loader } from 'lucide-react';

interface EmailVerificationReminderProps {
  email: string;
}

export default function EmailVerificationReminder({ email }: EmailVerificationReminderProps) {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [success, setSuccess] = useState(false);

  const handleResendVerification = async () => {
    if (!auth.currentUser) return;

    setLoading(true);
    setError('');
    setSuccess(false);

    try {
      const continueUrl = `${window.location.origin}/?verified=true`;
      await sendEmailVerification(auth.currentUser, {
        url: continueUrl,
        handleCodeInApp: false
      });
      setSuccess(true);
    } catch (err: any) {
      console.error('Resend verification error:', err);
      if (err.code === 'auth/too-many-requests') {
        setError('Too many requests. Please wait a few minutes before trying again.');
      } else {
        setError('Failed to send verification email. Please try again later.');
      }
    } finally {
      setLoading(false);
    }
  };

  const handleSignOut = async () => {
    await signOut(auth);
  };

  return (
    <div className="min-h-screen bg-gradient-to-br from-amber-50 via-orange-50 to-red-50 flex items-center justify-center p-4">
      <Card className="w-full max-w-2xl">
        <CardHeader className="space-y-1 text-center">
          <div className="mx-auto mb-2 flex h-16 w-16 items-center justify-center rounded-full bg-amber-100">
            <AlertCircle className="h-10 w-10 text-amber-600" />
          </div>
          <CardTitle className="text-2xl font-bold">Email Verification Required</CardTitle>
        </CardHeader>
        <CardContent className="space-y-6">
          <div className="rounded-lg bg-amber-50 border-2 border-amber-300 p-6">
            <div className="flex items-start gap-3">
              <Mail className="h-6 w-6 text-amber-600 mt-0.5 flex-shrink-0" />
              <div className="flex-1">
                <h3 className="text-lg font-semibold text-amber-900 mb-2">
                  Please verify your email address
                </h3>
                <p className="text-sm text-amber-800 mb-4">
                  We sent a verification email to <strong className="font-semibold">{email}</strong>
                </p>
                <div className="bg-white rounded-md p-4 mb-4">
                  <p className="text-sm font-medium text-slate-900 mb-2">To verify your account:</p>
                  <ol className="text-sm text-slate-700 space-y-2 list-decimal list-inside">
                    <li>Check your email inbox (and spam folder)</li>
                    <li>Click the verification link in the email</li>
                    <li>Return to this page and refresh</li>
                  </ol>
                </div>
                <p className="text-xs text-amber-700">
                  The verification link expires in 24 hours. If you haven't received the email after a few minutes, you can request a new one below.
                </p>
              </div>
            </div>
          </div>

          {/* Success message */}
          {success && (
            <div className="rounded-lg bg-green-50 border border-green-300 p-4">
              <div className="flex items-start gap-3">
                <CheckCircle className="h-5 w-5 text-green-600 mt-0.5 flex-shrink-0" />
                <div>
                  <p className="text-sm font-medium text-green-900">Verification email sent!</p>
                  <p className="text-sm text-green-700 mt-1">
                    Please check your inbox and click the verification link.
                  </p>
                </div>
              </div>
            </div>
          )}

          {/* Error message */}
          {error && (
            <div className="rounded-lg bg-red-50 border border-red-300 p-4">
              <div className="flex items-start gap-3">
                <AlertCircle className="h-5 w-5 text-red-600 mt-0.5 flex-shrink-0" />
                <div>
                  <p className="text-sm font-medium text-red-900">Error</p>
                  <p className="text-sm text-red-700 mt-1">{error}</p>
                </div>
              </div>
            </div>
          )}

          {/* Warning about access */}
          <div className="rounded-lg bg-blue-50 border border-blue-300 p-4">
            <div className="flex items-start gap-3">
              <AlertCircle className="h-5 w-5 text-blue-600 mt-0.5 flex-shrink-0" />
              <div>
                <p className="text-sm font-medium text-blue-900">Access Restricted</p>
                <p className="text-sm text-blue-700 mt-1">
                  You cannot access the ATGCFLOW dashboard until your email is verified. This is a security measure to ensure the integrity of our platform.
                </p>
              </div>
            </div>
          </div>

          {/* Action buttons */}
          <div className="flex flex-col sm:flex-row gap-3">
            <Button
              onClick={handleResendVerification}
              disabled={loading || success}
              className="flex-1"
              variant="default"
            >
              {loading ? (
                <>
                  <Loader className="mr-2 h-4 w-4 animate-spin" />
                  Sending...
                </>
              ) : success ? (
                'Email Sent'
              ) : (
                'Resend Verification Email'
              )}
            </Button>

            <Button
              onClick={handleSignOut}
              variant="outline"
              className="flex-1"
            >
              Sign Out
            </Button>
          </div>

          {/* Additional help */}
          <div className="text-center pt-4 border-t border-slate-200">
            <p className="text-sm text-slate-600 mb-2">
              Need help? Contact us at{' '}
              <a href="mailto:support@atgcflow.com" className="text-blue-600 hover:underline font-medium">
                support@atgcflow.com
              </a>
            </p>
          </div>
        </CardContent>
      </Card>
    </div>
  );
}
