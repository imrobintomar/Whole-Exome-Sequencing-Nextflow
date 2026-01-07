'use client';

import { useState, useEffect } from 'react';
import { billingApi, UserSubscription, UsageStats } from '../lib/api';

export default function BillingPage() {
  const [subscription, setSubscription] = useState<UserSubscription | null>(null);
  const [usage, setUsage] = useState<UsageStats | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [actionLoading, setActionLoading] = useState(false);

  useEffect(() => {
    loadBillingData();
  }, []);

  const loadBillingData = async () => {
    try {
      const [subData, usageData] = await Promise.all([
        billingApi.getSubscription(),
        billingApi.getUsage(),
      ]);

      setSubscription(subData);
      setUsage(usageData);
      setError(null);
    } catch (err) {
      setError('Failed to load billing information');
      console.error('Billing error:', err);
    } finally {
      setLoading(false);
    }
  };

  const handleManageSubscription = async () => {
    setActionLoading(true);
    try {
      const { portal_url } = await billingApi.createPortalSession();
      window.location.href = portal_url;
    } catch (err) {
      alert('Failed to open customer portal. Please try again.');
      setActionLoading(false);
    }
  };

  if (loading) {
    return (
      <div className="flex items-center justify-center min-h-screen">
        <div className="text-center">
          <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-600 mx-auto mb-4"></div>
          <p className="text-gray-600">Loading billing information...</p>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="flex items-center justify-center min-h-screen">
        <div className="bg-red-50 border border-red-200 rounded-lg p-6 max-w-md">
          <h3 className="text-red-800 font-semibold mb-2">Error</h3>
          <p className="text-red-600">{error}</p>
          <button
            onClick={loadBillingData}
            className="mt-4 px-4 py-2 bg-red-600 text-white rounded hover:bg-red-700"
          >
            Retry
          </button>
        </div>
      </div>
    );
  }

  if (!subscription || !usage) return null;

  const isFreePlan = subscription.plan.name === 'Free';
  const hasActiveSubscription = subscription.subscription && subscription.subscription.status === 'active';
  const usagePercentage = (usage.jobs_executed / usage.jobs_limit) * 100;

  return (
    <div className="min-h-screen bg-gray-50">
      {/* Header */}
      <div className="bg-white border-b">
        <div className="max-w-4xl mx-auto px-4 py-6">
          <h1 className="text-3xl font-bold text-gray-900">Billing & Subscription</h1>
          <p className="text-gray-600 mt-1">Manage your subscription and view usage</p>
        </div>
      </div>

      <div className="max-w-4xl mx-auto px-4 py-8 space-y-6">
        {/* Current Plan Card */}
        <div className="bg-white rounded-lg shadow">
          <div className="px-6 py-4 border-b border-gray-200">
            <h2 className="text-lg font-semibold">Current Plan</h2>
          </div>
          <div className="p-6">
            <div className="flex items-start justify-between">
              <div>
                <div className="flex items-center gap-3 mb-2">
                  <h3 className="text-2xl font-bold text-gray-900">{subscription.plan.name}</h3>
                  {hasActiveSubscription && (
                    <span className="px-3 py-1 bg-green-100 text-green-800 text-sm font-medium rounded-full">
                      Active
                    </span>
                  )}
                  {subscription.subscription?.cancel_at_period_end && (
                    <span className="px-3 py-1 bg-orange-100 text-orange-800 text-sm font-medium rounded-full">
                      Cancels {new Date(subscription.subscription.current_period_end).toLocaleDateString()}
                    </span>
                  )}
                </div>
                <p className="text-3xl font-bold text-gray-900 mb-1">
                  ${(subscription.plan.price_cents / 100).toFixed(2)}
                  <span className="text-lg text-gray-600 font-normal">/month</span>
                </p>
                <p className="text-gray-600">
                  {subscription.plan.monthly_jobs_limit} WES jobs per month
                </p>
                {subscription.plan.chat_support && (
                  <div className="mt-2 flex items-center text-sm text-gray-600">
                    <svg className="h-4 w-4 text-green-500 mr-1" fill="currentColor" viewBox="0 0 20 20">
                      <path
                        fillRule="evenodd"
                        d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z"
                        clipRule="evenodd"
                      />
                    </svg>
                    Live chat support included
                  </div>
                )}
              </div>

              {!isFreePlan && hasActiveSubscription && (
                <button
                  onClick={handleManageSubscription}
                  disabled={actionLoading}
                  className="px-4 py-2 bg-gray-100 text-gray-700 rounded-lg hover:bg-gray-200 font-medium transition-colors disabled:opacity-50"
                >
                  {actionLoading ? 'Loading...' : 'Manage Subscription'}
                </button>
              )}

              {isFreePlan && (
                <a
                  href="/pricing"
                  className="px-4 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700 font-medium transition-colors"
                >
                  Upgrade Plan
                </a>
              )}
            </div>

            {hasActiveSubscription && subscription.subscription && (
              <div className="mt-6 pt-6 border-t border-gray-200">
                <div className="grid grid-cols-2 gap-4 text-sm">
                  <div>
                    <span className="text-gray-600">Current Period:</span>
                    <p className="font-medium text-gray-900">
                      {new Date(subscription.subscription.current_period_start).toLocaleDateString()} -{' '}
                      {new Date(subscription.subscription.current_period_end).toLocaleDateString()}
                    </p>
                  </div>
                  <div>
                    <span className="text-gray-600">Next Billing Date:</span>
                    <p className="font-medium text-gray-900">
                      {new Date(subscription.subscription.current_period_end).toLocaleDateString()}
                    </p>
                  </div>
                </div>
              </div>
            )}

            {subscription.plan.features && subscription.plan.features.length > 0 && (
              <div className="mt-6 pt-6 border-t border-gray-200">
                <h4 className="text-sm font-semibold text-gray-900 mb-3">Plan Features:</h4>
                <ul className="space-y-2">
                  {subscription.plan.features.map((feature, index) => (
                    <li key={index} className="flex items-start text-sm text-gray-600">
                      <svg className="h-5 w-5 text-green-500 mr-2 flex-shrink-0" fill="currentColor" viewBox="0 0 20 20">
                        <path
                          fillRule="evenodd"
                          d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z"
                          clipRule="evenodd"
                        />
                      </svg>
                      {feature}
                    </li>
                  ))}
                </ul>
              </div>
            )}
          </div>
        </div>

        {/* Usage Card */}
        <div className="bg-white rounded-lg shadow">
          <div className="px-6 py-4 border-b border-gray-200">
            <h2 className="text-lg font-semibold">Current Usage</h2>
          </div>
          <div className="p-6">
            <div className="mb-4">
              <div className="flex items-center justify-between mb-2">
                <span className="text-sm text-gray-600">Jobs Used This Month</span>
                <span className="text-sm font-semibold text-gray-900">
                  {usage.jobs_executed} / {usage.jobs_limit}
                </span>
              </div>
              <div className="w-full bg-gray-200 rounded-full h-3">
                <div
                  className={`h-3 rounded-full transition-all ${
                    usagePercentage >= 100 ? 'bg-red-500' : usagePercentage >= 80 ? 'bg-yellow-500' : 'bg-green-500'
                  }`}
                  style={{ width: `${Math.min(usagePercentage, 100)}%` }}
                ></div>
              </div>
              <p className="text-xs text-gray-500 mt-2">
                {usage.jobs_remaining} jobs remaining • {usage.usage_percent.toFixed(0)}% used
              </p>
            </div>

            {usagePercentage >= 100 && (
              <div className="bg-red-50 border border-red-200 rounded-lg p-4">
                <div className="flex items-start">
                  <svg className="h-5 w-5 text-red-600 mr-2 flex-shrink-0 mt-0.5" fill="currentColor" viewBox="0 0 20 20">
                    <path
                      fillRule="evenodd"
                      d="M10 18a8 8 0 100-16 8 8 0 000 16zM8.707 7.293a1 1 0 00-1.414 1.414L8.586 10l-1.293 1.293a1 1 0 101.414 1.414L10 11.414l1.293 1.293a1 1 0 001.414-1.414L11.414 10l1.293-1.293a1 1 0 00-1.414-1.414L10 8.586 8.707 7.293z"
                      clipRule="evenodd"
                    />
                  </svg>
                  <div>
                    <p className="text-sm font-medium text-red-800">Monthly limit reached</p>
                    <p className="text-sm text-red-700 mt-1">
                      You've used all your jobs for this month. Upgrade your plan to submit more jobs.
                    </p>
                  </div>
                </div>
              </div>
            )}

            {usagePercentage >= 80 && usagePercentage < 100 && (
              <div className="bg-yellow-50 border border-yellow-200 rounded-lg p-4">
                <div className="flex items-start">
                  <svg className="h-5 w-5 text-yellow-600 mr-2 flex-shrink-0 mt-0.5" fill="currentColor" viewBox="0 0 20 20">
                    <path
                      fillRule="evenodd"
                      d="M8.257 3.099c.765-1.36 2.722-1.36 3.486 0l5.58 9.92c.75 1.334-.213 2.98-1.742 2.98H4.42c-1.53 0-2.493-1.646-1.743-2.98l5.58-9.92zM11 13a1 1 0 11-2 0 1 1 0 012 0zm-1-8a1 1 0 00-1 1v3a1 1 0 002 0V6a1 1 0 00-1-1z"
                      clipRule="evenodd"
                    />
                  </svg>
                  <div>
                    <p className="text-sm font-medium text-yellow-800">Approaching limit</p>
                    <p className="text-sm text-yellow-700 mt-1">
                      You're using {usage.usage_percent.toFixed(0)}% of your monthly allocation. Consider upgrading if you need more jobs.
                    </p>
                  </div>
                </div>
              </div>
            )}

            <div className="mt-6 pt-6 border-t border-gray-200 grid grid-cols-2 gap-4 text-sm">
              <div>
                <span className="text-gray-600">Plan:</span>
                <p className="font-medium text-gray-900">{usage.plan_name}</p>
              </div>
              <div>
                <span className="text-gray-600">Billing Cycle:</span>
                <p className="font-medium text-gray-900">
                  {new Date(usage.month.toString().substring(0, 4) + '-' + usage.month.toString().substring(4, 6) + '-01').toLocaleString(
                    'default',
                    { month: 'long', year: 'numeric' }
                  )}
                </p>
              </div>
            </div>
          </div>
        </div>

        {/* Upgrade Options (for free tier) */}
        {isFreePlan && (
          <div className="bg-gradient-to-r from-blue-600 to-purple-600 rounded-lg shadow-lg p-8 text-white">
            <h3 className="text-2xl font-bold mb-3">Unlock More Power</h3>
            <p className="text-blue-100 mb-6">
              Upgrade to a paid plan for more monthly jobs, priority processing, and live chat support.
            </p>
            <div className="grid md:grid-cols-2 gap-4">
              <div className="bg-white/10 backdrop-blur rounded-lg p-4">
                <h4 className="font-semibold mb-2">Basic Plan</h4>
                <p className="text-3xl font-bold mb-1">$29<span className="text-lg font-normal">/mo</span></p>
                <p className="text-blue-100 text-sm mb-3">10 jobs per month</p>
                <a
                  href="/pricing"
                  className="block text-center bg-white text-blue-600 px-4 py-2 rounded-lg font-medium hover:bg-blue-50 transition-colors"
                >
                  Learn More
                </a>
              </div>
              <div className="bg-white/10 backdrop-blur rounded-lg p-4">
                <h4 className="font-semibold mb-2">Pro Plan</h4>
                <p className="text-3xl font-bold mb-1">$99<span className="text-lg font-normal">/mo</span></p>
                <p className="text-blue-100 text-sm mb-3">50 jobs per month</p>
                <a
                  href="/pricing"
                  className="block text-center bg-white text-purple-600 px-4 py-2 rounded-lg font-medium hover:bg-purple-50 transition-colors"
                >
                  Learn More
                </a>
              </div>
            </div>
          </div>
        )}

        {/* Help Section */}
        <div className="bg-blue-50 border border-blue-200 rounded-lg p-6">
          <div className="flex items-start">
            <svg className="h-6 w-6 text-blue-600 mr-3 flex-shrink-0" fill="none" viewBox="0 0 24 24" stroke="currentColor">
              <path
                strokeLinecap="round"
                strokeLinejoin="round"
                strokeWidth={2}
                d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z"
              />
            </svg>
            <div>
              <h3 className="text-lg font-semibold text-blue-900 mb-2">Need Help?</h3>
              <p className="text-blue-800 text-sm mb-3">
                Have questions about your subscription or billing? Our support team is here to help.
              </p>
              {usage.chat_support_enabled ? (
                <a
                  href="/support"
                  className="inline-block bg-blue-600 text-white px-4 py-2 rounded-lg text-sm font-medium hover:bg-blue-700 transition-colors"
                >
                  Open Live Chat
                </a>
              ) : (
                <a
                  href="mailto:support@wesplatform.com"
                  className="inline-block bg-blue-600 text-white px-4 py-2 rounded-lg text-sm font-medium hover:bg-blue-700 transition-colors"
                >
                  Email Support
                </a>
              )}
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
