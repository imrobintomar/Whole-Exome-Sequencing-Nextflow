'use client';

import { useState, useEffect } from 'react';
import { billingApi, SubscriptionPlan } from '../lib/api';

interface PricingPageProps {
  onSignIn?: () => void;
  onSelectPlan?: (planId: number) => void;
  isAuthenticated?: boolean;
}

export default function PricingPage({ onSignIn, onSelectPlan, isAuthenticated }: PricingPageProps) {
  const [plans, setPlans] = useState<SubscriptionPlan[]>([]);
  const [loading, setLoading] = useState(true);
  const [selectedPlanId, setSelectedPlanId] = useState<number | null>(null);

  useEffect(() => {
    loadPlans();
  }, []);

  const loadPlans = async () => {
    try {
      const plansData = await billingApi.getPlans();
      setPlans(plansData);
    } catch (error) {
      console.error('Failed to load plans:', error);
    } finally {
      setLoading(false);
    }
  };

  const handleSelectPlan = async (planId: number) => {
    if (!isAuthenticated) {
      onSignIn?.();
      return;
    }

    setSelectedPlanId(planId);
    if (onSelectPlan) {
      onSelectPlan(planId);
    } else {
      // Default behavior: create checkout session
      try {
        const { checkout_url } = await billingApi.createCheckout(planId);
        window.location.href = checkout_url;
      } catch (error) {
        console.error('Failed to create checkout:', error);
        setSelectedPlanId(null);
        alert('Failed to start checkout. Please try again.');
      }
    }
  };

  if (loading) {
    return (
      <div className="flex items-center justify-center min-h-screen">
        <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-600"></div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 via-white to-purple-50 py-12 px-4">
      <div className="max-w-7xl mx-auto">
        {/* Header */}
        <div className="text-center mb-16">
          <h1 className="text-4xl md:text-5xl font-bold text-gray-900 mb-4">
            Choose Your Plan
          </h1>
          <p className="text-xl text-gray-600 max-w-2xl mx-auto">
            Unlock the full power of whole exome sequencing analysis with our flexible pricing plans
          </p>
        </div>

        {/* Pricing Cards */}
        <div className="grid md:grid-cols-3 gap-8 max-w-6xl mx-auto">
          {plans.map((plan) => (
            <PricingCard
              key={plan.id}
              plan={plan}
              isPopular={plan.name === 'Basic'}
              onSelect={() => handleSelectPlan(plan.id)}
              loading={selectedPlanId === plan.id}
              isAuthenticated={isAuthenticated}
            />
          ))}
        </div>

        {/* Features Comparison */}
        <div className="mt-20">
          <h2 className="text-3xl font-bold text-center text-gray-900 mb-12">Feature Comparison</h2>
          <div className="bg-white rounded-lg shadow-lg overflow-hidden">
            <table className="w-full">
              <thead className="bg-gray-50">
                <tr>
                  <th className="px-6 py-4 text-left text-sm font-semibold text-gray-900">Feature</th>
                  {plans.map((plan) => (
                    <th key={plan.id} className="px-6 py-4 text-center text-sm font-semibold text-gray-900">
                      {plan.name}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody className="divide-y divide-gray-200">
                <FeatureRow
                  feature="Monthly WES Jobs"
                  values={plans.map((p) => p.monthly_jobs_limit.toString())}
                />
                <FeatureRow
                  feature="Analysis Pipeline"
                  values={plans.map(() => '✓')}
                />
                <FeatureRow
                  feature="ACMG Classification"
                  values={plans.map(() => '✓')}
                />
                <FeatureRow
                  feature="Gene Panel Filtering"
                  values={plans.map(() => '✓')}
                />
                <FeatureRow
                  feature="IGV Browser"
                  values={plans.map(() => '✓')}
                />
                <FeatureRow
                  feature="Result Downloads"
                  values={plans.map(() => '✓')}
                />
                <FeatureRow
                  feature="Live Chat Support"
                  values={plans.map((p) => (p.chat_support ? '✓' : '-'))}
                />
                <FeatureRow
                  feature="Priority Processing"
                  values={['-', '✓', '✓']}
                />
                <FeatureRow
                  feature="Turnaround Time"
                  values={['48 hours', '24 hours', '12 hours']}
                />
                <FeatureRow
                  feature="API Access"
                  values={['-', '-', '✓']}
                />
                <FeatureRow
                  feature="Dedicated Support"
                  values={['-', '-', '✓']}
                />
              </tbody>
            </table>
          </div>
        </div>

        {/* FAQ Section */}
        <div className="mt-20 max-w-3xl mx-auto">
          <h2 className="text-3xl font-bold text-center text-gray-900 mb-12">Frequently Asked Questions</h2>
          <div className="space-y-6">
            <FAQ
              question="Can I change my plan later?"
              answer="Yes! You can upgrade or downgrade your plan at any time from your account settings. Changes take effect immediately, and we'll prorate any charges."
            />
            <FAQ
              question="What happens if I exceed my monthly job limit?"
              answer="Once you reach your monthly limit, you'll need to wait until the next billing cycle or upgrade to a higher plan to submit more jobs."
            />
            <FAQ
              question="Do unused jobs roll over to the next month?"
              answer="No, job credits reset at the beginning of each monthly billing cycle."
            />
            <FAQ
              question="Can I cancel my subscription?"
              answer="Yes, you can cancel anytime from your account settings. Your subscription will remain active until the end of your current billing period."
            />
            <FAQ
              question="What file formats do you support?"
              answer="We support standard FASTQ files (.fastq.gz or .fq.gz) for paired-end sequencing data."
            />
          </div>
        </div>

        {/* CTA Section */}
        {!isAuthenticated && (
          <div className="mt-20 text-center bg-blue-600 rounded-2xl p-12">
            <h2 className="text-3xl font-bold text-white mb-4">Ready to get started?</h2>
            <p className="text-blue-100 text-lg mb-8 max-w-2xl mx-auto">
              Create your free account today and start analyzing your WES data with our powerful platform
            </p>
            <button
              onClick={onSignIn}
              className="bg-white text-blue-600 px-8 py-4 rounded-lg font-semibold text-lg hover:bg-blue-50 transition-colors"
            >
              Sign Up Free
            </button>
          </div>
        )}
      </div>
    </div>
  );
}

// Helper Components
function PricingCard({ plan, isPopular, onSelect, loading, isAuthenticated }: {
  plan: SubscriptionPlan;
  isPopular: boolean;
  onSelect: () => void;
  loading: boolean;
  isAuthenticated?: boolean;
}) {
  const isFree = plan.price_usd === 0;

  return (
    <div
      className={`relative bg-white rounded-2xl shadow-lg overflow-hidden ${
        isPopular ? 'ring-2 ring-blue-600 scale-105' : ''
      }`}
    >
      {isPopular && (
        <div className="absolute top-0 right-0 bg-blue-600 text-white px-4 py-1 text-sm font-semibold rounded-bl-lg">
          Most Popular
        </div>
      )}

      <div className="p-8">
        <h3 className="text-2xl font-bold text-gray-900 mb-2">{plan.name}</h3>
        <div className="mb-6">
          <span className="text-5xl font-bold text-gray-900">${plan.price_usd}</span>
          <span className="text-gray-600">/month</span>
        </div>

        <div className="mb-8">
          <div className="text-gray-600 mb-4">
            <span className="text-2xl font-bold text-gray-900">{plan.monthly_jobs_limit}</span> WES jobs per month
          </div>

          <ul className="space-y-3">
            {plan.features.map((feature, index) => (
              <li key={index} className="flex items-start">
                <svg className="h-5 w-5 text-green-500 mr-2 flex-shrink-0 mt-0.5" fill="currentColor" viewBox="0 0 20 20">
                  <path
                    fillRule="evenodd"
                    d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z"
                    clipRule="evenodd"
                  />
                </svg>
                <span className="text-gray-600 text-sm">{feature}</span>
              </li>
            ))}
          </ul>
        </div>

        <button
          onClick={onSelect}
          disabled={loading}
          className={`w-full py-3 px-6 rounded-lg font-semibold transition-colors ${
            isFree
              ? 'bg-gray-100 text-gray-900 hover:bg-gray-200'
              : isPopular
              ? 'bg-blue-600 text-white hover:bg-blue-700'
              : 'bg-purple-600 text-white hover:bg-purple-700'
          } ${loading ? 'opacity-50 cursor-not-allowed' : ''}`}
        >
          {loading ? 'Processing...' : isFree ? 'Get Started Free' : isAuthenticated ? 'Subscribe Now' : 'Sign Up'}
        </button>
      </div>
    </div>
  );
}

function FeatureRow({ feature, values }: { feature: string; values: string[] }) {
  return (
    <tr>
      <td className="px-6 py-4 text-sm text-gray-900 font-medium">{feature}</td>
      {values.map((value, index) => (
        <td key={index} className="px-6 py-4 text-sm text-center">
          {value === '✓' ? (
            <span className="text-green-600 font-bold">✓</span>
          ) : value === '-' ? (
            <span className="text-gray-300">-</span>
          ) : (
            <span className="text-gray-900">{value}</span>
          )}
        </td>
      ))}
    </tr>
  );
}

function FAQ({ question, answer }: { question: string; answer: string }) {
  const [isOpen, setIsOpen] = useState(false);

  return (
    <div className="bg-white rounded-lg shadow-md overflow-hidden">
      <button
        onClick={() => setIsOpen(!isOpen)}
        className="w-full px-6 py-4 text-left flex items-center justify-between hover:bg-gray-50 transition-colors"
      >
        <span className="font-semibold text-gray-900">{question}</span>
        <svg
          className={`w-5 h-5 text-gray-500 transform transition-transform ${isOpen ? 'rotate-180' : ''}`}
          fill="none"
          viewBox="0 0 24 24"
          stroke="currentColor"
        >
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </button>
      {isOpen && (
        <div className="px-6 pb-4 text-gray-600">
          {answer}
        </div>
      )}
    </div>
  );
}
