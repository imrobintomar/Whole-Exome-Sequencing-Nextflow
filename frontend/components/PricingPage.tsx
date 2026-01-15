'use client';

import { useState, useEffect } from 'react';
import { billingApi, SubscriptionPlan } from '../lib/api';
import { Badge } from '@/components/ui/badge';

interface PricingPageProps {
  onSignIn?: () => void;
  onSelectPlan?: (planId: number) => void;
  isAuthenticated?: boolean;
}

export default function PricingPage({ onSignIn, onSelectPlan, isAuthenticated }: PricingPageProps) {
  const [plans, setPlans] = useState<SubscriptionPlan[]>([]);
  const [loading, setLoading] = useState(true);
  const [selectedPlanId, setSelectedPlanId] = useState<number | null>(null);
  const [billingPeriod, setBillingPeriod] = useState<'monthly' | 'annual'>('monthly');

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
        <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-purple-primary"></div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-gradient-to-br from-slate-50 via-white to-purple-50/30 py-16 px-4">
      <div className="max-w-7xl mx-auto">
        {/* Header */}
        <div className="text-center mb-12">
          <Badge variant="cyan" className="mb-6">
            Flexible Pricing
          </Badge>
          <h1 className="text-4xl md:text-5xl font-bold text-slate-900 mb-4">
            Choose Your Plan
          </h1>
          <p className="text-xl text-slate-600 max-w-2xl mx-auto">
            Unlock the full power of whole exome sequencing analysis with our flexible pricing plans
          </p>

          {/* Billing Period Toggle */}
          <div className="flex items-center justify-center gap-4 mt-8">
            <span className={`text-sm font-medium ${billingPeriod === 'monthly' ? 'text-slate-900' : 'text-slate-500'}`}>
              Monthly
            </span>
            <button
              onClick={() => setBillingPeriod(billingPeriod === 'monthly' ? 'annual' : 'monthly')}
              className="relative inline-flex h-8 w-14 items-center rounded-full bg-slate-200 transition-colors focus:outline-none focus:ring-2 focus:ring-cyan focus:ring-offset-2"
            >
              <span
                className={`inline-block h-6 w-6 transform rounded-full bg-purple-primary transition-transform ${
                  billingPeriod === 'annual' ? 'translate-x-7' : 'translate-x-1'
                }`}
              />
            </button>
            <span className={`text-sm font-medium ${billingPeriod === 'annual' ? 'text-slate-900' : 'text-slate-500'}`}>
              Annual
            </span>
            <Badge variant="success" className="ml-2">Save 20%</Badge>
          </div>
        </div>

        {/* Pricing Cards */}
        <div className="grid md:grid-cols-3 gap-8 max-w-6xl mx-auto mb-20">
          {plans.map((plan) => (
            <PricingCard
              key={plan.id}
              plan={plan}
              billingPeriod={billingPeriod}
              isPopular={plan.name === 'Basic'}
              onSelect={() => handleSelectPlan(plan.id)}
              loading={selectedPlanId === plan.id}
              isAuthenticated={isAuthenticated}
            />
          ))}
        </div>

        {/* Enterprise Section */}
        <div className="max-w-6xl mx-auto mb-20">
          <div className="bg-gradient-to-br from-purple-primary via-purple-light to-purple-dark rounded-2xl p-12 text-center text-white relative overflow-hidden">
            <div className="absolute inset-0 opacity-10">
              <div className="absolute top-10 right-20 w-96 h-96 bg-cyan rounded-full blur-3xl"></div>
            </div>
            <div className="relative">
              <Badge variant="cyan" className="mb-4">
                Enterprise
              </Badge>
              <h2 className="text-3xl font-bold mb-4">Need More? Contact Sales</h2>
              <p className="text-blue-100 text-lg mb-8 max-w-2xl mx-auto">
                Custom plans for institutions with high-volume sequencing needs
              </p>
              <div className="grid md:grid-cols-4 gap-6 mb-8">
                {[
                  'Unlimited jobs',
                  'Dedicated support',
                  'Custom SLAs',
                  'API access'
                ].map((feature, i) => (
                  <div key={i} className="flex items-center gap-2 justify-center">
                    <div className="w-2 h-2 bg-cyan rounded-full"></div>
                    <span className="text-sm">{feature}</span>
                  </div>
                ))}
              </div>
              <button
                onClick={onSignIn}
                className="bg-white text-purple-primary px-8 py-4 rounded-lg font-semibold text-lg hover:bg-blue-50 transition-colors"
              >
                Contact Sales
              </button>
            </div>
          </div>
        </div>

        {/* Features Comparison */}
        <div className="mt-20">
          <h2 className="text-3xl font-bold text-center text-slate-900 mb-12">Feature Comparison</h2>
          <div className="bg-white rounded-2xl shadow-lg overflow-hidden border-2 border-slate-100">
            <table className="w-full">
              <thead className="bg-gradient-to-r from-slate-50 to-purple-50/30">
                <tr>
                  <th className="px-6 py-4 text-left text-sm font-semibold text-slate-900">Feature</th>
                  {plans.map((plan) => (
                    <th key={plan.id} className="px-6 py-4 text-center text-sm font-semibold text-slate-900">
                      {plan.name}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody className="divide-y divide-slate-200">
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
          <h2 className="text-3xl font-bold text-center text-slate-900 mb-12">Frequently Asked Questions</h2>
          <div className="space-y-4">
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
            <FAQ
              question="How does annual billing work?"
              answer="With annual billing, you pay for 12 months upfront and save 20%. You'll have access to your plan benefits for the entire year."
            />
          </div>
        </div>

        {/* CTA Section */}
        {!isAuthenticated && (
          <div className="mt-20 text-center bg-gradient-to-br from-purple-primary via-purple-light to-purple-dark rounded-2xl p-12 relative overflow-hidden">
            <div className="absolute inset-0 opacity-10">
              <div className="absolute bottom-10 left-20 w-96 h-96 bg-cyan rounded-full blur-3xl"></div>
            </div>
            <div className="relative">
              <h2 className="text-3xl font-bold text-white mb-4">Ready to get started?</h2>
              <p className="text-blue-100 text-lg mb-8 max-w-2xl mx-auto">
                Create your free account today and start analyzing your WES data with our powerful platform
              </p>
              <button
                onClick={onSignIn}
                className="bg-white text-purple-primary px-8 py-4 rounded-lg font-semibold text-lg hover:bg-blue-50 transition-colors"
              >
                Sign Up Free
              </button>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}

// Helper Components
function PricingCard({ plan, billingPeriod, isPopular, onSelect, loading, isAuthenticated }: {
  plan: SubscriptionPlan;
  billingPeriod: 'monthly' | 'annual';
  isPopular: boolean;
  onSelect: () => void;
  loading: boolean;
  isAuthenticated?: boolean;
}) {
  const isFree = plan.price_usd === 0;
  const displayPrice = billingPeriod === 'annual' ? Math.round(plan.price_usd * 0.8) : plan.price_usd;

  return (
    <div
      className={`relative bg-white rounded-2xl shadow-lg overflow-hidden transition-all duration-300 ${
        isPopular ? 'ring-2 ring-cyan scale-105 shadow-cyan/20' : 'border-2 border-slate-100'
      }`}
    >
      {isPopular && (
        <div className="absolute top-0 right-0 bg-cyan text-white px-4 py-1 text-sm font-semibold rounded-bl-lg">
          Most Popular
        </div>
      )}

      <div className="p-8">
        <h3 className="text-2xl font-bold text-slate-900 mb-2">{plan.name}</h3>
        <div className="mb-6">
          <span className="text-5xl font-bold text-slate-900">${displayPrice}</span>
          <span className="text-slate-600">/month</span>
          {billingPeriod === 'annual' && !isFree && (
            <div className="text-sm text-teal mt-1">
              Save ${plan.price_usd * 12 - displayPrice * 12}/year
            </div>
          )}
        </div>

        <div className="mb-8">
          <div className="text-slate-600 mb-4">
            <span className="text-2xl font-bold text-slate-900">{plan.monthly_jobs_limit}</span> WES jobs per month
          </div>

          <ul className="space-y-3">
            {plan.features.map((feature, index) => (
              <li key={index} className="flex items-start">
                <svg className="h-5 w-5 text-cyan mr-2 flex-shrink-0 mt-0.5" fill="currentColor" viewBox="0 0 20 20">
                  <path
                    fillRule="evenodd"
                    d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z"
                    clipRule="evenodd"
                  />
                </svg>
                <span className="text-slate-600 text-sm">{feature}</span>
              </li>
            ))}
          </ul>
        </div>

        <button
          onClick={onSelect}
          disabled={loading}
          className={`w-full py-3 px-6 rounded-lg font-semibold transition-colors ${
            isFree
              ? 'bg-slate-100 text-slate-900 hover:bg-slate-200'
              : isPopular
              ? 'bg-cyan text-white hover:bg-cyan/90'
              : 'bg-purple-primary text-white hover:bg-purple-light'
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
    <tr className="hover:bg-slate-50/50 transition-colors">
      <td className="px-6 py-4 text-sm text-slate-900 font-medium">{feature}</td>
      {values.map((value, index) => (
        <td key={index} className="px-6 py-4 text-sm text-center">
          {value === '✓' ? (
            <span className="text-cyan font-bold text-lg">✓</span>
          ) : value === '-' ? (
            <span className="text-slate-300">-</span>
          ) : (
            <span className="text-slate-900">{value}</span>
          )}
        </td>
      ))}
    </tr>
  );
}

function FAQ({ question, answer }: { question: string; answer: string }) {
  const [isOpen, setIsOpen] = useState(false);

  return (
    <div className="bg-white rounded-xl shadow-md overflow-hidden border border-slate-200">
      <button
        onClick={() => setIsOpen(!isOpen)}
        className="w-full px-6 py-4 text-left flex items-center justify-between hover:bg-slate-50 transition-colors"
      >
        <span className="font-semibold text-slate-900">{question}</span>
        <svg
          className={`w-5 h-5 text-cyan transform transition-transform ${isOpen ? 'rotate-180' : ''}`}
          fill="none"
          viewBox="0 0 24 24"
          stroke="currentColor"
        >
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
        </svg>
      </button>
      {isOpen && (
        <div className="px-6 pb-4 text-slate-600 border-t border-slate-100 pt-4">
          {answer}
        </div>
      )}
    </div>
  );
}
