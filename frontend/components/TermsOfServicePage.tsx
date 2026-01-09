'use client';

import { ArrowLeft } from 'lucide-react';
import { Button } from '@/components/ui/button';

interface TermsOfServicePageProps {
  onBack: () => void;
}

export default function TermsOfServicePage({ onBack }: TermsOfServicePageProps) {
  return (
    <div className="min-h-screen bg-gradient-to-br from-slate-50 via-white to-slate-50">
      <div className="max-w-4xl mx-auto px-4 py-12 sm:px-6 lg:px-8">
        <Button
          variant="ghost"
          onClick={onBack}
          className="mb-8 flex items-center gap-2 text-slate-600 hover:text-slate-900"
        >
          <ArrowLeft className="h-4 w-4" />
          Back to Home
        </Button>

        <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-8 md:p-12">
          <h1 className="text-4xl font-bold text-slate-900 mb-4">Terms of Service</h1>
          <p className="text-sm text-slate-600 mb-8">Last Updated: January 9, 2026</p>

          <div className="prose prose-slate max-w-none">
            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">1. Acceptance of Terms</h2>
              <p className="text-slate-700 mb-4">
                Welcome to ATGCFLOW. By accessing or using our whole exome sequencing analysis platform, you agree to be bound by these Terms of Service ("Terms"). If you do not agree to these Terms, please do not use our platform.
              </p>
              <p className="text-slate-700">
                These Terms constitute a legally binding agreement between you and ATGCFLOW. We reserve the right to modify these Terms at any time, and your continued use of the platform constitutes acceptance of any changes.
              </p>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">2. Research Use Only</h2>
              <div className="bg-red-50 border-2 border-red-300 rounded-lg p-6 mb-4">
                <p className="text-red-900 font-bold text-lg mb-3">CRITICAL NOTICE - NOT FOR CLINICAL USE</p>
                <p className="text-red-800 mb-2">
                  ATGCFLOW IS A RESEARCH PLATFORM ONLY. This platform is:
                </p>
                <ul className="list-disc pl-6 text-red-800 space-y-1">
                  <li>NOT approved for clinical diagnosis or treatment decisions</li>
                  <li>NOT a substitute for professional medical advice</li>
                  <li>NOT validated for clinical laboratory use</li>
                  <li>NOT intended to diagnose, treat, cure, or prevent any disease</li>
                </ul>
                <p className="text-red-800 mt-3">
                  All results must be validated through certified clinical laboratories before any clinical application.
                </p>
              </div>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">3. User Accounts and Registration</h2>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">3.1 Account Creation</h3>
              <p className="text-slate-700 mb-4">
                To use ATGCFLOW, you must create an account by providing accurate and complete information. You agree to:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>Provide accurate, current, and complete information during registration</li>
                <li>Maintain and promptly update your account information</li>
                <li>Maintain the security of your password and account credentials</li>
                <li>Accept responsibility for all activities under your account</li>
                <li>Notify us immediately of any unauthorized use of your account</li>
              </ul>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">3.2 User Eligibility</h3>
              <p className="text-slate-700 mb-4">
                You must be at least 18 years old to use ATGCFLOW. By using the platform, you represent and warrant that you meet this age requirement and have the legal capacity to enter into these Terms.
              </p>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">4. Acceptable Use</h2>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">4.1 Permitted Uses</h3>
              <p className="text-slate-700 mb-4">
                You may use ATGCFLOW for:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>Research and educational purposes</li>
                <li>Whole exome sequencing data analysis</li>
                <li>Variant annotation and classification</li>
                <li>Gene panel analysis for research studies</li>
              </ul>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">4.2 Prohibited Activities</h3>
              <p className="text-slate-700 mb-4">
                You agree NOT to:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>Use the platform for clinical diagnosis or patient care decisions</li>
                <li>Upload data without proper authorization and consent</li>
                <li>Attempt to reverse engineer, decompile, or disassemble the platform</li>
                <li>Use the platform to transmit malicious code, viruses, or harmful data</li>
                <li>Violate any applicable laws, regulations, or third-party rights</li>
                <li>Share your account credentials with unauthorized users</li>
                <li>Attempt to gain unauthorized access to other users' data or accounts</li>
                <li>Use automated tools to scrape or download data without permission</li>
                <li>Interfere with or disrupt the platform's operation or servers</li>
              </ul>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">5. Data Ownership and Usage Rights</h2>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">5.1 Your Data</h3>
              <p className="text-slate-700 mb-4">
                You retain all ownership rights to the genomic data you upload to ATGCFLOW. By uploading data, you grant us a limited license to:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>Process and analyze your data to provide our services</li>
                <li>Store your data securely on our servers</li>
                <li>Generate analysis results, annotations, and classifications</li>
              </ul>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">5.2 Data Responsibility</h3>
              <p className="text-slate-700 mb-4">
                You represent and warrant that:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>You have all necessary rights, permissions, and consents to upload the data</li>
                <li>The data complies with all applicable laws and regulations (e.g., HIPAA, GDPR)</li>
                <li>You have obtained appropriate informed consent from data subjects</li>
                <li>The data does not violate any third-party intellectual property rights</li>
              </ul>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">6. Service Availability and Modifications</h2>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">6.1 Service Availability</h3>
              <p className="text-slate-700 mb-4">
                We strive to maintain high availability of ATGCFLOW but do not guarantee uninterrupted access. The platform may be unavailable due to:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>Scheduled maintenance and updates</li>
                <li>Technical issues or system failures</li>
                <li>Third-party service disruptions</li>
                <li>Force majeure events</li>
              </ul>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">6.2 Service Modifications</h3>
              <p className="text-slate-700 mb-4">
                We reserve the right to modify, suspend, or discontinue any aspect of ATGCFLOW at any time, with or without notice. We are not liable for any modifications, suspensions, or discontinuations of the service.
              </p>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">7. Intellectual Property Rights</h2>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">7.1 Platform Ownership</h3>
              <p className="text-slate-700 mb-4">
                ATGCFLOW, including all software, code, algorithms, designs, text, graphics, and other materials, is owned by us and protected by copyright, trademark, and other intellectual property laws.
              </p>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">7.2 License Grant</h3>
              <p className="text-slate-700 mb-4">
                We grant you a limited, non-exclusive, non-transferable, revocable license to access and use ATGCFLOW for research purposes only, subject to these Terms.
              </p>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">7.3 Third-Party Tools</h3>
              <p className="text-slate-700 mb-4">
                ATGCFLOW integrates various third-party bioinformatics tools (GATK, BWA, ANNOVAR, etc.). Your use of these tools is subject to their respective licenses and terms.
              </p>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">8. Disclaimers and Limitations of Liability</h2>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">8.1 "AS IS" Service</h3>
              <div className="bg-amber-50 border border-amber-200 rounded-lg p-4 mb-4">
                <p className="text-amber-900 mb-2">
                  ATGCFLOW IS PROVIDED "AS IS" AND "AS AVAILABLE" WITHOUT WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO:
                </p>
                <ul className="list-disc pl-6 text-amber-800 space-y-1">
                  <li>Accuracy, completeness, or reliability of analysis results</li>
                  <li>Fitness for any particular purpose (especially clinical use)</li>
                  <li>Non-infringement of third-party rights</li>
                  <li>Uninterrupted or error-free operation</li>
                </ul>
              </div>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">8.2 Limitation of Liability</h3>
              <p className="text-slate-700 mb-4">
                TO THE MAXIMUM EXTENT PERMITTED BY LAW, WE SHALL NOT BE LIABLE FOR ANY:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>Indirect, incidental, special, consequential, or punitive damages</li>
                <li>Loss of profits, data, use, goodwill, or other intangible losses</li>
                <li>Damages resulting from reliance on analysis results</li>
                <li>Unauthorized access to or alteration of your data</li>
                <li>Errors, mistakes, or inaccuracies in analysis results</li>
              </ul>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">9. Indemnification</h2>
              <p className="text-slate-700 mb-4">
                You agree to indemnify, defend, and hold harmless ATGCFLOW and its affiliates, officers, directors, employees, and agents from any claims, liabilities, damages, losses, costs, or expenses arising from:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>Your use or misuse of the platform</li>
                <li>Your violation of these Terms</li>
                <li>Your violation of any rights of third parties</li>
                <li>Your data uploads or lack of proper consent</li>
              </ul>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">10. Termination</h2>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">10.1 Termination by You</h3>
              <p className="text-slate-700 mb-4">
                You may terminate your account at any time by contacting us or using the account deletion feature. Upon termination, you may request deletion of your data.
              </p>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">10.2 Termination by Us</h3>
              <p className="text-slate-700 mb-4">
                We reserve the right to suspend or terminate your account immediately, without notice, if you:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>Violate these Terms</li>
                <li>Use the platform for clinical purposes</li>
                <li>Engage in prohibited activities</li>
                <li>Pose a security or legal risk</li>
              </ul>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">11. Governing Law and Dispute Resolution</h2>
              <p className="text-slate-700 mb-4">
                These Terms shall be governed by and construed in accordance with the laws of the applicable jurisdiction, without regard to conflict of law principles. Any disputes arising from these Terms or your use of ATGCFLOW shall be resolved through binding arbitration or in the courts of the applicable jurisdiction.
              </p>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">12. Miscellaneous</h2>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">12.1 Entire Agreement</h3>
              <p className="text-slate-700 mb-4">
                These Terms, together with our Privacy Policy, constitute the entire agreement between you and ATGCFLOW.
              </p>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">12.2 Severability</h3>
              <p className="text-slate-700 mb-4">
                If any provision of these Terms is found to be unenforceable, the remaining provisions will remain in full force and effect.
              </p>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">12.3 No Waiver</h3>
              <p className="text-slate-700 mb-4">
                Our failure to enforce any right or provision of these Terms shall not be deemed a waiver of such right or provision.
              </p>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">13. Contact Information</h2>
              <p className="text-slate-700 mb-4">
                If you have questions about these Terms of Service, please contact us:
              </p>
              <div className="bg-slate-50 border border-slate-200 rounded-lg p-4">
                <p className="text-slate-700">
                  <strong>Email:</strong> legal@atgcflow.com<br />
                  <strong>Support:</strong> support@atgcflow.com<br />
                  <strong>Platform:</strong> ATGCFLOW<br />
                  <strong>Type:</strong> Research Platform - Not for Clinical Use
                </p>
              </div>
            </section>

            <div className="bg-blue-50 border border-blue-200 rounded-lg p-6 mt-8">
              <p className="text-blue-900 font-semibold mb-2">
                By using ATGCFLOW, you acknowledge that you have read, understood, and agree to be bound by these Terms of Service.
              </p>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
