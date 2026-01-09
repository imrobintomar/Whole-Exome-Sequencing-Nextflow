'use client';

import { ArrowLeft } from 'lucide-react';
import { Button } from '@/components/ui/button';

interface PrivacyPolicyPageProps {
  onBack: () => void;
}

export default function PrivacyPolicyPage({ onBack }: PrivacyPolicyPageProps) {
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
          <h1 className="text-4xl font-bold text-slate-900 mb-4">Privacy Policy</h1>
          <p className="text-sm text-slate-600 mb-8">Last Updated: January 9, 2026</p>

          <div className="prose prose-slate max-w-none">
            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">1. Introduction</h2>
              <p className="text-slate-700 mb-4">
                Welcome to ATGCFLOW ("we," "our," or "us"). This Privacy Policy explains how we collect, use, disclose, and safeguard your information when you use our whole exome sequencing analysis platform. Please read this privacy policy carefully.
              </p>
              <p className="text-slate-700">
                By using ATGCFLOW, you agree to the collection and use of information in accordance with this policy. If you do not agree with our policies and practices, please do not use our platform.
              </p>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">2. Information We Collect</h2>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">2.1 Personal Information</h3>
              <p className="text-slate-700 mb-4">
                We collect information that you provide directly to us, including:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>Account information (name, email address, password)</li>
                <li>Profile information</li>
                <li>Contact information</li>
                <li>Communication preferences</li>
              </ul>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">2.2 Genomic Data</h3>
              <p className="text-slate-700 mb-4">
                When you upload sequencing files for analysis, we process:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>FASTQ files (raw sequencing data)</li>
                <li>VCF files (variant call data)</li>
                <li>BAM files (aligned sequencing data)</li>
                <li>Analysis results and annotations</li>
              </ul>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">2.3 Usage Information</h3>
              <p className="text-slate-700 mb-4">
                We automatically collect certain information about your device and usage patterns:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>Log data (IP address, browser type, operating system)</li>
                <li>Usage statistics (features accessed, time spent on platform)</li>
                <li>Performance metrics (pipeline execution times, error rates)</li>
              </ul>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">3. How We Use Your Information</h2>
              <p className="text-slate-700 mb-4">
                We use the information we collect for the following purposes:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>To provide, maintain, and improve our sequencing analysis services</li>
                <li>To process your genomic data and generate analysis results</li>
                <li>To communicate with you about your account and analysis results</li>
                <li>To monitor and analyze usage patterns and platform performance</li>
                <li>To detect, prevent, and address technical issues and security threats</li>
                <li>To comply with legal obligations and protect our rights</li>
              </ul>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">4. Data Security</h2>
              <p className="text-slate-700 mb-4">
                We implement appropriate technical and organizational measures to protect your data:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>Encryption of data in transit and at rest</li>
                <li>Secure authentication using Firebase Authentication</li>
                <li>Regular security audits and vulnerability assessments</li>
                <li>Access controls and role-based permissions</li>
                <li>Secure data storage and backup procedures</li>
              </ul>
              <p className="text-slate-700 mb-4">
                However, no method of transmission over the Internet or electronic storage is 100% secure. While we strive to protect your data, we cannot guarantee absolute security.
              </p>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">5. Data Retention</h2>
              <p className="text-slate-700 mb-4">
                We retain your information for as long as necessary to fulfill the purposes outlined in this Privacy Policy:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>Account information: Until you request account deletion</li>
                <li>Genomic data: Until you delete the analysis or request data removal</li>
                <li>Usage logs: Typically retained for 90 days for troubleshooting purposes</li>
              </ul>
              <p className="text-slate-700">
                You may request deletion of your data at any time by contacting us.
              </p>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">6. Data Sharing and Disclosure</h2>
              <p className="text-slate-700 mb-4">
                We do not sell, trade, or rent your personal or genomic data to third parties. We may share your information only in the following circumstances:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li><strong>With Your Consent:</strong> When you explicitly authorize us to share your data</li>
                <li><strong>Service Providers:</strong> With third-party vendors who assist in operating our platform (e.g., cloud hosting)</li>
                <li><strong>Legal Requirements:</strong> When required by law, regulation, or legal process</li>
                <li><strong>Protection of Rights:</strong> To protect our rights, property, or safety, or that of our users</li>
              </ul>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">7. Your Rights</h2>
              <p className="text-slate-700 mb-4">
                You have the following rights regarding your personal information:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li><strong>Access:</strong> Request access to your personal data</li>
                <li><strong>Correction:</strong> Request correction of inaccurate data</li>
                <li><strong>Deletion:</strong> Request deletion of your data</li>
                <li><strong>Data Portability:</strong> Request a copy of your data in a structured format</li>
                <li><strong>Withdraw Consent:</strong> Withdraw consent for data processing at any time</li>
              </ul>
              <p className="text-slate-700">
                To exercise these rights, please contact us using the information provided below.
              </p>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">8. Research Use Disclaimer</h2>
              <div className="bg-amber-50 border border-amber-200 rounded-lg p-4 mb-4">
                <p className="text-amber-900 font-semibold mb-2">IMPORTANT NOTICE:</p>
                <p className="text-amber-800">
                  ATGCFLOW is a research platform and is NOT intended for clinical use, diagnosis, or treatment decisions. All analysis results should be validated through certified clinical laboratories before any clinical application.
                </p>
              </div>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">9. Children's Privacy</h2>
              <p className="text-slate-700 mb-4">
                ATGCFLOW is not intended for use by individuals under the age of 18. We do not knowingly collect personal information from children. If you become aware that a child has provided us with personal information, please contact us.
              </p>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">10. International Data Transfers</h2>
              <p className="text-slate-700 mb-4">
                Your information may be transferred to and processed in countries other than your country of residence. These countries may have different data protection laws. By using ATGCFLOW, you consent to such transfers.
              </p>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">11. Changes to This Privacy Policy</h2>
              <p className="text-slate-700 mb-4">
                We may update this Privacy Policy from time to time. We will notify you of any changes by posting the new Privacy Policy on this page and updating the "Last Updated" date. You are advised to review this Privacy Policy periodically for any changes.
              </p>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">12. Contact Us</h2>
              <p className="text-slate-700 mb-4">
                If you have any questions about this Privacy Policy or our data practices, please contact us:
              </p>
              <div className="bg-slate-50 border border-slate-200 rounded-lg p-4">
                <p className="text-slate-700">
                  <strong>Email:</strong> privacy@atgcflow.com<br />
                  <strong>Platform:</strong> ATGCFLOW<br />
                  <strong>Type:</strong> Research Platform - Not for Clinical Use
                </p>
              </div>
            </section>
          </div>
        </div>
      </div>
    </div>
  );
}
