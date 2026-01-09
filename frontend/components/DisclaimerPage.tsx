'use client';

import { ArrowLeft, AlertTriangle, ShieldAlert } from 'lucide-react';
import { Button } from '@/components/ui/button';

interface DisclaimerPageProps {
  onBack: () => void;
}

export default function DisclaimerPage({ onBack }: DisclaimerPageProps) {
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
          <div className="flex items-center gap-3 mb-4">
            <ShieldAlert className="h-10 w-10 text-red-600" />
            <h1 className="text-4xl font-bold text-slate-900">Disclaimer</h1>
          </div>
          <p className="text-sm text-slate-600 mb-8">Last Updated: January 9, 2026</p>

          <div className="prose prose-slate max-w-none">
            {/* Critical Warning */}
            <div className="bg-red-50 border-2 border-red-400 rounded-xl p-6 mb-8">
              <div className="flex items-start gap-3">
                <AlertTriangle className="h-6 w-6 text-red-600 flex-shrink-0 mt-1" />
                <div>
                  <h2 className="text-2xl font-bold text-red-900 mb-3">CRITICAL WARNING - NOT FOR CLINICAL USE</h2>
                  <p className="text-red-800 font-semibold text-lg mb-3">
                    ATGCFLOW IS A RESEARCH PLATFORM ONLY AND IS NOT INTENDED FOR:
                  </p>
                  <ul className="list-disc pl-6 text-red-800 space-y-2 mb-3">
                    <li className="font-semibold">Clinical diagnosis or treatment decisions</li>
                    <li className="font-semibold">Patient care or medical advice</li>
                    <li className="font-semibold">Diagnostic testing in a clinical laboratory</li>
                    <li className="font-semibold">Replacing professional medical consultation</li>
                  </ul>
                  <p className="text-red-900 font-bold text-lg">
                    DO NOT make any medical decisions based on results from this platform. All findings must be validated through certified clinical laboratories before any clinical application.
                  </p>
                </div>
              </div>
            </div>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">1. General Disclaimer</h2>
              <p className="text-slate-700 mb-4">
                The information provided by ATGCFLOW is for research and educational purposes only. While we strive to provide accurate and reliable analysis results using industry-standard bioinformatics tools and best practices, we make no representations or warranties of any kind regarding:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>The accuracy, completeness, or reliability of analysis results</li>
                <li>The suitability of results for any particular purpose</li>
                <li>The interpretation or clinical significance of variants</li>
                <li>The quality or validity of uploaded sequencing data</li>
              </ul>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">2. Research Use Only</h2>
              <div className="bg-amber-50 border border-amber-300 rounded-lg p-4 mb-4">
                <p className="text-amber-900 font-semibold mb-2">Research Platform Status:</p>
                <p className="text-amber-800 mb-2">
                  ATGCFLOW is designed and maintained exclusively for research purposes. It has NOT been:
                </p>
                <ul className="list-disc pl-6 text-amber-800 space-y-1">
                  <li>Validated for clinical use or diagnostic purposes</li>
                  <li>Approved by any regulatory authority (FDA, CLIA, CAP, etc.)</li>
                  <li>Certified for use in clinical decision-making</li>
                  <li>Subjected to clinical validation studies</li>
                </ul>
              </div>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">3. No Medical Advice</h2>
              <p className="text-slate-700 mb-4">
                ATGCFLOW does NOT provide medical advice, diagnosis, or treatment recommendations. The platform:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>Does not replace consultation with qualified healthcare professionals</li>
                <li>Should not be used to diagnose or treat any medical condition</li>
                <li>Does not constitute a doctor-patient relationship</li>
                <li>Cannot be relied upon for clinical decision-making</li>
              </ul>
              <p className="text-slate-700 mb-4 font-semibold">
                Always seek the advice of qualified healthcare providers with any questions regarding medical conditions or genetic findings. Never disregard professional medical advice or delay seeking it because of information from ATGCFLOW.
              </p>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">4. Analysis Limitations</h2>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">4.1 Technical Limitations</h3>
              <p className="text-slate-700 mb-4">
                Whole exome sequencing and bioinformatics analysis have inherent limitations:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li><strong>Coverage Limitations:</strong> Not all genomic regions are equally covered</li>
                <li><strong>Variant Detection:</strong> Some variant types may not be detected (e.g., large structural variants, repeat expansions)</li>
                <li><strong>False Positives/Negatives:</strong> Analysis may produce false positive or false negative results</li>
                <li><strong>Annotation Accuracy:</strong> Variant annotations depend on database completeness and accuracy</li>
                <li><strong>Interpretation Complexity:</strong> Genetic variants have complex interpretations that evolve with scientific knowledge</li>
              </ul>

              <h3 className="text-xl font-semibold text-slate-800 mb-3">4.2 ACMG Classification Limitations</h3>
              <p className="text-slate-700 mb-4">
                The ACMG variant classification provided by ATGCFLOW:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>Is automated and may not capture all relevant clinical context</li>
                <li>Should be reviewed by certified genetic counselors or clinical geneticists</li>
                <li>May change as new scientific evidence emerges</li>
                <li>Does not replace manual curation by experts</li>
                <li>May not follow the latest ACMG/AMP guidelines updates</li>
              </ul>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">5. Third-Party Tools and Databases</h2>
              <p className="text-slate-700 mb-4">
                ATGCFLOW integrates various third-party bioinformatics tools and databases:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li><strong>GATK (Genome Analysis Toolkit):</strong> For variant calling and processing</li>
                <li><strong>BWA (Burrows-Wheeler Aligner):</strong> For read alignment</li>
                <li><strong>ANNOVAR:</strong> For variant annotation</li>
                <li><strong>gnomAD, ClinVar, dbSNP:</strong> For population frequencies and clinical annotations</li>
                <li><strong>Prediction Tools:</strong> CADD, REVEL, SIFT, PolyPhen</li>
              </ul>
              <p className="text-slate-700 mb-4">
                We are not responsible for the accuracy, reliability, or availability of these third-party tools and databases. Each tool has its own limitations and should be used according to its respective documentation and licenses.
              </p>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">6. Data Quality and User Responsibility</h2>
              <p className="text-slate-700 mb-4">
                Analysis results are highly dependent on input data quality. Users are responsible for:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>Ensuring sequencing data meets quality standards</li>
                <li>Verifying proper sample preparation and sequencing protocols</li>
                <li>Confirming appropriate reference genome usage (hg38)</li>
                <li>Validating results through orthogonal methods</li>
                <li>Obtaining all necessary consents and permissions for data use</li>
              </ul>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">7. No Warranty</h2>
              <div className="bg-slate-100 border border-slate-300 rounded-lg p-4 mb-4">
                <p className="text-slate-800 mb-2">
                  ATGCFLOW IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF:
                </p>
                <ul className="list-disc pl-6 text-slate-700 space-y-1">
                  <li>Merchantability</li>
                  <li>Fitness for a particular purpose</li>
                  <li>Non-infringement</li>
                  <li>Accuracy or completeness</li>
                  <li>Reliability or availability</li>
                </ul>
              </div>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">8. Limitation of Liability</h2>
              <p className="text-slate-700 mb-4">
                TO THE MAXIMUM EXTENT PERMITTED BY LAW, ATGCFLOW AND ITS OPERATORS SHALL NOT BE LIABLE FOR ANY DAMAGES ARISING FROM:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>Use or inability to use the platform</li>
                <li>Reliance on analysis results</li>
                <li>Errors, inaccuracies, or omissions in results</li>
                <li>Data loss or corruption</li>
                <li>Service interruptions or downtime</li>
                <li>Any medical decisions or actions taken based on platform results</li>
              </ul>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">9. Regulatory Compliance</h2>
              <p className="text-slate-700 mb-4">
                Users are responsible for ensuring their use of ATGCFLOW complies with all applicable regulations:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li><strong>HIPAA (USA):</strong> If handling protected health information</li>
                <li><strong>GDPR (EU):</strong> If processing EU citizens' genetic data</li>
                <li><strong>CLIA/CAP (USA):</strong> Clinical laboratories must use certified platforms</li>
                <li><strong>IRB Approval:</strong> Research studies may require institutional review board approval</li>
                <li><strong>Informed Consent:</strong> Proper consent must be obtained for genetic testing</li>
              </ul>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">10. Updates and Changes</h2>
              <p className="text-slate-700 mb-4">
                This disclaimer may be updated periodically to reflect:
              </p>
              <ul className="list-disc pl-6 mb-4 text-slate-700 space-y-2">
                <li>Changes in platform capabilities</li>
                <li>Updates to integrated tools and databases</li>
                <li>New regulatory requirements</li>
                <li>Evolving best practices in bioinformatics</li>
              </ul>
              <p className="text-slate-700">
                Continued use of ATGCFLOW after updates constitutes acceptance of the revised disclaimer.
              </p>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">11. Validation Requirements</h2>
              <div className="bg-blue-50 border border-blue-300 rounded-lg p-4 mb-4">
                <p className="text-blue-900 font-semibold mb-2">Before Clinical Application:</p>
                <p className="text-blue-800 mb-2">
                  Any finding from ATGCFLOW that may have clinical implications MUST be:
                </p>
                <ul className="list-disc pl-6 text-blue-800 space-y-1">
                  <li>Confirmed through a CLIA-certified clinical laboratory</li>
                  <li>Validated using orthogonal sequencing methods (e.g., Sanger sequencing)</li>
                  <li>Interpreted by board-certified genetic counselors or clinical geneticists</li>
                  <li>Considered in the context of patient phenotype and family history</li>
                  <li>Discussed with the patient's healthcare provider</li>
                </ul>
              </div>
            </section>

            <section className="mb-8">
              <h2 className="text-2xl font-semibold text-slate-900 mb-4">12. Contact Information</h2>
              <p className="text-slate-700 mb-4">
                If you have questions about this disclaimer, please contact us:
              </p>
              <div className="bg-slate-50 border border-slate-200 rounded-lg p-4">
                <p className="text-slate-700">
                  <strong>Email:</strong> support@atgcflow.com<br />
                  <strong>Platform:</strong> ATGCFLOW<br />
                  <strong>Type:</strong> Research Platform - Not for Clinical Use
                </p>
              </div>
            </section>

            {/* Final Warning */}
            <div className="bg-red-50 border-2 border-red-400 rounded-xl p-6 mt-8">
              <p className="text-red-900 font-bold text-center text-lg">
                BY USING ATGCFLOW, YOU ACKNOWLEDGE AND AGREE THAT THIS IS A RESEARCH PLATFORM NOT INTENDED FOR CLINICAL USE. YOU ASSUME ALL RISKS ASSOCIATED WITH USE OF THE PLATFORM AND AGREE NOT TO MAKE MEDICAL DECISIONS BASED ON ANALYSIS RESULTS WITHOUT PROPER VALIDATION THROUGH CERTIFIED CLINICAL LABORATORIES.
              </p>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
