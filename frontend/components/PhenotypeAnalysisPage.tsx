'use client';

import { Button } from './ui/button';
import { ArrowLeft } from 'lucide-react';
import PhenotypeAnalysisPanel from './PhenotypeAnalysisPanel';

interface PhenotypeAnalysisPageProps {
  jobId: string;
  sampleName: string;
  onBack: () => void;
}

export default function PhenotypeAnalysisPage({ jobId, sampleName, onBack }: PhenotypeAnalysisPageProps) {
  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-center gap-4">
        <Button onClick={onBack} variant="outline" size="sm">
          <ArrowLeft className="h-4 w-4 mr-2" />
          Back to Jobs
        </Button>
        <div>
          <h1 className="text-3xl font-bold text-gray-900">Phenotype-Driven Analysis</h1>
          <p className="text-gray-600 mt-1">
            Sample: <span className="font-semibold">{sampleName}</span>
          </p>
        </div>
      </div>

      {/* Phenotype Analysis Panel */}
      <PhenotypeAnalysisPanel jobId={jobId} sampleName={sampleName} />

      {/* Info Section */}
      <div className="bg-blue-50 border border-blue-200 rounded-lg p-6">
        <h3 className="text-lg font-semibold text-blue-900 mb-3">About Phenotype-Driven Analysis</h3>
        <div className="space-y-2 text-blue-800">
          <p>
            <strong>Exomiser Integration:</strong> This feature uses Exomiser to prioritize variants based on patient phenotype
            described using Human Phenotype Ontology (HPO) terms.
          </p>
          <p>
            <strong>How it works:</strong> Exomiser analyzes your original VCF file along with HPO terms to rank genes and variants
            by their likelihood of causing the observed phenotype.
          </p>
          <p>
            <strong>Output:</strong> Results are appended to your existing ANNOVAR annotation file with 6 additional columns:
          </p>
          <ul className="list-disc list-inside ml-4 mt-2 space-y-1">
            <li><code className="bg-white px-2 py-0.5 rounded">HPO_GENE</code> - Gene symbol</li>
            <li><code className="bg-white px-2 py-0.5 rounded">HPO_DISEASE</code> - Associated disease</li>
            <li><code className="bg-white px-2 py-0.5 rounded">HPO_SCORE</code> - Exomiser combined score</li>
            <li><code className="bg-white px-2 py-0.5 rounded">HPO_PHENO_SCORE</code> - Phenotype match score</li>
            <li><code className="bg-white px-2 py-0.5 rounded">HPO_VARIANT_SCORE</code> - Variant pathogenicity score</li>
            <li><code className="bg-white px-2 py-0.5 rounded">HPO_INHERITANCE</code> - Predicted inheritance mode</li>
          </ul>
        </div>
      </div>

      {/* Resources */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
        <div className="bg-white border border-gray-200 rounded-lg p-4">
          <h4 className="font-semibold text-gray-900 mb-2">HPO Browser</h4>
          <p className="text-sm text-gray-600 mb-3">
            Search for phenotype terms to describe patient symptoms
          </p>
          <a
            href="https://hpo.jax.org/app/"
            target="_blank"
            rel="noopener noreferrer"
            className="text-cyan hover:underline text-sm font-medium"
          >
            Browse HPO Terms →
          </a>
        </div>
        <div className="bg-white border border-gray-200 rounded-lg p-4">
          <h4 className="font-semibold text-gray-900 mb-2">Exomiser Documentation</h4>
          <p className="text-sm text-gray-600 mb-3">
            Learn more about phenotype-driven variant prioritization
          </p>
          <a
            href="https://exomiser.github.io/Exomiser/"
            target="_blank"
            rel="noopener noreferrer"
            className="text-cyan hover:underline text-sm font-medium"
          >
            View Documentation →
          </a>
        </div>
      </div>
    </div>
  );
}
