'use client';

import { useState, useEffect } from 'react';
import { Card, CardHeader, CardTitle, CardContent } from './ui/card';
import { Button } from './ui/button';
import { Input } from './ui/input';
import { Label } from './ui/label';
import { phenotypeApi, PhenotypeStatusResponse } from '../lib/api_phenotype';
import { Badge } from './ui/badge';

interface PhenotypeAnalysisPanelProps {
  jobId: string;
  sampleName: string;
}

export default function PhenotypeAnalysisPanel({ jobId, sampleName }: PhenotypeAnalysisPanelProps) {
  const [hpoInput, setHpoInput] = useState('');
  const [hpoTerms, setHpoTerms] = useState<string[]>([]);
  const [status, setStatus] = useState<PhenotypeStatusResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [running, setRunning] = useState(false);

  useEffect(() => {
    loadStatus();
    const interval = setInterval(loadStatus, 5000);
    return () => clearInterval(interval);
  }, [jobId]);

  const loadStatus = async () => {
    try {
      const data = await phenotypeApi.getPhenotypeStatus(jobId);
      setStatus(data);
      if (data.has_phenotype_analysis) {
        setRunning(false);
      }
    } catch (err: any) {
      // Silently handle 404 errors (job not found or no phenotype analysis yet)
      // Only log other errors
      if (err.response?.status !== 404) {
        console.error('Failed to load phenotype status:', err);
      }
    }
  };

  const addHpoTerm = () => {
    const term = hpoInput.trim().toUpperCase();
    if (!term) return;

    const hpoPattern = /^HP:\d{7}$/;
    if (!hpoPattern.test(term)) {
      setError('Invalid HPO term format. Use HP:0000000');
      return;
    }

    if (hpoTerms.includes(term)) {
      setError('HPO term already added');
      return;
    }

    setHpoTerms([...hpoTerms, term]);
    setHpoInput('');
    setError(null);
  };

  const removeTerm = (term: string) => {
    setHpoTerms(hpoTerms.filter(t => t !== term));
  };

  const runAnalysis = async () => {
    if (hpoTerms.length === 0) {
      setError('Add at least one HPO term');
      return;
    }

    try {
      setLoading(true);
      setRunning(true);
      setError(null);
      await phenotypeApi.runPhenotypeAnalysis(jobId, hpoTerms);
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to start analysis');
      setRunning(false);
    } finally {
      setLoading(false);
    }
  };

  const downloadResults = async () => {
    try {
      await phenotypeApi.downloadPhenotypeFile(jobId);
    } catch (err: any) {
      setError('Failed to download results');
    }
  };

  return (
    <Card className="border-2 border-purple-primary/20">
      <CardHeader className="bg-gradient-to-r from-purple-primary/5 to-cyan/5">
        <CardTitle className="flex items-center gap-2">
          <span className="text-2xl">🧬</span>
          Phenotype-Driven Analysis
        </CardTitle>
      </CardHeader>
      <CardContent className="space-y-6 pt-6">
        {/* Status Display */}
        {status?.has_phenotype_analysis && (
          <div className="p-4 bg-green-50 border border-green-200 rounded-lg">
            <div className="flex items-center justify-between">
              <div>
                <p className="font-semibold text-green-800">Analysis Complete</p>
                <p className="text-sm text-green-600">
                  {status.hpo_terms.length} HPO terms analyzed
                </p>
              </div>
              <Button onClick={downloadResults} variant="outline" size="sm">
                Download Results
              </Button>
            </div>
          </div>
        )}

        {running && !status?.has_phenotype_analysis && (
          <div className="p-4 bg-blue-50 border border-blue-200 rounded-lg">
            <div className="flex items-center gap-2">
              <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-blue-600"></div>
              <p className="text-blue-800">Running Exomiser analysis...</p>
            </div>
          </div>
        )}

        {/* HPO Term Input */}
        <div className="space-y-4">
          <div>
            <Label htmlFor="hpo-input">HPO Terms</Label>
            <div className="flex gap-2 mt-1">
              <Input
                id="hpo-input"
                placeholder="HP:0001250 (e.g., Seizure)"
                value={hpoInput}
                onChange={(e) => setHpoInput(e.target.value)}
                onKeyPress={(e) => e.key === 'Enter' && addHpoTerm()}
                disabled={running}
              />
              <Button onClick={addHpoTerm} disabled={running}>Add</Button>
            </div>
            <p className="text-xs text-gray-600 mt-1">
              Enter HPO terms in format HP:0000000. Search terms at{' '}
              <a
                href="https://hpo.jax.org/app/"
                target="_blank"
                rel="noopener noreferrer"
                className="text-cyan hover:underline"
              >
                hpo.jax.org
              </a>
            </p>
          </div>

          {/* HPO Terms List */}
          {hpoTerms.length > 0 && (
            <div className="flex flex-wrap gap-2">
              {hpoTerms.map((term) => (
                <Badge key={term} variant="secondary" className="flex items-center gap-1">
                  {term}
                  {!running && (
                    <button
                      onClick={() => removeTerm(term)}
                      className="ml-1 hover:text-red-600"
                    >
                      ✕
                    </button>
                  )}
                </Badge>
              ))}
            </div>
          )}
        </div>

        {/* Error Display */}
        {error && (
          <div className="p-3 bg-red-50 border border-red-200 rounded text-red-800 text-sm">
            {error}
          </div>
        )}

        {/* Action Buttons */}
        <div className="flex gap-2">
          <Button
            onClick={runAnalysis}
            disabled={loading || running || hpoTerms.length === 0}
            className="flex-1"
          >
            {loading ? 'Starting...' : 'Run Phenotype Analysis'}
          </Button>
          {status?.has_phenotype_analysis && (
            <Button onClick={downloadResults} variant="outline">
              Download
            </Button>
          )}
        </div>

        {/* Info Box */}
        <div className="p-4 bg-purple-50 border border-purple-200 rounded-lg text-sm">
          <p className="font-semibold text-purple-900 mb-2">How it works:</p>
          <ul className="space-y-1 text-purple-800 list-disc list-inside">
            <li>Enter HPO terms describing patient phenotype</li>
            <li>Exomiser runs on original VCF with phenotype prioritization</li>
            <li>Results appended to ANNOVAR file: HPO_GENE, HPO_DISEASE, HPO_SCORE</li>
            <li>Download augmented file with phenotype-driven variant ranking</li>
          </ul>
        </div>
      </CardContent>
    </Card>
  );
}
