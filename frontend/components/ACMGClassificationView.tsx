'use client';

import { useState } from 'react';
import { acmgApi, JobClassificationResult } from '@/lib/api';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from './ui/card';
import { Button } from './ui/button';
import { Badge } from './ui/badge';
import ACMGBadge from './ACMGBadge';
import ACMGEvidenceMatrix from './ACMGEvidenceMatrix';
import { Loader, Download, PieChart, List } from 'lucide-react';
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from './ui/table';
import {
  Dialog,
  DialogContent,
  DialogDescription,
  DialogHeader,
  DialogTitle,
} from './ui/dialog';

interface ACMGClassificationViewProps {
  jobId: string;
  sampleName: string;
}

export default function ACMGClassificationView({ jobId, sampleName }: ACMGClassificationViewProps) {
  const [result, setResult] = useState<JobClassificationResult | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [selectedVariant, setSelectedVariant] = useState<any>(null);
  const [viewMode, setViewMode] = useState<'summary' | 'table'>('summary');

  const handleClassify = async () => {
    setLoading(true);
    setError('');
    try {
      const classificationResult = await acmgApi.classifyJobVariants(jobId);
      setResult(classificationResult);
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Classification failed');
      console.error('Classification error:', err);
    } finally {
      setLoading(false);
    }
  };

  const exportResults = () => {
    if (!result) return;

    const csv = [
      ['Position', 'Gene', 'Consequence', 'Classification', 'Met Criteria'].join(','),
      ...result.classifications.map(c =>
        [c.position, c.gene, c.consequence, c.classification, c.met_criteria.join(';')].join(',')
      )
    ].join('\n');

    const blob = new Blob([csv], { type: 'text/csv' });
    const url = window.URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = `${sampleName}_acmg_classifications.csv`;
    link.click();
    window.URL.revokeObjectURL(url);
  };

  if (!result) {
    return (
      <Card>
        <CardHeader>
          <CardTitle>ACMG/AMP Classification</CardTitle>
          <CardDescription>
            Classify Variants according to 2015 ACMG/AMP guidelines
          </CardDescription>
        </CardHeader>
        <CardContent>
          {error && (
            <div className="mb-4 p-3 bg-destructive/10 text-destructive rounded-md text-sm">
              {error}
            </div>
          )}
          <Button onClick={handleClassify} disabled={loading} className="w-full">
            {loading ? (
              <>
                <Loader className="mr-2 h-4 w-4 animate-spin" />
                Classifying variants...
              </>
            ) : (
              'Run ACMG Classification'
            )}
          </Button>
          <p className="text-xs text-muted-foreground mt-2">
            This will analyze all variants in the TSV file and apply ACMG criteria
          </p>
        </CardContent>
      </Card>
    );
  }

  return (
    <>
      <Card>
        <CardHeader>
          <div className="flex items-center justify-between">
            <div>
              <CardTitle>ACMG Classification Results</CardTitle>
              <CardDescription>
                {result.total_variants} variants classified for {result.sample_name}
              </CardDescription>
            </div>
            <div className="flex gap-2">
              <Button
                variant="outline"
                size="sm"
                onClick={() => setViewMode(viewMode === 'summary' ? 'table' : 'summary')}
              >
                {viewMode === 'summary' ? <List className="h-4 w-4" /> : <PieChart className="h-4 w-4" />}
              </Button>
              <Button variant="outline" size="sm" onClick={exportResults}>
                <Download className="mr-2 h-4 w-4" />
                Export CSV
              </Button>
              <Button variant="outline" size="sm" onClick={handleClassify} disabled={loading}>
                Re-classify
              </Button>
            </div>
          </div>
        </CardHeader>
        <CardContent>
          {viewMode === 'summary' ? (
            <div className="space-y-6">
              {/* Summary Cards */}
              <div className="grid grid-cols-5 gap-4">
                <Card className="border-red-200 bg-red-50">
                  <CardContent className="pt-6">
                    <div className="text-center">
                      <div className="text-3xl font-bold text-red-600">
                        {result.summary.pathogenic}
                      </div>
                      <ACMGBadge classification="Pathogenic" size="sm" className="mt-2" />
                    </div>
                  </CardContent>
                </Card>

                <Card className="border-orange-200 bg-orange-50">
                  <CardContent className="pt-6">
                    <div className="text-center">
                      <div className="text-3xl font-bold text-orange-600">
                        {result.summary.likely_pathogenic}
                      </div>
                      <ACMGBadge classification="Likely Pathogenic" size="sm" className="mt-2" />
                    </div>
                  </CardContent>
                </Card>

                <Card className="border-yellow-200 bg-yellow-50">
                  <CardContent className="pt-6">
                    <div className="text-center">
                      <div className="text-3xl font-bold text-yellow-600">
                        {result.summary.vus}
                      </div>
                      <ACMGBadge classification="Uncertain Significance" size="sm" className="mt-2" />
                    </div>
                  </CardContent>
                </Card>

                <Card className="border-green-200 bg-green-50">
                  <CardContent className="pt-6">
                    <div className="text-center">
                      <div className="text-3xl font-bold text-green-600">
                        {result.summary.likely_benign}
                      </div>
                      <ACMGBadge classification="Likely Benign" size="sm" className="mt-2" />
                    </div>
                  </CardContent>
                </Card>

                <Card className="border-green-300 bg-green-100">
                  <CardContent className="pt-6">
                    <div className="text-center">
                      <div className="text-3xl font-bold text-green-700">
                        {result.summary.benign}
                      </div>
                      <ACMGBadge classification="Benign" size="sm" className="mt-2" />
                    </div>
                  </CardContent>
                </Card>
              </div>

              {/* P/LP Variants Table */}
              {(result.summary.pathogenic + result.summary.likely_pathogenic) > 0 && (
                <div>
                  <h3 className="text-lg font-semibold mb-3">
                    Pathogenic & Likely Pathogenic Variants
                  </h3>
                  <Table>
                    <TableHeader>
                      <TableRow>
                        <TableHead>Position</TableHead>
                        <TableHead>Gene</TableHead>
                        <TableHead>Consequence</TableHead>
                        <TableHead>Classification</TableHead>
                        <TableHead>Evidence</TableHead>
                        <TableHead>Actions</TableHead>
                      </TableRow>
                    </TableHeader>
                    <TableBody>
                      {result.classifications
                        .filter(c => c.classification === 'Pathogenic' || c.classification === 'Likely Pathogenic')
                        .map((variant, idx) => (
                          <TableRow key={idx}>
                            <TableCell className="font-mono text-sm">{variant.position}</TableCell>
                            <TableCell className="font-semibold">{variant.gene}</TableCell>
                            <TableCell className="text-sm">{variant.consequence}</TableCell>
                            <TableCell>
                              <ACMGBadge
                                classification={variant.classification as any}
                                size="sm"
                              />
                            </TableCell>
                            <TableCell>
                              <div className="flex gap-1 flex-wrap">
                                {variant.met_criteria.slice(0, 3).map(c => (
                                  <Badge key={c} variant="outline" className="text-xs">
                                    {c}
                                  </Badge>
                                ))}
                                {variant.met_criteria.length > 3 && (
                                  <Badge variant="secondary" className="text-xs">
                                    +{variant.met_criteria.length - 3}
                                  </Badge>
                                )}
                              </div>
                            </TableCell>
                            <TableCell>
                              <Button
                                variant="ghost"
                                size="sm"
                                onClick={() => setSelectedVariant(variant)}
                              >
                                Details
                              </Button>
                            </TableCell>
                          </TableRow>
                        ))}
                    </TableBody>
                  </Table>
                </div>
              )}
            </div>
          ) : (
            /* Full Table View */
            <Table>
              <TableHeader>
                <TableRow>
                  <TableHead>Variant</TableHead>
                  <TableHead>Gene</TableHead>
                  <TableHead>Consequence</TableHead>
                  <TableHead>Classification</TableHead>
                  <TableHead>Evidence Count</TableHead>
                  <TableHead>Actions</TableHead>
                </TableRow>
              </TableHeader>
              <TableBody>
                {result.classifications.map((variant, idx) => (
                  <TableRow key={idx}>
                    <TableCell className="font-mono text-sm">{variant.position}</TableCell>
                    <TableCell className="font-semibold">{variant.gene}</TableCell>
                    <TableCell className="text-sm">{variant.consequence}</TableCell>
                    <TableCell>
                      <ACMGBadge classification={variant.classification as any} size="sm" />
                    </TableCell>
                    <TableCell>
                      <div className="text-sm">
                        {Object.entries(variant.evidence_summary)
                          .filter(([_, count]) => count > 0)
                          .map(([key, count]) => `${key}:${count}`)
                          .join(', ') || 'None'}
                      </div>
                    </TableCell>
                    <TableCell>
                      <Button
                        variant="ghost"
                        size="sm"
                        onClick={() => setSelectedVariant(variant)}
                      >
                        Details
                      </Button>
                    </TableCell>
                  </TableRow>
                ))}
              </TableBody>
            </Table>
          )}
        </CardContent>
      </Card>

      {/* Variant Detail Dialog */}
      <Dialog open={selectedVariant !== null} onOpenChange={() => setSelectedVariant(null)}>
        <DialogContent className="max-w-4xl max-h-[90vh] overflow-y-auto">
          <DialogHeader>
            <DialogTitle>
              {selectedVariant?.gene} - {selectedVariant?.position}
            </DialogTitle>
            <DialogDescription>
              {selectedVariant?.consequence}
            </DialogDescription>
          </DialogHeader>
          {selectedVariant && (
            <div className="space-y-4">
              <div className="flex items-center gap-4">
                <div>
                  <div className="text-sm text-muted-foreground">Classification</div>
                  <ACMGBadge classification={selectedVariant.classification} size="lg" />
                </div>
                <div>
                  <div className="text-sm text-muted-foreground">Met Criteria</div>
                  <div className="flex gap-1 flex-wrap mt-1">
                    {selectedVariant.met_criteria.map((c: string) => (
                      <Badge key={c}>{c}</Badge>
                    ))}
                  </div>
                </div>
              </div>

              <div>
                <h4 className="font-semibold mb-2">Evidence Summary</h4>
                <div className="grid grid-cols-7 gap-2">
                  {Object.entries(selectedVariant.evidence_summary).map(([key, count]) => (
                    <div key={key} className="text-center p-2 bg-muted rounded">
                      <div className="text-2xl font-bold">{count as number}</div>
                      <div className="text-xs text-muted-foreground">{key}</div>
                    </div>
                  ))}
                </div>
              </div>
            </div>
          )}
        </DialogContent>
      </Dialog>
    </>
  );
}
