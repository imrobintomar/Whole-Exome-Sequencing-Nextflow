'use client';

import { useEffect, useState } from 'react';
import { acmgApi, Job, jobApi } from '@/lib/api';
import { formatJobId } from '@/lib/utils';
import {
  ArrowLeft,
  Download,
  Loader,
  AlertCircle,
  CheckCircle2,
  AlertTriangle,
  Info,
  Dna,
  FileText,
  Filter,
  Search,
  ChevronDown,
  ChevronUp,
} from 'lucide-react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from './ui/card';
import { Button } from './ui/button';
import { Badge } from './ui/badge';
import { Input } from './ui/input';
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from './ui/table';
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from './ui/select';

interface ACMGClassificationPageProps {
  jobId: string;
  sampleName: string;
  onBack: () => void;
}

interface VariantClassification {
  position: string;
  gene: string;
  consequence: string;
  classification: string;
  evidence_summary: string | Record<string, number>;
  met_criteria: string[];
}

interface ClassificationResult {
  job_id: string;
  sample_name: string;
  total_variants: number;
  classifications: VariantClassification[];
  summary: {
    pathogenic: number;
    likely_pathogenic: number;
    vus: number;
    likely_benign: number;
    benign: number;
  };
}

export default function ACMGClassificationPage({ jobId, sampleName, onBack }: ACMGClassificationPageProps) {
  const [loading, setLoading] = useState(true);
  const [classifying, setClassifying] = useState(false);
  const [error, setError] = useState('');
  const [result, setResult] = useState<ClassificationResult | null>(null);
  const [searchTerm, setSearchTerm] = useState('');
  const [filterClassification, setFilterClassification] = useState<string>('all');
  const [expandedRows, setExpandedRows] = useState<Set<number>>(new Set());

  useEffect(() => {
    classifyVariants();
  }, [jobId]);

  const classifyVariants = async () => {
    try {
      setClassifying(true);
      setError('');
      const data = await acmgApi.classifyJobVariants(jobId);
      setResult(data);
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to classify variants');
    } finally {
      setLoading(false);
      setClassifying(false);
    }
  };

  const getClassificationColor = (classification: string) => {
    switch (classification.toLowerCase()) {
      case 'pathogenic':
        return 'bg-red-500/10 text-red-500 border-red-500/20';
      case 'likely pathogenic':
        return 'bg-orange-500/10 text-orange-500 border-orange-500/20';
      case 'uncertain significance':
        return 'bg-yellow-500/10 text-yellow-500 border-yellow-500/20';
      case 'likely benign':
        return 'bg-blue-500/10 text-blue-500 border-blue-500/20';
      case 'benign':
        return 'bg-green-500/10 text-green-500 border-green-500/20';
      default:
        return 'bg-gray-500/10 text-gray-500 border-gray-500/20';
    }
  };

  const getClassificationIcon = (classification: string) => {
    switch (classification.toLowerCase()) {
      case 'pathogenic':
      case 'likely pathogenic':
        return <AlertCircle className="h-4 w-4" />;
      case 'uncertain significance':
        return <AlertTriangle className="h-4 w-4" />;
      case 'likely benign':
      case 'benign':
        return <CheckCircle2 className="h-4 w-4" />;
      default:
        return <Info className="h-4 w-4" />;
    }
  };

  const exportToCSV = () => {
    if (!result) return;

    const headers = ['UniqueID', 'Gene', 'Consequence', 'Classification', 'Evidence Summary', 'Met Criteria'];
    const rows = filteredVariants.map(v => [
      v.position,
      v.gene,
      v.consequence,
      v.classification,
      typeof v.evidence_summary === 'string' ? v.evidence_summary : JSON.stringify(v.evidence_summary),
      v.met_criteria.join('; '),
    ]);

    const csv = [headers, ...rows].map(row => row.map(cell => `"${cell}"`).join(',')).join('\n');
    const blob = new Blob([csv], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${sampleName}_ACMG_classification.csv`;
    a.click();
    URL.revokeObjectURL(url);
  };

  const toggleRow = (index: number) => {
    setExpandedRows(prev => {
      const next = new Set(prev);
      if (next.has(index)) {
        next.delete(index);
      } else {
        next.add(index);
      }
      return next;
    });
  };

  if (loading || classifying) {
    return (
      <div className="space-y-6">
        <Button onClick={onBack} variant="ghost" size="sm">
          <ArrowLeft className="h-4 w-4 mr-2" />
          Back
        </Button>
        <Card>
          <CardContent className="flex items-center justify-center h-96">
            <div className="text-center space-y-4">
              <Loader className="h-12 w-12 animate-spin text-primary mx-auto" />
              <div>
                <p className="text-lg font-medium">Classifying Variants</p>
                <p className="text-sm text-muted-foreground mt-1">
                  Applying ACMG/AMP 2015 guidelines to {sampleName}...
                </p>
              </div>
            </div>
          </CardContent>
        </Card>
      </div>
    );
  }

  if (error || !result) {
    return (
      <div className="space-y-6">
        <Button onClick={onBack} variant="ghost" size="sm">
          <ArrowLeft className="h-4 w-4 mr-2" />
          Back
        </Button>
        <Card className="border-red-500/50">
          <CardHeader>
            <CardTitle className="flex items-center gap-2 text-red-500">
              <AlertCircle className="h-5 w-5" />
              Classification Error
            </CardTitle>
          </CardHeader>
          <CardContent>
            <p className="text-muted-foreground mb-4">{error}</p>
            <Button onClick={classifyVariants} variant="outline">
              Retry Classification
            </Button>
          </CardContent>
        </Card>
      </div>
    );
  }

  const filteredVariants = result.classifications.filter(v => {
    const matchesSearch =
      searchTerm === '' ||
      v.gene.toLowerCase().includes(searchTerm.toLowerCase()) ||
      v.position.toLowerCase().includes(searchTerm.toLowerCase()) ||
      v.consequence.toLowerCase().includes(searchTerm.toLowerCase());

    const matchesFilter =
      filterClassification === 'all' || v.classification.toLowerCase() === filterClassification.toLowerCase();

    return matchesSearch && matchesFilter;
  });

  return (
    <div className="space-y-6">
      {/* Header */}
      <div>
        <Button onClick={onBack} variant="ghost" size="sm" className="mb-2">
          <ArrowLeft className="h-4 w-4 mr-2" />
          Back to Job Details
        </Button>
        <div className="flex items-start justify-between">
          <div>
            <h1 className="text-3xl font-bold tracking-tight flex items-center gap-3">
              <Dna className="h-8 w-8 text-primary" />
              ACMG Variant Classification
            </h1>
            <p className="text-muted-foreground mt-1">{sampleName}</p>
          </div>
          <Button onClick={exportToCSV} variant="outline">
            <Download className="h-4 w-4 mr-2" />
            Export CSV
          </Button>
        </div>
      </div>

      {/* Summary Cards */}
      <div className="grid gap-4 md:grid-cols-5">
        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Total Variants</CardTitle>
            <FileText className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold">{result.total_variants}</div>
            <p className="text-xs text-muted-foreground mt-1">Analyzed</p>
          </CardContent>
        </Card>

        <Card className="border-red-500/20">
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Pathogenic</CardTitle>
            <AlertCircle className="h-4 w-4 text-red-500" />
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold text-red-500">
              {result.summary.pathogenic + result.summary.likely_pathogenic}
            </div>
            <p className="text-xs text-muted-foreground mt-1">
              {result.summary.pathogenic} P / {result.summary.likely_pathogenic} LP
            </p>
          </CardContent>
        </Card>

        <Card className="border-yellow-500/20">
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">VUS</CardTitle>
            <AlertTriangle className="h-4 w-4 text-yellow-500" />
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold text-yellow-500">{result.summary.vus}</div>
            <p className="text-xs text-muted-foreground mt-1">Uncertain</p>
          </CardContent>
        </Card>

        <Card className="border-green-500/20">
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Benign</CardTitle>
            <CheckCircle2 className="h-4 w-4 text-green-500" />
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold text-green-500">
              {result.summary.benign + result.summary.likely_benign}
            </div>
            <p className="text-xs text-muted-foreground mt-1">
              {result.summary.likely_benign} LB / {result.summary.benign} B
            </p>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Actionable</CardTitle>
            <Dna className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold">
              {result.summary.pathogenic + result.summary.likely_pathogenic}
            </div>
            <p className="text-xs text-muted-foreground mt-1">
              {(((result.summary.pathogenic + result.summary.likely_pathogenic) / result.total_variants) * 100).toFixed(1)}%
            </p>
          </CardContent>
        </Card>
      </div>

      {/* ACMG Guidelines Info */}
      <Card className="bg-blue-500/5 border-blue-500/20">
        <CardHeader>
          <CardTitle className="flex items-center gap-2 text-blue-500">
            <Info className="h-5 w-5" />
            ACMG/AMP 2015 Guidelines
          </CardTitle>
        </CardHeader>
        <CardContent className="text-sm text-muted-foreground space-y-2">
          <p>
            Variants are classified according to the{' '}
            <strong>American College of Medical Genetics and Genomics (ACMG)</strong> 2015 standards.
          </p>
          <p>
            Evidence codes: <strong>PVS</strong> (very strong), <strong>PS</strong> (strong), <strong>PM</strong> (moderate),{' '}
            <strong>PP</strong> (supporting), <strong>BA/BS/BP</strong> (benign evidence)
          </p>
        </CardContent>
      </Card>

      {/* Variant Table with Filters */}
      <Card>
        <CardHeader>
          <CardTitle>Variant Table</CardTitle>
          <CardDescription>
            Classified variants with evidence summary
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="flex gap-4">
            <div className="relative flex-1">
              <Search className="absolute left-3 top-1/2 -translate-y-1/2 h-4 w-4 text-muted-foreground" />
              <Input
                placeholder="Search by gene, position, or consequence..."
                value={searchTerm}
                onChange={e => setSearchTerm(e.target.value)}
                className="pl-9"
              />
            </div>
            <Select value={filterClassification} onValueChange={setFilterClassification}>
              <SelectTrigger className="w-[200px]">
                <Filter className="h-4 w-4 mr-2" />
                <SelectValue placeholder="Filter by classification" />
              </SelectTrigger>
              <SelectContent>
                <SelectItem value="all">All Classifications</SelectItem>
                <SelectItem value="pathogenic">Pathogenic</SelectItem>
                <SelectItem value="likely pathogenic">Likely Pathogenic</SelectItem>
                <SelectItem value="uncertain significance">VUS</SelectItem>
                <SelectItem value="likely benign">Likely Benign</SelectItem>
                <SelectItem value="benign">Benign</SelectItem>
              </SelectContent>
            </Select>
          </div>

          <div className="text-sm text-muted-foreground">
            Showing {filteredVariants.length} of {result.total_variants} variants
          </div>

          {/* Variants Table */}
          <div className="border rounded-lg">
            <Table>
              <TableHeader>
                <TableRow>
                  <TableHead className="w-[50px]"></TableHead>
                  <TableHead>UniqueID</TableHead>
                  <TableHead>Gene</TableHead>
                  <TableHead>Consequence</TableHead>
                  <TableHead>Classification</TableHead>
                  <TableHead className="w-[100px]">Evidence</TableHead>
                </TableRow>
              </TableHeader>
              <TableBody>
                {filteredVariants.length === 0 ? (
                  <TableRow>
                    <TableCell colSpan={6} className="text-center py-8 text-muted-foreground">
                      No variants found matching your filters
                    </TableCell>
                  </TableRow>
                ) : (
                  filteredVariants.map((variant, index) => (
                    <>
                      <TableRow key={index} className="cursor-pointer hover:bg-accent" onClick={() => toggleRow(index)}>
                        <TableCell>
                          {expandedRows.has(index) ? (
                            <ChevronUp className="h-4 w-4 text-muted-foreground" />
                          ) : (
                            <ChevronDown className="h-4 w-4 text-muted-foreground" />
                          )}
                        </TableCell>
                        <TableCell className="font-mono text-sm">{variant.position}</TableCell>
                        <TableCell className="font-medium">{variant.gene || 'N/A'}</TableCell>
                        <TableCell>{variant.consequence || 'N/A'}</TableCell>
                        <TableCell>
                          <Badge className={`${getClassificationColor(variant.classification)} flex items-center gap-1.5 w-fit`}>
                            {getClassificationIcon(variant.classification)}
                            <span>{variant.classification}</span>
                          </Badge>
                        </TableCell>
                        <TableCell>
                          <div className="flex flex-wrap gap-1">
                            {variant.met_criteria.slice(0, 3).map(criteria => (
                              <Badge key={criteria} variant="outline" className="text-xs">
                                {criteria}
                              </Badge>
                            ))}
                            {variant.met_criteria.length > 3 && (
                              <Badge variant="outline" className="text-xs">
                                +{variant.met_criteria.length - 3}
                              </Badge>
                            )}
                          </div>
                        </TableCell>
                      </TableRow>
                      {expandedRows.has(index) && (
                        <TableRow>
                          <TableCell colSpan={6} className="bg-muted/50">
                            <div className="p-4 space-y-3">
                              <div>
                                <p className="text-sm font-medium mb-1">Evidence Summary</p>
                                <p className="text-sm text-muted-foreground">
                                  {typeof variant.evidence_summary === 'string'
                                    ? variant.evidence_summary || 'No evidence summary available'
                                    : JSON.stringify(variant.evidence_summary, null, 2) || 'No evidence summary available'}
                                </p>
                              </div>
                              <div>
                                <p className="text-sm font-medium mb-2">Met ACMG Criteria ({variant.met_criteria.length})</p>
                                <div className="flex flex-wrap gap-2">
                                  {variant.met_criteria.map(criteria => (
                                    <Badge key={criteria} variant="secondary">
                                      {criteria}
                                    </Badge>
                                  ))}
                                </div>
                              </div>
                            </div>
                          </TableCell>
                        </TableRow>
                      )}
                    </>
                  ))
                )}
              </TableBody>
            </Table>
          </div>
        </CardContent>
      </Card>
    </div>
  );
}
