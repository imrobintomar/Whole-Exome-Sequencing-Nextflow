'use client';

import { useEffect, useState } from 'react';
import { variantApi } from '@/lib/api';
import { ArrowLeft, Loader, AlertCircle, BarChart3, PieChart, Download, Filter, TrendingUp } from 'lucide-react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from './ui/card';
import { Button } from './ui/button';
import { Badge } from './ui/badge';
import { Tabs, TabsContent, TabsList, TabsTrigger } from './ui/tabs';
import { Switch } from './ui/switch';
import { Label } from './ui/label';
import {
  BarChart,
  Bar,
  PieChart as RechartsPieChart,
  Pie,
  Cell,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer,
} from 'recharts';

interface VariantVisualizationPageProps {
  jobId: string;
  sampleName: string;
  onBack: () => void;
}

interface VariantMetrics {
  job_id: string;
  sample_name: string;
  metrics: {
    headline: {
      total_variants: number;
      total_genes: number;
      snvs: number;
      indels: number;
      high_confidence: number;
    };
    chromosome_distribution: Array<{
      chromosome: string;
      count: number;
      type: string;
    }>;
    gene_distribution: Array<{
      gene: string;
      count: number;
    }>;
    functional_impact: {
      categories: Record<string, number>;
      exonic_subcategories: Record<string, number>;
    };
  };
}

const COLORS = {
  primary: ['#3b82f6', '#60a5fa', '#93c5fd', '#bfdbfe', '#dbeafe'],
  categorical: ['#3b82f6', '#10b981', '#f59e0b', '#ef4444', '#8b5cf6', '#ec4899', '#14b8a6', '#f97316'],
  impact: {
    high: '#ef4444',
    moderate: '#f59e0b',
    low: '#10b981',
    modifier: '#94a3b8',
  },
};

export default function VariantVisualizationPage({ jobId, sampleName, onBack }: VariantVisualizationPageProps) {
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState('');
  const [metrics, setMetrics] = useState<VariantMetrics | null>(null);
  const [normalizeChromosomes, setNormalizeChromosomes] = useState(false);
  const [proteinAlteringOnly, setProteinAlteringOnly] = useState(false);

  useEffect(() => {
    fetchMetrics();
  }, [jobId]);

  const fetchMetrics = async () => {
    try {
      setLoading(true);
      setError('');
      const data = await variantApi.getVariantMetrics(jobId);
      setMetrics(data);
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to fetch variant metrics');
    } finally {
      setLoading(false);
    }
  };

  const exportData = (data: any, filename: string) => {
    const json = JSON.stringify(data, null, 2);
    const blob = new Blob([json], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(url);
  };

  if (loading) {
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
                <p className="text-lg font-medium">Analyzing Variants</p>
                <p className="text-sm text-muted-foreground mt-1">Computing metrics for {sampleName}...</p>
              </div>
            </div>
          </CardContent>
        </Card>
      </div>
    );
  }

  if (error || !metrics) {
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
              Visualization Error
            </CardTitle>
          </CardHeader>
          <CardContent>
            <p className="text-muted-foreground mb-4">{error}</p>
            <Button onClick={fetchMetrics} variant="outline">
              Retry Loading
            </Button>
          </CardContent>
        </Card>
      </div>
    );
  }

  const { headline, chromosome_distribution, gene_distribution, functional_impact } = metrics.metrics;

  // Process chromosome data
  const chromosomeData = chromosome_distribution
    .reduce((acc: any[], curr) => {
      const existing = acc.find(item => item.chromosome === curr.chromosome);
      if (existing) {
        existing.count += curr.count;
      } else {
        acc.push({ chromosome: curr.chromosome, count: curr.count });
      }
      return acc;
    }, [])
    .sort((a, b) => {
      const chrA = a.chromosome.replace('chr', '');
      const chrB = b.chromosome.replace('chr', '');
      if (chrA === 'X') return 1;
      if (chrB === 'X') return -1;
      if (chrA === 'Y') return 1;
      if (chrB === 'Y') return -1;
      if (chrA === 'M') return 1;
      if (chrB === 'M') return -1;
      return parseInt(chrA) - parseInt(chrB);
    });

  // Process functional impact data
  const functionalImpactData = Object.entries(functional_impact.categories).map(([name, value]) => ({
    name,
    value,
  }));

  const exonicSubcategoriesData = Object.entries(functional_impact.exonic_subcategories).map(([name, value]) => ({
    name,
    value,
  }));

  // Top genes
  const topGenes = gene_distribution
    .filter(g => proteinAlteringOnly ? g.count >= 2 : true)
    .sort((a, b) => b.count - a.count)
    .slice(0, 20);

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
              <BarChart3 className="h-8 w-8 text-primary" />
              Variant Analytics
            </h1>
            <p className="text-muted-foreground mt-1">{sampleName}</p>
          </div>
          <Button onClick={() => exportData(metrics, `${sampleName}_metrics.json`)} variant="outline">
            <Download className="h-4 w-4 mr-2" />
            Export Data
          </Button>
        </div>
      </div>

      {/* Summary Cards */}
      <div className="grid gap-4 md:grid-cols-5">
        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Total Variants</CardTitle>
            <TrendingUp className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold">{headline.total_variants.toLocaleString()}</div>
            <p className="text-xs text-muted-foreground mt-1">Identified</p>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Unique Genes</CardTitle>
            <BarChart3 className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold">{headline.total_genes.toLocaleString()}</div>
            <p className="text-xs text-muted-foreground mt-1">Affected</p>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">SNVs</CardTitle>
            <Badge variant="secondary" className="text-xs">
              {((headline.snvs / headline.total_variants) * 100).toFixed(1)}%
            </Badge>
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold text-blue-500">{headline.snvs.toLocaleString()}</div>
            <p className="text-xs text-muted-foreground mt-1">Single nucleotide</p>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Indels</CardTitle>
            <Badge variant="secondary" className="text-xs">
              {((headline.indels / headline.total_variants) * 100).toFixed(1)}%
            </Badge>
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold text-green-500">{headline.indels.toLocaleString()}</div>
            <p className="text-xs text-muted-foreground mt-1">Insertions/Deletions</p>
          </CardContent>
        </Card>

        <Card className="border-yellow-500/20">
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">High Confidence</CardTitle>
            <Badge variant="secondary" className="text-xs">
              {((headline.high_confidence / headline.total_variants) * 100).toFixed(1)}%
            </Badge>
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold text-yellow-500">{headline.high_confidence.toLocaleString()}</div>
            <p className="text-xs text-muted-foreground mt-1">PASS filter</p>
          </CardContent>
        </Card>
      </div>

      {/* Visualizations */}
      <Tabs defaultValue="chromosomes" className="space-y-4">
        <TabsList className="grid w-full grid-cols-3">
          <TabsTrigger value="chromosomes">Chromosome Distribution</TabsTrigger>
          <TabsTrigger value="genes">Gene Distribution</TabsTrigger>
          <TabsTrigger value="impact">Functional Impact</TabsTrigger>
        </TabsList>

        {/* Chromosome Distribution */}
        <TabsContent value="chromosomes" className="space-y-4">
          <Card>
            <CardHeader>
              <div className="flex items-center justify-between">
                <div>
                  <CardTitle>Variants by Chromosome</CardTitle>
                  <CardDescription>Distribution across the genome</CardDescription>
                </div>
                <div className="flex items-center gap-2">
                  <Switch id="normalize-chr" checked={normalizeChromosomes} onCheckedChange={setNormalizeChromosomes} />
                  <Label htmlFor="normalize-chr" className="text-sm">
                    Normalize by size
                  </Label>
                </div>
              </div>
            </CardHeader>
            <CardContent>
              <ResponsiveContainer width="100%" height={400}>
                <BarChart data={chromosomeData}>
                  <CartesianGrid strokeDasharray="3 3" className="stroke-muted" />
                  <XAxis dataKey="chromosome" className="text-xs" />
                  <YAxis className="text-xs" />
                  <Tooltip
                    contentStyle={{ backgroundColor: 'hsl(var(--background))', border: '1px solid hsl(var(--border))' }}
                  />
                  <Legend />
                  <Bar dataKey="count" fill="hsl(var(--primary))" name="Variants" />
                </BarChart>
              </ResponsiveContainer>
            </CardContent>
          </Card>
        </TabsContent>

        {/* Gene Distribution */}
        <TabsContent value="genes" className="space-y-4">
          <Card>
            <CardHeader>
              <div className="flex items-center justify-between">
                <div>
                  <CardTitle>Top 20 Affected Genes</CardTitle>
                  <CardDescription>Genes with the most variants</CardDescription>
                </div>
                <div className="flex items-center gap-2">
                  <Switch id="protein-altering" checked={proteinAlteringOnly} onCheckedChange={setProteinAlteringOnly} />
                  <Label htmlFor="protein-altering" className="text-sm">
                    Protein-altering only
                  </Label>
                </div>
              </div>
            </CardHeader>
            <CardContent>
              <ResponsiveContainer width="100%" height={500}>
                <BarChart data={topGenes} layout="horizontal" margin={{ left: 80 }}>
                  <CartesianGrid strokeDasharray="3 3" className="stroke-muted" />
                  <XAxis type="number" className="text-xs" />
                  <YAxis dataKey="gene" type="category" className="text-xs" width={80} />
                  <Tooltip
                    contentStyle={{ backgroundColor: 'hsl(var(--background))', border: '1px solid hsl(var(--border))' }}
                  />
                  <Legend />
                  <Bar dataKey="count" fill="hsl(var(--primary))" name="Variants" />
                </BarChart>
              </ResponsiveContainer>
            </CardContent>
          </Card>
        </TabsContent>

        {/* Functional Impact */}
        <TabsContent value="impact" className="space-y-4">
          <div className="grid gap-4 md:grid-cols-2">
            <Card>
              <CardHeader>
                <CardTitle>Functional Categories</CardTitle>
                <CardDescription>Variant annotations by region</CardDescription>
              </CardHeader>
              <CardContent>
                <ResponsiveContainer width="100%" height={300}>
                  <RechartsPieChart>
                    <Pie
                      data={functionalImpactData}
                      cx="50%"
                      cy="50%"
                      labelLine={false}
                      label={({ name, percent }) => `${name} (${percent ? (percent * 100).toFixed(0) : '0'}%)`}
                      outerRadius={80}
                      fill="#8884d8"
                      dataKey="value"
                    >
                      {functionalImpactData.map((entry, index) => (
                        <Cell key={`cell-${index}`} fill={COLORS.categorical[index % COLORS.categorical.length]} />
                      ))}
                    </Pie>
                    <Tooltip />
                  </RechartsPieChart>
                </ResponsiveContainer>
              </CardContent>
            </Card>

            <Card>
              <CardHeader>
                <CardTitle>Exonic Subcategories</CardTitle>
                <CardDescription>Detailed coding region impact</CardDescription>
              </CardHeader>
              <CardContent>
                <ResponsiveContainer width="100%" height={300}>
                  <RechartsPieChart>
                    <Pie
                      data={exonicSubcategoriesData}
                      cx="50%"
                      cy="50%"
                      labelLine={false}
                      label={({ name, percent }) => `${name} (${percent ? (percent * 100).toFixed(0) : '0'}%)`}
                      outerRadius={80}
                      fill="#8884d8"
                      dataKey="value"
                    >
                      {exonicSubcategoriesData.map((entry, index) => (
                        <Cell key={`cell-${index}`} fill={COLORS.categorical[index % COLORS.categorical.length]} />
                      ))}
                    </Pie>
                    <Tooltip />
                  </RechartsPieChart>
                </ResponsiveContainer>
              </CardContent>
            </Card>
          </div>

          {/* Impact Summary Table */}
          <Card>
            <CardHeader>
              <CardTitle>Impact Summary</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="space-y-2">
                {functionalImpactData.map(({ name, value }) => (
                  <div key={name} className="flex items-center justify-between p-3 border rounded-lg">
                    <span className="font-medium">{name}</span>
                    <div className="flex items-center gap-3">
                      <div className="w-32 h-2 bg-muted rounded-full overflow-hidden">
                        <div
                          className="h-full bg-primary"
                          style={{ width: `${(value / headline.total_variants) * 100}%` }}
                        />
                      </div>
                      <span className="text-sm font-mono w-20 text-right">{value.toLocaleString()}</span>
                      <Badge variant="outline" className="w-16 justify-center">
                        {((value / headline.total_variants) * 100).toFixed(1)}%
                      </Badge>
                    </div>
                  </div>
                ))}
              </div>
            </CardContent>
          </Card>
        </TabsContent>
      </Tabs>
    </div>
  );
}
