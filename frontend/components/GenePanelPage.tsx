'use client';

import { useEffect, useState } from 'react';
import { panelApi, jobApi, Job, GenePanelResult } from '@/lib/api';
import { formatJobId } from '@/lib/utils';
import {
  ArrowLeft,
  Download,
  Search,
  Loader,
  AlertCircle,
  FileText,
  BarChart3,
  Filter,
  CheckCircle2,
  Database,
  Plus,
  X
} from 'lucide-react';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from '@/components/ui/table';
import { BarChart, Bar, PieChart, Pie, Cell, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';

interface GenePanelPageProps {
  jobId: string;
  onBack: () => void;
}

interface PanelSearchResult {
  id: number;
  name: string;
  disease_group: string;
  version: string;
  genes_count: number;
}

const COLORS = {
  categorical: ['#0088FE', '#00C49F', '#FFBB28', '#FF8042', '#8884D8', '#82CA9D', '#FFC658', '#FF6B9D'],
};

export default function GenePanelPage({ jobId, onBack }: GenePanelPageProps) {
  const [job, setJob] = useState<Job | null>(null);
  const [loading, setLoading] = useState(true);
  const [applying, setApplying] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Panel selection state
  const [selectedGenes, setSelectedGenes] = useState<string[]>([]);
  const [customGene, setCustomGene] = useState('');
  const [searchQuery, setSearchQuery] = useState('');
  const [searchResults, setSearchResults] = useState<PanelSearchResult[]>([]);
  const [searching, setSearching] = useState(false);

  // Results state
  const [result, setResult] = useState<GenePanelResult | null>(null);

  useEffect(() => {
    loadJob();
  }, [jobId]);

  const loadJob = async () => {
    try {
      setLoading(true);
      const jobData = await jobApi.getJob(jobId);
      setJob(jobData);
    } catch (err) {
      setError('Failed to load job details');
      console.error(err);
    } finally {
      setLoading(false);
    }
  };

  const handleSearchPanels = async () => {
    if (!searchQuery.trim()) return;

    try {
      setSearching(true);
      const results = await panelApi.searchPanels(searchQuery);
      setSearchResults(results.slice(0, 10)); // Limit to top 10
    } catch (err) {
      console.error('Panel search failed:', err);
    } finally {
      setSearching(false);
    }
  };

  const handleSelectPanel = async (panelId: number) => {
    try {
      setSearching(true);
      const panelGenes = await panelApi.getPanelGenes(panelId);
      const newGenes = panelGenes.genes.filter(g => !selectedGenes.includes(g));
      setSelectedGenes([...selectedGenes, ...newGenes]);
      setSearchResults([]);
      setSearchQuery('');
    } catch (err) {
      console.error('Failed to load panel genes:', err);
    } finally {
      setSearching(false);
    }
  };

  const handleLoadAcmgSF = async () => {
    try {
      setSearching(true);
      const acmgData = await panelApi.getAcmgSecondaryFindings();
      const newGenes = acmgData.genes.filter(g => !selectedGenes.includes(g));
      setSelectedGenes([...selectedGenes, ...newGenes]);
    } catch (err) {
      console.error('Failed to load ACMG SF genes:', err);
    } finally {
      setSearching(false);
    }
  };

  const handleAddCustomGene = () => {
    const gene = customGene.trim().toUpperCase();
    if (gene && !selectedGenes.includes(gene)) {
      setSelectedGenes([...selectedGenes, gene]);
      setCustomGene('');
    }
  };

  const handleRemoveGene = (gene: string) => {
    setSelectedGenes(selectedGenes.filter(g => g !== gene));
  };

  const handleApplyPanel = async () => {
    if (selectedGenes.length === 0) {
      setError('Please select at least one gene');
      return;
    }

    try {
      setApplying(true);
      setError(null);
      const panelResult = await panelApi.applyGenePanel(jobId, selectedGenes);
      setResult(panelResult);
    } catch (err) {
      setError('Failed to apply gene panel');
      console.error(err);
    } finally {
      setApplying(false);
    }
  };

  const handleDownload = async () => {
    if (selectedGenes.length === 0) return;

    try {
      await panelApi.downloadFilteredVariants(jobId, selectedGenes);
    } catch (err) {
      console.error('Download failed:', err);
    }
  };

  const handleClearAll = () => {
    setSelectedGenes([]);
    setResult(null);
  };

  if (loading) {
    return (
      <div className="flex items-center justify-center py-12">
        <Loader className="w-8 h-8 animate-spin text-primary" />
      </div>
    );
  }

  if (!job) {
    return (
      <div className="flex items-center gap-2 text-destructive">
        <AlertCircle className="w-5 h-5" />
        <span>Job not found</span>
      </div>
    );
  }

  // Prepare chart data
  const chrData = result ? Object.entries(result.statistics.chromosome_distribution).map(([chr, count]) => ({
    chromosome: chr,
    count
  })).sort((a, b) => {
    const aNum = parseInt(a.chromosome.replace(/[^0-9]/g, '')) || 999;
    const bNum = parseInt(b.chromosome.replace(/[^0-9]/g, '')) || 999;
    return aNum - bNum;
  }) : [];

  const funcData = result ? Object.entries(result.statistics.functional_impact).map(([name, value]) => ({
    name,
    value
  })) : [];

  const geneData = result ? Object.entries(result.statistics.gene_distribution).slice(0, 15).map(([gene, count]) => ({
    gene,
    count
  })) : [];

  return (
    <div className="space-y-6">
      {/* Header */}
      <div>
        <Button
          onClick={onBack}
          variant="ghost"
          className="mb-4"
        >
          <ArrowLeft className="w-4 h-4 mr-2" />
          Back to Job Details
        </Button>
        <div className="flex items-center justify-between">
          <div>
            <h1 className="text-3xl font-bold tracking-tight">Gene Panel Analysis</h1>
            <p className="text-sm text-muted-foreground mt-1">
              {job.sample_name} • Job ID: {formatJobId(job.job_id)}
            </p>
          </div>
          {result && (
            <Button onClick={handleDownload} className="gap-2">
              <Download className="w-4 h-4" />
              Download Filtered TSV
            </Button>
          )}
        </div>
      </div>

      {error && (
        <Card className="border-destructive">
          <CardContent className="flex items-center gap-2 pt-6 text-destructive">
            <AlertCircle className="w-5 h-5" />
            <span>{error}</span>
          </CardContent>
        </Card>
      )}

      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        {/* Left Panel - Gene Selection */}
        <div className="lg:col-span-1 space-y-4">
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <Filter className="w-5 h-5" />
                Select Genes
              </CardTitle>
              <CardDescription>Choose genes to filter variants</CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              {/* Quick Actions */}
              <div className="space-y-2">
                <Label>Quick Select</Label>
                <Button
                  onClick={handleLoadAcmgSF}
                  variant="outline"
                  size="sm"
                  className="w-full justify-start"
                  disabled={searching}
                >
                  <Database className="w-4 h-4 mr-2" />
                  ACMG Secondary Findings (81 genes)
                </Button>
              </div>

              {/* Search Panels */}
              <div className="space-y-2">
                <Label>Search PanelApp</Label>
                <div className="flex gap-2">
                  <Input
                    placeholder="Search gene panels..."
                    value={searchQuery}
                    onChange={(e) => setSearchQuery(e.target.value)}
                    onKeyPress={(e) => e.key === 'Enter' && handleSearchPanels()}
                  />
                  <Button
                    onClick={handleSearchPanels}
                    size="icon"
                    disabled={searching || !searchQuery.trim()}
                  >
                    <Search className="w-4 h-4" />
                  </Button>
                </div>
                {searching && (
                  <div className="flex items-center gap-2 text-sm text-muted-foreground">
                    <Loader className="w-4 h-4 animate-spin" />
                    Searching...
                  </div>
                )}
                {searchResults.length > 0 && (
                  <div className="border rounded-md max-h-48 overflow-y-auto">
                    {searchResults.map(panel => (
                      <button
                        key={panel.id}
                        onClick={() => handleSelectPanel(panel.id)}
                        className="w-full text-left p-2 hover:bg-accent text-sm border-b last:border-b-0"
                      >
                        <div className="font-medium">{panel.name}</div>
                        <div className="text-xs text-muted-foreground">
                          {panel.disease_group} • {panel.genes_count} genes
                        </div>
                      </button>
                    ))}
                  </div>
                )}
              </div>

              {/* Add Custom Gene */}
              <div className="space-y-2">
                <Label>Add Custom Gene</Label>
                <div className="flex gap-2">
                  <Input
                    placeholder="Gene symbol (e.g., BRCA1)"
                    value={customGene}
                    onChange={(e) => setCustomGene(e.target.value.toUpperCase())}
                    onKeyPress={(e) => e.key === 'Enter' && handleAddCustomGene()}
                  />
                  <Button
                    onClick={handleAddCustomGene}
                    size="icon"
                    disabled={!customGene.trim()}
                  >
                    <Plus className="w-4 h-4" />
                  </Button>
                </div>
              </div>

              {/* Selected Genes */}
              <div className="space-y-2">
                <div className="flex items-center justify-between">
                  <Label>Selected Genes ({selectedGenes.length})</Label>
                  {selectedGenes.length > 0 && (
                    <Button
                      onClick={handleClearAll}
                      variant="ghost"
                      size="sm"
                      className="h-auto p-1 text-xs"
                    >
                      Clear All
                    </Button>
                  )}
                </div>
                <div className="border rounded-md p-2 max-h-48 overflow-y-auto">
                  {selectedGenes.length === 0 ? (
                    <p className="text-sm text-muted-foreground text-center py-4">
                      No genes selected
                    </p>
                  ) : (
                    <div className="flex flex-wrap gap-1">
                      {selectedGenes.map(gene => (
                        <Badge key={gene} variant="secondary" className="gap-1">
                          {gene}
                          <button
                            onClick={() => handleRemoveGene(gene)}
                            className="hover:bg-destructive/20 rounded-full"
                          >
                            <X className="w-3 h-3" />
                          </button>
                        </Badge>
                      ))}
                    </div>
                  )}
                </div>
              </div>

              {/* Apply Button */}
              <Button
                onClick={handleApplyPanel}
                className="w-full"
                disabled={selectedGenes.length === 0 || applying}
              >
                {applying ? (
                  <>
                    <Loader className="w-4 h-4 mr-2 animate-spin" />
                    Applying Filter...
                  </>
                ) : (
                  <>
                    <CheckCircle2 className="w-4 h-4 mr-2" />
                    Apply Gene Panel
                  </>
                )}
              </Button>
            </CardContent>
          </Card>
        </div>

        {/* Right Panel - Results */}
        <div className="lg:col-span-2 space-y-4">
          {!result ? (
            <Card>
              <CardContent className="flex flex-col items-center justify-center py-12 text-center">
                <FileText className="w-12 h-12 text-muted-foreground mb-4" />
                <h3 className="text-lg font-semibold mb-2">No Results Yet</h3>
                <p className="text-sm text-muted-foreground max-w-md">
                  Select genes from the panel on the left and click "Apply Gene Panel" to filter variants
                  and view detailed statistics.
                </p>
              </CardContent>
            </Card>
          ) : (
            <>
              {/* Summary Cards */}
              <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                <Card>
                  <CardHeader className="pb-3">
                    <CardDescription>Total Variants</CardDescription>
                  </CardHeader>
                  <CardContent>
                    <div className="text-2xl font-bold">{result.total_variants.toLocaleString()}</div>
                  </CardContent>
                </Card>
                <Card>
                  <CardHeader className="pb-3">
                    <CardDescription>Genes Applied</CardDescription>
                  </CardHeader>
                  <CardContent>
                    <div className="text-2xl font-bold">{result.gene_count}</div>
                  </CardContent>
                </Card>
                <Card>
                  <CardHeader className="pb-3">
                    <CardDescription>Chromosomes</CardDescription>
                  </CardHeader>
                  <CardContent>
                    <div className="text-2xl font-bold">
                      {Object.keys(result.statistics.chromosome_distribution).length}
                    </div>
                  </CardContent>
                </Card>
                <Card>
                  <CardHeader className="pb-3">
                    <CardDescription>Unique Genes</CardDescription>
                  </CardHeader>
                  <CardContent>
                    <div className="text-2xl font-bold">
                      {Object.keys(result.statistics.gene_distribution).length}
                    </div>
                  </CardContent>
                </Card>
              </div>

              {/* Tabs for different views */}
              <Tabs defaultValue="charts" className="space-y-4">
                <TabsList className="grid w-full grid-cols-2">
                  <TabsTrigger value="charts">Visualizations</TabsTrigger>
                  <TabsTrigger value="variants">Variant Table</TabsTrigger>
                </TabsList>

                <TabsContent value="charts" className="space-y-4">
                  <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                    {/* Chromosome Distribution */}
                    <Card>
                      <CardHeader>
                        <CardTitle>Chromosome Distribution</CardTitle>
                        <CardDescription>Variants per chromosome</CardDescription>
                      </CardHeader>
                      <CardContent>
                        <ResponsiveContainer width="100%" height={300}>
                          <BarChart data={chrData}>
                            <CartesianGrid strokeDasharray="3 3" />
                            <XAxis dataKey="chromosome" />
                            <YAxis />
                            <Tooltip />
                            <Bar dataKey="count" fill="#0088FE" />
                          </BarChart>
                        </ResponsiveContainer>
                      </CardContent>
                    </Card>

                    {/* Functional Impact */}
                    <Card>
                      <CardHeader>
                        <CardTitle>Functional Impact</CardTitle>
                        <CardDescription>Variant categories</CardDescription>
                      </CardHeader>
                      <CardContent>
                        <ResponsiveContainer width="100%" height={300}>
                          <PieChart>
                            <Pie
                              data={funcData}
                              cx="50%"
                              cy="50%"
                              labelLine={false}
                              label={({ name, percent }) => `${name} (${percent ? (percent * 100).toFixed(0) : '0'}%)`}
                              outerRadius={80}
                              fill="#8884d8"
                              dataKey="value"
                            >
                              {funcData.map((entry, index) => (
                                <Cell key={`cell-${index}`} fill={COLORS.categorical[index % COLORS.categorical.length]} />
                              ))}
                            </Pie>
                            <Tooltip />
                          </PieChart>
                        </ResponsiveContainer>
                      </CardContent>
                    </Card>

                    {/* Gene Distribution */}
                    <Card className="md:col-span-2">
                      <CardHeader>
                        <CardTitle>Top Genes by Variant Count</CardTitle>
                        <CardDescription>Genes with most variants (top 15)</CardDescription>
                      </CardHeader>
                      <CardContent>
                        <ResponsiveContainer width="100%" height={300}>
                          <BarChart data={geneData} layout="horizontal">
                            <CartesianGrid strokeDasharray="3 3" />
                            <XAxis type="number" />
                            <YAxis dataKey="gene" type="category" width={80} />
                            <Tooltip />
                            <Bar dataKey="count" fill="#00C49F" />
                          </BarChart>
                        </ResponsiveContainer>
                      </CardContent>
                    </Card>
                  </div>
                </TabsContent>

                <TabsContent value="variants">
                  <Card>
                    <CardHeader>
                      <CardTitle>Filtered Variants</CardTitle>
                      <CardDescription>
                        Showing {result.showing_variants} of {result.total_variants} variants
                        {result.total_variants > result.showing_variants && ' (download full results)'}
                      </CardDescription>
                    </CardHeader>
                    <CardContent>
                      <div className="border rounded-md max-h-[600px] overflow-auto">
                        <Table>
                          <TableHeader>
                            <TableRow>
                              <TableHead>Position</TableHead>
                              <TableHead>Gene</TableHead>
                              <TableHead>Ref</TableHead>
                              <TableHead>Alt</TableHead>
                              <TableHead>Function</TableHead>
                              <TableHead>Consequence</TableHead>
                            </TableRow>
                          </TableHeader>
                          <TableBody>
                            {result.variants.map((variant, index) => (
                              <TableRow key={index}>
                                <TableCell className="font-mono text-xs">
                                  {variant.Chr}:{variant.Start}
                                </TableCell>
                                <TableCell className="font-medium">
                                  {variant['Gene.refGeneWithVer'] || '-'}
                                </TableCell>
                                <TableCell className="font-mono text-xs">{variant.Ref}</TableCell>
                                <TableCell className="font-mono text-xs">{variant.Alt}</TableCell>
                                <TableCell className="text-xs">
                                  {variant['Func.refGeneWithVer'] || '-'}
                                </TableCell>
                                <TableCell className="text-xs">
                                  {variant['ExonicFunc.refGeneWithVer'] || '-'}
                                </TableCell>
                              </TableRow>
                            ))}
                          </TableBody>
                        </Table>
                      </div>
                    </CardContent>
                  </Card>
                </TabsContent>
              </Tabs>
            </>
          )}
        </div>
      </div>
    </div>
  );
}
