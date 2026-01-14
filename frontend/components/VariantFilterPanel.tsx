import React, { useState } from 'react';
import { Card, CardContent, CardHeader, CardTitle } from './ui/card';
import { Button } from './ui/button';
import { Checkbox } from './ui/checkbox';
import { Slider } from './ui/slider';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from './ui/select';
import { Label } from './ui/label';
import { Input } from './ui/input';
import { Loader2, Download, BarChart3, FileText, ChevronDown, ChevronUp } from 'lucide-react';
import axios from 'axios';
import { auth } from '@/lib/firebase';
import { Alert, AlertDescription } from './ui/alert';

const API_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

// Helper function to get auth headers
async function getAuthHeaders() {
  const headers: any = { 'Content-Type': 'application/json' };
  try {
    const user = auth.currentUser;
    if (user) {
      const token = await user.getIdToken(true);
      headers.Authorization = `Bearer ${token}`;
    }
  } catch (error) {
    console.error('Error getting auth token:', error);
  }
  return headers;
}

interface FilterConfig {
  include_exonic: boolean;
  include_splicing: boolean;
  exclude_synonymous: boolean;
  exclude_benign: boolean;
  exclude_likely_benign: boolean;
  max_gnomad_af: number;
  gnomad_af_column: string;
  min_depth: number;
  require_pass: boolean;
}

interface FilterStatistics {
  total_variants: number;
  filtered_variants: number;
  removed_variants: number;
  reduction_percentage: number;
  top_genes: Array<{ gene: string; count: number }>;
  functional_distribution: Record<string, number>;
}

interface VariantFilterPanelProps {
  jobId: string;
  onClassifyResults?: (results: any) => void;
  onVisualizeResults?: (metrics: any) => void;
}

export function VariantFilterPanel({ jobId, onClassifyResults, onVisualizeResults }: VariantFilterPanelProps) {
  const [filterConfig, setFilterConfig] = useState<FilterConfig>({
    include_exonic: true,
    include_splicing: true,
    exclude_synonymous: true,
    exclude_benign: true,
    exclude_likely_benign: true,
    max_gnomad_af: 0.05,
    gnomad_af_column: 'gnomad40_exome_AF',
    min_depth: 5,
    require_pass: false,
  });

  const [previewStats, setPreviewStats] = useState<FilterStatistics | null>(null);
  const [loading, setLoading] = useState<{
    preview: boolean;
    download: boolean;
    classify: boolean;
    visualize: boolean;
  }>({
    preview: false,
    download: false,
    classify: false,
    visualize: false,
  });
  const [error, setError] = useState<string | null>(null);
  const [expandedSections, setExpandedSections] = useState({
    functional: true,
    clinical: true,
    population: true,
    quality: true,
  });

  const toggleSection = (section: keyof typeof expandedSections) => {
    setExpandedSections(prev => ({ ...prev, [section]: !prev[section] }));
  };

  const handlePreview = async () => {
    setLoading(prev => ({ ...prev, preview: true }));
    setError(null);

    try {
      const headers = await getAuthHeaders();
      const response = await axios.post(
        `${API_URL}/jobs/${jobId}/filter/preview`,
        filterConfig,
        { headers }
      );
      setPreviewStats(response.data.statistics);
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to preview filters');
    } finally {
      setLoading(prev => ({ ...prev, preview: false }));
    }
  };

  const handleDownload = async () => {
    setLoading(prev => ({ ...prev, download: true }));
    setError(null);

    try {
      const headers = await getAuthHeaders();
      const response = await axios.post(
        `${API_URL}/jobs/${jobId}/filter/download`,
        filterConfig,
        {
          headers,
          responseType: 'blob',
        }
      );

      // Create download link
      const url = window.URL.createObjectURL(new Blob([response.data]));
      const link = document.createElement('a');
      link.href = url;
      const filename = response.headers['content-disposition']?.split('filename=')[1] || 'filtered_variants.csv';
      link.setAttribute('download', filename);
      document.body.appendChild(link);
      link.click();
      link.remove();
      window.URL.revokeObjectURL(url);
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to download filtered variants');
    } finally {
      setLoading(prev => ({ ...prev, download: false }));
    }
  };

  const handleClassify = async () => {
    setLoading(prev => ({ ...prev, classify: true }));
    setError(null);

    try {
      const headers = await getAuthHeaders();
      const response = await axios.post(
        `${API_URL}/jobs/${jobId}/filter/classify`,
        filterConfig,
        { headers }
      );
      if (onClassifyResults) {
        onClassifyResults(response.data);
      }
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to classify filtered variants');
    } finally {
      setLoading(prev => ({ ...prev, classify: false }));
    }
  };

  const handleVisualize = async () => {
    setLoading(prev => ({ ...prev, visualize: true }));
    setError(null);

    try {
      const headers = await getAuthHeaders();
      const response = await axios.post(
        `${API_URL}/jobs/${jobId}/filter/visualize`,
        filterConfig,
        { headers }
      );
      if (onVisualizeResults) {
        onVisualizeResults(response.data);
      }
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to generate visualizations');
    } finally {
      setLoading(prev => ({ ...prev, visualize: false }));
    }
  };

  return (
    <div className="space-y-4">
      <Card>
        <CardHeader>
          <CardTitle>Variant Filtering</CardTitle>
          <p className="text-sm text-muted-foreground">
            Apply custom filters to annotated variants. Output: CSV file (comma-delimited, 287 columns)
          </p>
        </CardHeader>
        <CardContent className="space-y-6">
          {error && (
            <Alert variant="destructive">
              <AlertDescription>{error}</AlertDescription>
            </Alert>
          )}

          {/* Functional Impact Filters */}
          <div className="space-y-3">
            <div
              className="flex items-center justify-between cursor-pointer"
              onClick={() => toggleSection('functional')}
            >
              <h3 className="text-lg font-semibold flex items-center gap-2">
                📊 Functional Impact Filters
              </h3>
              {expandedSections.functional ? <ChevronUp size={20} /> : <ChevronDown size={20} />}
            </div>
            {expandedSections.functional && (
              <div className="space-y-2 pl-6">
                <div className="flex items-center space-x-2">
                  <Checkbox
                    id="include_exonic"
                    checked={filterConfig.include_exonic}
                    onCheckedChange={(checked) =>
                      setFilterConfig(prev => ({ ...prev, include_exonic: checked as boolean }))
                    }
                  />
                  <Label htmlFor="include_exonic">Include exonic variants</Label>
                </div>
                <div className="flex items-center space-x-2">
                  <Checkbox
                    id="include_splicing"
                    checked={filterConfig.include_splicing}
                    onCheckedChange={(checked) =>
                      setFilterConfig(prev => ({ ...prev, include_splicing: checked as boolean }))
                    }
                  />
                  <Label htmlFor="include_splicing">Include splicing variants</Label>
                </div>
                <div className="flex items-center space-x-2">
                  <Checkbox
                    id="exclude_synonymous"
                    checked={filterConfig.exclude_synonymous}
                    onCheckedChange={(checked) =>
                      setFilterConfig(prev => ({ ...prev, exclude_synonymous: checked as boolean }))
                    }
                  />
                  <Label htmlFor="exclude_synonymous">Exclude synonymous SNVs</Label>
                </div>
              </div>
            )}
          </div>

          {/* Clinical Significance Filters */}
          <div className="space-y-3">
            <div
              className="flex items-center justify-between cursor-pointer"
              onClick={() => toggleSection('clinical')}
            >
              <h3 className="text-lg font-semibold flex items-center gap-2">
                🏥 Clinical Significance Filters
              </h3>
              {expandedSections.clinical ? <ChevronUp size={20} /> : <ChevronDown size={20} />}
            </div>
            {expandedSections.clinical && (
              <div className="space-y-2 pl-6">
                <div className="flex items-center space-x-2">
                  <Checkbox
                    id="exclude_benign"
                    checked={filterConfig.exclude_benign}
                    onCheckedChange={(checked) =>
                      setFilterConfig(prev => ({ ...prev, exclude_benign: checked as boolean }))
                    }
                  />
                  <Label htmlFor="exclude_benign">Exclude benign variants</Label>
                </div>
                <div className="flex items-center space-x-2">
                  <Checkbox
                    id="exclude_likely_benign"
                    checked={filterConfig.exclude_likely_benign}
                    onCheckedChange={(checked) =>
                      setFilterConfig(prev => ({ ...prev, exclude_likely_benign: checked as boolean }))
                    }
                  />
                  <Label htmlFor="exclude_likely_benign">Exclude likely benign variants</Label>
                </div>
              </div>
            )}
          </div>

          {/* Population Frequency Filters */}
          <div className="space-y-3">
            <div
              className="flex items-center justify-between cursor-pointer"
              onClick={() => toggleSection('population')}
            >
              <h3 className="text-lg font-semibold flex items-center gap-2">
                🧬 Population Frequency Filters
              </h3>
              {expandedSections.population ? <ChevronUp size={20} /> : <ChevronDown size={20} />}
            </div>
            {expandedSections.population && (
              <div className="space-y-4 pl-6">
                <div className="space-y-2">
                  <Label htmlFor="max_gnomad_af">
                    Max gnomAD AF: {filterConfig.max_gnomad_af.toFixed(3)}
                  </Label>
                  <Slider
                    id="max_gnomad_af"
                    min={0}
                    max={0.1}
                    step={0.001}
                    value={[filterConfig.max_gnomad_af]}
                    onValueChange={([value]) =>
                      setFilterConfig(prev => ({ ...prev, max_gnomad_af: value }))
                    }
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="gnomad_af_column">gnomAD Column</Label>
                  <Select
                    value={filterConfig.gnomad_af_column}
                    onValueChange={(value) =>
                      setFilterConfig(prev => ({ ...prev, gnomad_af_column: value }))
                    }
                  >
                    <SelectTrigger id="gnomad_af_column">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="gnomad40_exome_AF">gnomad40_exome_AF</SelectItem>
                      <SelectItem value="gnomad40_genome_AF">gnomad40_genome_AF</SelectItem>
                      <SelectItem value="ExAC_ALL">ExAC_ALL</SelectItem>
                    </SelectContent>
                  </Select>
                </div>
              </div>
            )}
          </div>

          {/* Quality Filters */}
          <div className="space-y-3">
            <div
              className="flex items-center justify-between cursor-pointer"
              onClick={() => toggleSection('quality')}
            >
              <h3 className="text-lg font-semibold flex items-center gap-2">
                ⭐ Quality Filters
              </h3>
              {expandedSections.quality ? <ChevronUp size={20} /> : <ChevronDown size={20} />}
            </div>
            {expandedSections.quality && (
              <div className="space-y-4 pl-6">
                <div className="space-y-2">
                  <Label htmlFor="min_depth">Min Read Depth (DP)</Label>
                  <Input
                    id="min_depth"
                    type="number"
                    min={0}
                    value={filterConfig.min_depth}
                    onChange={(e) =>
                      setFilterConfig(prev => ({ ...prev, min_depth: parseInt(e.target.value) || 0 }))
                    }
                  />
                </div>
                <div className="flex items-center space-x-2">
                  <Checkbox
                    id="require_pass"
                    checked={filterConfig.require_pass}
                    onCheckedChange={(checked) =>
                      setFilterConfig(prev => ({ ...prev, require_pass: checked as boolean }))
                    }
                  />
                  <Label htmlFor="require_pass">Require PASS filter only</Label>
                </div>
              </div>
            )}
          </div>

          {/* Preview Statistics */}
          {previewStats && (
            <Card className="bg-muted/50">
              <CardHeader>
                <CardTitle className="text-base">Filter Preview</CardTitle>
              </CardHeader>
              <CardContent className="space-y-2">
                <div className="grid grid-cols-2 gap-4">
                  <div>
                    <p className="text-sm font-medium">Total Variants</p>
                    <p className="text-2xl font-bold">{previewStats.total_variants.toLocaleString()}</p>
                  </div>
                  <div>
                    <p className="text-sm font-medium">Filtered Variants</p>
                    <p className="text-2xl font-bold text-primary">
                      {previewStats.filtered_variants.toLocaleString()}
                    </p>
                  </div>
                  <div>
                    <p className="text-sm font-medium">Removed Variants</p>
                    <p className="text-2xl font-bold text-destructive">
                      {previewStats.removed_variants.toLocaleString()}
                    </p>
                  </div>
                  <div>
                    <p className="text-sm font-medium">Reduction</p>
                    <p className="text-2xl font-bold">{previewStats.reduction_percentage}%</p>
                  </div>
                </div>
                {previewStats.top_genes.length > 0 && (
                  <div className="mt-4">
                    <p className="text-sm font-medium mb-2">Top Genes</p>
                    <div className="text-xs space-y-1">
                      {previewStats.top_genes.slice(0, 5).map((gene, idx) => (
                        <div key={idx} className="flex justify-between">
                          <span>{gene.gene}</span>
                          <span className="font-medium">{gene.count}</span>
                        </div>
                      ))}
                    </div>
                  </div>
                )}
              </CardContent>
            </Card>
          )}

          {/* Action Buttons */}
          <div className="flex flex-wrap gap-3">
            <Button onClick={handlePreview} disabled={loading.preview}>
              {loading.preview ? (
                <>
                  <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                  Previewing...
                </>
              ) : (
                <>
                  <FileText className="mr-2 h-4 w-4" />
                  Preview Results
                </>
              )}
            </Button>

            <Button onClick={handleDownload} disabled={loading.download} variant="default">
              {loading.download ? (
                <>
                  <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                  Downloading...
                </>
              ) : (
                <>
                  <Download className="mr-2 h-4 w-4" />
                  Download CSV
                </>
              )}
            </Button>

            <Button onClick={handleClassify} disabled={loading.classify} variant="secondary">
              {loading.classify ? (
                <>
                  <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                  Classifying...
                </>
              ) : (
                <>
                  <FileText className="mr-2 h-4 w-4" />
                  Classify with ACMG
                </>
              )}
            </Button>

            <Button onClick={handleVisualize} disabled={loading.visualize} variant="outline">
              {loading.visualize ? (
                <>
                  <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                  Generating...
                </>
              ) : (
                <>
                  <BarChart3 className="mr-2 h-4 w-4" />
                  Generate Charts
                </>
              )}
            </Button>
          </div>
        </CardContent>
      </Card>
    </div>
  );
}
