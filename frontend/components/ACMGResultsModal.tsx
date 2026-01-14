import React, { useState, useMemo } from 'react';
import {
  Dialog,
  DialogContent,
  DialogDescription,
  DialogHeader,
  DialogTitle,
} from './ui/dialog';
import { Card, CardContent, CardHeader, CardTitle } from './ui/card';
import { Button } from './ui/button';
import { Input } from './ui/input';
import { Label } from './ui/label';
import { Badge } from './ui/badge';
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from './ui/select';
import { Download, Search } from 'lucide-react';
import { PieChart, Pie, Cell, ResponsiveContainer, Legend, Tooltip } from 'recharts';

interface ACMGVariant {
  gene: string;
  variant: string;
  classification: string;
  evidence: string[];
  score: number;
}

interface ACMGResults {
  job_id: string;
  sample_name: string;
  total_variants: number;
  classifications: {
    pathogenic: number;
    likely_pathogenic: number;
    vus: number;
    likely_benign: number;
    benign: number;
  };
  variants: ACMGVariant[];
}

interface ACMGResultsModalProps {
  open: boolean;
  onClose: () => void;
  results: ACMGResults | null;
}

const CLASSIFICATION_COLORS: Record<string, string> = {
  pathogenic: '#ef4444',          // red-500
  likely_pathogenic: '#f97316',   // orange-500
  vus: '#eab308',                  // yellow-500
  likely_benign: '#84cc16',        // lime-500
  benign: '#22c55e',               // green-500
};

const CLASSIFICATION_LABELS: Record<string, string> = {
  pathogenic: 'Pathogenic (P)',
  likely_pathogenic: 'Likely Pathogenic (LP)',
  vus: 'VUS',
  likely_benign: 'Likely Benign (LB)',
  benign: 'Benign (B)',
};

export function ACMGResultsModal({ open, onClose, results }: ACMGResultsModalProps) {
  const [searchTerm, setSearchTerm] = useState('');
  const [filterClassification, setFilterClassification] = useState<string>('all');

  const chartData = useMemo(() => {
    if (!results) return [];

    return Object.entries(results.classifications)
      .filter(([_, count]) => count > 0)
      .map(([classification, count]) => ({
        name: CLASSIFICATION_LABELS[classification] || classification,
        value: count,
        color: CLASSIFICATION_COLORS[classification],
      }));
  }, [results]);

  const filteredVariants = useMemo(() => {
    if (!results) return [];

    let filtered = results.variants;

    // Filter by classification
    if (filterClassification !== 'all') {
      filtered = filtered.filter(v =>
        v.classification.toLowerCase().replace(' ', '_') === filterClassification
      );
    }

    // Filter by search term
    if (searchTerm) {
      const term = searchTerm.toLowerCase();
      filtered = filtered.filter(v =>
        v.gene?.toLowerCase().includes(term) ||
        v.variant?.toLowerCase().includes(term) ||
        v.classification?.toLowerCase().includes(term)
      );
    }

    return filtered;
  }, [results, searchTerm, filterClassification]);

  const handleDownloadCSV = () => {
    if (!results) return;

    const headers = ['Gene', 'Variant', 'Classification', 'Score', 'Evidence'];
    const rows = results.variants.map(v => [
      v.gene || '',
      v.variant || '',
      v.classification || '',
      v.score?.toString() || '',
      v.evidence?.join('; ') || '',
    ]);

    const csv = [
      headers.join(','),
      ...rows.map(row => row.map(cell => `"${cell}"`).join(',')),
    ].join('\n');

    const blob = new Blob([csv], { type: 'text/csv' });
    const url = window.URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = `${results.sample_name}_acmg_classifications.csv`;
    document.body.appendChild(link);
    link.click();
    link.remove();
    window.URL.revokeObjectURL(url);
  };

  const getClassificationBadgeVariant = (classification: string): "default" | "destructive" | "secondary" | "outline" => {
    const normalized = classification.toLowerCase().replace(' ', '_');
    if (normalized === 'pathogenic') return 'destructive';
    if (normalized === 'likely_pathogenic') return 'destructive';
    if (normalized === 'vus') return 'secondary';
    return 'outline';
  };

  if (!results) return null;

  return (
    <Dialog open={open} onOpenChange={onClose}>
      <DialogContent className="max-w-6xl max-h-[90vh] overflow-y-auto">
        <DialogHeader>
          <DialogTitle>ACMG Classification Results</DialogTitle>
          <DialogDescription>
            Sample: {results.sample_name} | Total Variants: {results.total_variants}
          </DialogDescription>
        </DialogHeader>

        <div className="space-y-6">
          {/* Classification Summary Chart */}
          <Card>
            <CardHeader>
              <CardTitle className="text-base">Classification Summary</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                {/* Pie Chart */}
                <div className="h-64">
                  <ResponsiveContainer width="100%" height="100%">
                    <PieChart>
                      <Pie
                        data={chartData}
                        cx="50%"
                        cy="50%"
                        labelLine={false}
                        label={({ name, value }) => `${name}: ${value}`}
                        outerRadius={80}
                        fill="#8884d8"
                        dataKey="value"
                      >
                        {chartData.map((entry, index) => (
                          <Cell key={`cell-${index}`} fill={entry.color} />
                        ))}
                      </Pie>
                      <Tooltip />
                      <Legend />
                    </PieChart>
                  </ResponsiveContainer>
                </div>

                {/* Classification Counts */}
                <div className="space-y-2">
                  <h4 className="font-semibold mb-3">Classification Breakdown</h4>
                  {Object.entries(results.classifications).map(([classification, count]) => (
                    <div key={classification} className="flex justify-between items-center p-2 rounded border">
                      <span className="text-sm font-medium capitalize">
                        {CLASSIFICATION_LABELS[classification] || classification}
                      </span>
                      <Badge
                        variant={getClassificationBadgeVariant(classification)}
                        className="ml-2"
                      >
                        {count}
                      </Badge>
                    </div>
                  ))}
                </div>
              </div>
            </CardContent>
          </Card>

          {/* Filters and Search */}
          <Card>
            <CardContent className="pt-6">
              <div className="flex flex-col md:flex-row gap-4">
                <div className="flex-1">
                  <Label htmlFor="search">Search Variants</Label>
                  <div className="relative">
                    <Search className="absolute left-2 top-2.5 h-4 w-4 text-muted-foreground" />
                    <Input
                      id="search"
                      placeholder="Search by gene or variant..."
                      value={searchTerm}
                      onChange={(e) => setSearchTerm(e.target.value)}
                      className="pl-8"
                    />
                  </div>
                </div>
                <div className="w-full md:w-48">
                  <Label htmlFor="filter-classification">Filter by Classification</Label>
                  <Select value={filterClassification} onValueChange={setFilterClassification}>
                    <SelectTrigger id="filter-classification">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="all">All Classifications</SelectItem>
                      <SelectItem value="pathogenic">Pathogenic (P)</SelectItem>
                      <SelectItem value="likely_pathogenic">Likely Pathogenic (LP)</SelectItem>
                      <SelectItem value="vus">VUS</SelectItem>
                      <SelectItem value="likely_benign">Likely Benign (LB)</SelectItem>
                      <SelectItem value="benign">Benign (B)</SelectItem>
                    </SelectContent>
                  </Select>
                </div>
                <div className="flex items-end">
                  <Button onClick={handleDownloadCSV} variant="outline">
                    <Download className="mr-2 h-4 w-4" />
                    Download CSV
                  </Button>
                </div>
              </div>
            </CardContent>
          </Card>

          {/* Variants Table */}
          <Card>
            <CardHeader>
              <CardTitle className="text-base">
                Classified Variants ({filteredVariants.length})
              </CardTitle>
            </CardHeader>
            <CardContent>
              <div className="overflow-x-auto">
                <table className="w-full text-sm">
                  <thead>
                    <tr className="border-b">
                      <th className="text-left p-2">Gene</th>
                      <th className="text-left p-2">Variant</th>
                      <th className="text-left p-2">Classification</th>
                      <th className="text-left p-2">Score</th>
                      <th className="text-left p-2">Evidence</th>
                    </tr>
                  </thead>
                  <tbody>
                    {filteredVariants.length === 0 ? (
                      <tr>
                        <td colSpan={5} className="text-center p-4 text-muted-foreground">
                          No variants match the current filters
                        </td>
                      </tr>
                    ) : (
                      filteredVariants.map((variant, idx) => (
                        <tr key={idx} className="border-b hover:bg-muted/50">
                          <td className="p-2 font-medium">{variant.gene || '-'}</td>
                          <td className="p-2 font-mono text-xs">{variant.variant || '-'}</td>
                          <td className="p-2">
                            <Badge variant={getClassificationBadgeVariant(variant.classification)}>
                              {variant.classification}
                            </Badge>
                          </td>
                          <td className="p-2">{variant.score?.toFixed(1) || '-'}</td>
                          <td className="p-2 text-xs max-w-md truncate" title={variant.evidence?.join(', ')}>
                            {variant.evidence?.join(', ') || '-'}
                          </td>
                        </tr>
                      ))
                    )}
                  </tbody>
                </table>
              </div>
            </CardContent>
          </Card>
        </div>
      </DialogContent>
    </Dialog>
  );
}
