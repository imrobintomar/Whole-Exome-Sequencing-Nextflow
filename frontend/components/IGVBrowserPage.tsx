'use client';

import { useEffect, useState, useRef } from 'react';
import { jobApi } from '@/lib/api';
import { ArrowLeft, Search, Loader, AlertCircle, Download, ZoomIn, ZoomOut, RotateCcw, Info } from 'lucide-react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from './ui/card';
import { Button } from './ui/button';
import { Input } from './ui/input';
import { Badge } from './ui/badge';

interface IGVBrowserPageProps {
  jobId: string;
  sampleName: string;
  onBack: () => void;
}

export default function IGVBrowserPage({ jobId, sampleName, onBack }: IGVBrowserPageProps) {
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState('');
  const [searchLocus, setSearchLocus] = useState('chr1:1,211,606-1,211,700');
  const [igvBrowser, setIgvBrowser] = useState<any>(null);
  const igvContainer = useRef<HTMLDivElement>(null);

  useEffect(() => {
    loadIGV();
    return () => {
      // Cleanup IGV instance
      if (igvBrowser) {
        igvBrowser.remove();
      }
    };
  }, []);

  const loadIGV = async () => {
    try {
      setLoading(true);
      setError('');

      // Load IGV.js library dynamically
      if (typeof window !== 'undefined' && !(window as any).igv) {
        const script = document.createElement('script');
        script.src = 'https://cdn.jsdelivr.net/npm/igv@3.7.1/dist/igv.min.js';
        script.async = true;
        await new Promise((resolve, reject) => {
          script.onload = resolve;
          script.onerror = reject;
          document.head.appendChild(script);
        });
      }

      const igv = (window as any).igv;

      // Configure IGV options
      const options = {
        genome: 'hg38',
        locus: 'chr1:1,211,606-1,211,700',
        tracks: [
          {
            name: 'Genes',
            type: 'annotation',
            format: 'refgene',
            url: 'https://s3.amazonaws.com/igv.org.genomes/hg38/refGene.txt.gz',
            indexURL: 'https://s3.amazonaws.com/igv.org.genomes/hg38/refGene.txt.gz.tbi',
            order: 1000000,
            visibilityWindow: 300000000,
            displayMode: 'EXPANDED',
          },
        ],
      };

      // Try to add BAM track if available
      try {
        const bamResponse = await jobApi.checkFileAvailability(jobId, 'bam');
        if (bamResponse) {
          options.tracks.push({
            name: `${sampleName} Alignments`,
            type: 'alignment',
            format: 'bam',
            url: `${process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000'}/jobs/${jobId}/download/bam`,
            indexURL: `${process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000'}/jobs/${jobId}/download/bam.bai`,
            order: 1,
          } as any);
        }
      } catch (err) {
        console.warn('BAM file not available for IGV');
      }

      // Try to add VCF track if available
      try {
        const vcfResponse = await jobApi.checkFileAvailability(jobId, 'vcf');
        if (vcfResponse) {
          options.tracks.push({
            name: `${sampleName} Variants`,
            type: 'variant',
            format: 'vcf',
            url: `${process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000'}/jobs/${jobId}/download/raw_vcf`,
            indexURL: `${process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000'}/jobs/${jobId}/download/raw_vcf.tbi`,
            order: 2,
          } as any);
        }
      } catch (err) {
        console.warn('VCF file not available for IGV');
      }

      // Create IGV browser
      if (igvContainer.current) {
        const browser = await igv.createBrowser(igvContainer.current, options);
        setIgvBrowser(browser);
      }

      setLoading(false);
    } catch (err: any) {
      console.error('IGV initialization error:', err);
      setError(err.message || 'Failed to initialize genome browser');
      setLoading(false);
    }
  };

  const handleSearch = async () => {
    if (igvBrowser && searchLocus) {
      try {
        await igvBrowser.search(searchLocus);
      } catch (err) {
        alert('Invalid locus format. Use format: chr1:1,000,000-2,000,000 or gene name');
      }
    }
  };

  const handleZoomIn = () => {
    if (igvBrowser) {
      igvBrowser.zoomIn();
    }
  };

  const handleZoomOut = () => {
    if (igvBrowser) {
      igvBrowser.zoomOut();
    }
  };

  const handleReset = () => {
    if (igvBrowser) {
      igvBrowser.search('chr1:1,211,606-1,211,700');
      setSearchLocus('chr1:1,211,606-1,211,700');
    }
  };

  const handleKeyPress = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter') {
      handleSearch();
    }
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
                <p className="text-lg font-medium">Loading Genome Browser</p>
                <p className="text-sm text-muted-foreground mt-1">Initializing IGV.js for {sampleName}...</p>
              </div>
            </div>
          </CardContent>
        </Card>
      </div>
    );
  }

  if (error) {
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
              Browser Error
            </CardTitle>
          </CardHeader>
          <CardContent>
            <p className="text-muted-foreground mb-4">{error}</p>
            <Button onClick={loadIGV} variant="outline">
              Retry Loading
            </Button>
          </CardContent>
        </Card>
      </div>
    );
  }

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
            <h1 className="text-3xl font-bold tracking-tight">Integrative Genomics Viewer</h1>
            <p className="text-muted-foreground mt-1">{sampleName}</p>
          </div>
          <Badge variant="outline" className="flex items-center gap-1.5">
            <Info className="h-3 w-3" />
            Reference: hg38
          </Badge>
        </div>
      </div>

      {/* Info Card */}
      <Card className="bg-blue-500/5 border-blue-500/20">
        <CardHeader>
          <CardTitle className="flex items-center gap-2 text-blue-500 text-base">
            <Info className="h-4 w-4" />
            How to Use IGV Browser
          </CardTitle>
        </CardHeader>
        <CardContent className="text-sm text-muted-foreground space-y-2">
          <ul className="list-disc list-inside space-y-1">
            <li>Search by gene name (e.g., BRCA1) or genomic coordinates (e.g., chr17:41,196,312-41,277,500)</li>
            <li>Use the zoom controls to adjust the view resolution</li>
            <li>Click and drag to pan across the genome</li>
            <li>Right-click on tracks for additional options</li>
            <li>Hover over features for detailed information</li>
          </ul>
        </CardContent>
      </Card>

      {/* Controls */}
      <Card>
        <CardHeader>
          <CardTitle className="text-base">Navigation Controls</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="flex gap-3">
            <div className="relative flex-1">
              <Search className="absolute left-3 top-1/2 -translate-y-1/2 h-4 w-4 text-muted-foreground" />
              <Input
                placeholder="Enter gene name or locus (e.g., BRCA1 or chr1:1,000,000-2,000,000)"
                value={searchLocus}
                onChange={e => setSearchLocus(e.target.value)}
                onKeyPress={handleKeyPress}
                className="pl-9"
              />
            </div>
            <Button onClick={handleSearch} variant="default">
              <Search className="h-4 w-4 mr-2" />
              Go
            </Button>
            <Button onClick={handleZoomIn} variant="outline" size="icon" title="Zoom In">
              <ZoomIn className="h-4 w-4" />
            </Button>
            <Button onClick={handleZoomOut} variant="outline" size="icon" title="Zoom Out">
              <ZoomOut className="h-4 w-4" />
            </Button>
            <Button onClick={handleReset} variant="outline" size="icon" title="Reset View">
              <RotateCcw className="h-4 w-4" />
            </Button>
          </div>
        </CardContent>
      </Card>

      {/* IGV Container */}
      <Card>
        <CardHeader>
          <CardTitle className="text-base">Genome Browser</CardTitle>
          <CardDescription>Interactive visualization of alignments and variants</CardDescription>
        </CardHeader>
        <CardContent>
          <div
            ref={igvContainer}
            className="w-full border rounded-lg overflow-hidden bg-white dark:bg-gray-900"
            style={{ minHeight: '600px' }}
          />
        </CardContent>
      </Card>

      {/* Quick Links */}
      <Card>
        <CardHeader>
          <CardTitle className="text-base">Quick Navigation</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-2">
            {[
              { label: 'BRCA1', locus: 'chr17:41,196,312-41,277,500' },
              { label: 'BRCA2', locus: 'chr13:32,889,611-32,973,805' },
              { label: 'TP53', locus: 'chr17:7,661,779-7,687,550' },
              { label: 'EGFR', locus: 'chr7:55,019,032-55,207,338' },
              { label: 'KRAS', locus: 'chr12:25,205,246-25,250,936' },
              { label: 'PTEN', locus: 'chr10:87,863,113-87,971,930' },
              { label: 'APC', locus: 'chr5:112,707,498-112,846,239' },
              { label: 'CDKN2A', locus: 'chr9:21,967,752-21,995,301' },
            ].map(gene => (
              <Button
                key={gene.label}
                variant="outline"
                size="sm"
                onClick={() => {
                  setSearchLocus(gene.locus);
                  if (igvBrowser) {
                    igvBrowser.search(gene.locus);
                  }
                }}
                className="text-xs"
              >
                {gene.label}
              </Button>
            ))}
          </div>
        </CardContent>
      </Card>

      {/* Tracks Info */}
      <Card>
        <CardHeader>
          <CardTitle className="text-base">Available Tracks</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="space-y-2 text-sm">
            <div className="flex items-center gap-2">
              <Badge variant="secondary">RefSeq Genes</Badge>
              <span className="text-muted-foreground">Human gene annotations (hg38)</span>
            </div>
            <div className="flex items-center gap-2">
              <Badge variant="secondary">Alignments</Badge>
              <span className="text-muted-foreground">Aligned reads from BAM file</span>
            </div>
            <div className="flex items-center gap-2">
              <Badge variant="secondary">Variants</Badge>
              <span className="text-muted-foreground">Called variants from VCF file</span>
            </div>
          </div>
        </CardContent>
      </Card>
    </div>
  );
}
