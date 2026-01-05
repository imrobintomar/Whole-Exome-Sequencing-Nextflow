'use client';

import { useEffect, useRef, useState } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from './ui/card';
import { Button } from './ui/button';
import { Input } from './ui/input';
import { Label } from './ui/label';
import { Loader, AlertCircle } from 'lucide-react';
import { auth } from '@/lib/firebase';
import Script from 'next/script';

interface IGVBrowserProps {
  jobId: string;
  sampleName: string;
  className?: string;
}

export default function IGVBrowser({ jobId, sampleName, className }: IGVBrowserProps) {
  const igvContainerRef = useRef<HTMLDivElement>(null);
  const browserRef = useRef<any>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [locusInput, setLocusInput] = useState('chr17:43,044,295-43,125,364'); // BRCA1 default
  const [igvLoaded, setIgvLoaded] = useState(false);

  useEffect(() => {
    if (igvLoaded) {
      loadIGV();
    }

    return () => {
      // Cleanup IGV instance on unmount
      if (browserRef.current) {
        try {
          browserRef.current.remove();
        } catch (e) {
          console.error('Error cleaning up IGV:', e);
        }
      }
    };
  }, [jobId, igvLoaded]);

  const loadIGV = async () => {
    try {
      setLoading(true);
      setError(null);

      // IGV is loaded via script tag, accessible as window.igv
      const igv = (window as any).igv;

      if (!igvContainerRef.current) {
        throw new Error('IGV container not found');
      }

      if (!igv || typeof igv.createBrowser !== 'function') {
        throw new Error('IGV.js not loaded. Please refresh the page.');
      }

      // Get authentication token from Firebase
      const user = auth.currentUser;
      if (!user) {
        throw new Error('Not authenticated');
      }

      const token = await user.getIdToken();
      const API_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

      // IGV configuration
      const igvOptions: any = {
        genome: 'hg38',
        locus: 'chr17:43,044,295-43,125,364', // BRCA1
        loadDefaultGenome: true,
        tracks: [
          // BAM track
          {
            name: `${sampleName} - Aligned Reads`,
            type: 'alignment',
            format: 'bam',
            url: `${API_URL}/jobs/${jobId}/download/bam`,
            indexURL: `${API_URL}/jobs/${jobId}/download/bam?index=true`,
            height: 300,
            autoHeight: false,
            visibilityWindow: 1000000,
            oauthToken: token
          },
          // Raw VCF track (using raw VCF as it has a tabix index)
          {
            name: `${sampleName} - Variants`,
            type: 'variant',
            format: 'vcf',
            url: `${API_URL}/jobs/${jobId}/download/raw_vcf`,
            indexURL: `${API_URL}/jobs/${jobId}/download/raw_vcf?index=true`,
            displayMode: 'EXPANDED',
            oauthToken: token
          }
        ]
      };

      // Create IGV browser
      browserRef.current = await igv.createBrowser(igvContainerRef.current, igvOptions);

      setLoading(false);
    } catch (err: any) {
      console.error('IGV initialization error:', err);
      setError(err.message || 'Failed to load IGV browser');
      setLoading(false);
    }
  };

  const handleLocusChange = () => {
    if (browserRef.current && locusInput) {
      try {
        browserRef.current.search(locusInput);
      } catch (err: any) {
        setError(`Invalid locus: ${err.message}`);
      }
    }
  };

  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter') {
      handleLocusChange();
    }
  };

  return (
    <div className={className}>
      {/* Load IGV.js from CDN */}
      <Script
        src="https://cdn.jsdelivr.net/npm/igv@3.7.1/dist/igv.min.js"
        onLoad={() => setIgvLoaded(true)}
        onError={() => setError('Failed to load IGV.js library')}
      />
      <Card>
        <CardHeader>
          <CardTitle>Genome Browser - {sampleName}</CardTitle>
          <CardDescription>
            Interactive visualization of aligned reads (BAM) and variants (VCF)
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          {/* Locus navigation */}
          <div className="flex gap-2">
            <div className="flex-1">
              <Label htmlFor="locus">Navigate to locus</Label>
              <div className="flex gap-2 mt-1">
                <Input
                  id="locus"
                  type="text"
                  value={locusInput}
                  onChange={(e) => setLocusInput(e.target.value)}
                  onKeyDown={handleKeyDown}
                  placeholder="chr:start-end or gene name"
                  className="font-mono text-sm"
                />
                <Button onClick={handleLocusChange} disabled={loading}>
                  Go
                </Button>
              </div>
            </div>
          </div>

          {/* Quick navigation buttons */}
          <div className="flex gap-2 flex-wrap">
            <Button
              variant="outline"
              size="sm"
              onClick={() => {
                setLocusInput('chr17:43,044,295-43,125,364');
                browserRef.current?.search('chr17:43,044,295-43,125,364');
              }}
              disabled={loading}
            >
              BRCA1
            </Button>
            <Button
              variant="outline"
              size="sm"
              onClick={() => {
                setLocusInput('chr13:32,315,474-32,400,266');
                browserRef.current?.search('chr13:32,315,474-32,400,266');
              }}
              disabled={loading}
            >
              BRCA2
            </Button>
            <Button
              variant="outline"
              size="sm"
              onClick={() => {
                setLocusInput('chr17:7,668,402-7,687,490');
                browserRef.current?.search('chr17:7,668,402-7,687,490');
              }}
              disabled={loading}
            >
              TP53
            </Button>
            <Button
              variant="outline"
              size="sm"
              onClick={() => {
                setLocusInput('chr7:55,019,017-55,211,628');
                browserRef.current?.search('chr7:55,019,017-55,211,628');
              }}
              disabled={loading}
            >
              EGFR
            </Button>
          </div>

          {/* Error message */}
          {error && (
            <div className="flex items-center gap-2 p-3 bg-destructive/10 border border-destructive/30 rounded-lg">
              <AlertCircle className="h-4 w-4 text-destructive flex-shrink-0" />
              <p className="text-sm text-destructive">{error}</p>
              <Button
                variant="outline"
                size="sm"
                onClick={loadIGV}
                className="ml-auto"
              >
                Retry
              </Button>
            </div>
          )}

          {/* IGV Container */}
          <div className="relative">
            {loading && (
              <div className="flex items-center justify-center py-12 bg-muted rounded-lg">
                <div className="text-center space-y-2">
                  <Loader className="h-8 w-8 animate-spin text-primary mx-auto" />
                  <p className="text-sm text-muted-foreground">Loading genome browser...</p>
                </div>
              </div>
            )}
            <div
              ref={igvContainerRef}
              className={loading ? 'hidden' : 'border rounded-lg overflow-hidden'}
              style={{ minHeight: '500px' }}
            />
          </div>

          {/* Instructions */}
          <div className="text-xs text-muted-foreground space-y-1 p-3 bg-muted rounded">
            <p><strong>Navigation:</strong></p>
            <ul className="list-disc list-inside space-y-1 ml-2">
              <li>Zoom: Mouse wheel or +/- buttons</li>
              <li>Pan: Click and drag in the track area</li>
              <li>Click on reads/variants for detailed information</li>
              <li>Right-click for additional options</li>
            </ul>
          </div>
        </CardContent>
      </Card>
    </div>
  );
}
