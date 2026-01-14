'use client';

import { useEffect, useState } from 'react';
import { jobApi, Job } from '@/lib/api';
import { format } from 'date-fns';
import { formatJobId } from '@/lib/utils';
import {
  Download,
  RefreshCw,
  CheckCircle,
  XCircle,
  Loader,
  Clock,
  ArrowLeft,
  Calendar,
  FileText,
  Database,
  Activity,
  Dna,
  BarChart3,
  Eye,
  AlertCircle,
  PlayCircle,
  RotateCcw,
  XOctagon,
  Filter,
} from 'lucide-react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from './ui/card';
import { Button } from './ui/button';
import { Badge } from './ui/badge';
import { Separator } from './ui/separator';
import { Tabs, TabsContent, TabsList, TabsTrigger } from './ui/tabs';
import PipelineProgress from './PipelineProgress';
import { VariantFilterPanel } from './VariantFilterPanel';
import { ACMGResultsModal } from './ACMGResultsModal';

interface JobDetailsPageProps {
  jobId: string;
  onBack: () => void;
  onClassifyClick?: (jobId: string, sampleName: string) => void;
  onIGVClick?: (jobId: string, sampleName: string) => void;
  onVariantsClick?: (jobId: string) => void;
  onGenePanelClick?: (jobId: string) => void;
}

export default function JobDetailsPage({
  jobId,
  onBack,
  onClassifyClick,
  onIGVClick,
  onVariantsClick,
  onGenePanelClick,
}: JobDetailsPageProps) {
  const [job, setJob] = useState<Job | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState('');
  const [downloadingFiles, setDownloadingFiles] = useState<Set<string>>(new Set());
  const [processing, setProcessing] = useState(false);
  const [acmgResults, setAcmgResults] = useState<any>(null);
  const [showACMGModal, setShowACMGModal] = useState(false);
  const [filteredMetrics, setFilteredMetrics] = useState<any>(null);

  useEffect(() => {
    fetchJob();
    const interval = setInterval(fetchJob, 5000);
    return () => clearInterval(interval);
  }, [jobId]);

  const fetchJob = async () => {
    try {
      const data = await jobApi.getJob(jobId);
      setJob(data);
      setError('');
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to fetch job details');
    } finally {
      setLoading(false);
    }
  };

  const handleDownload = async (fileType: 'bam' | 'raw_vcf' | 'annotated_vcf' | 'filtered_tsv') => {
    const downloadKey = `${jobId}-${fileType}`;
    if (downloadingFiles.has(downloadKey)) return;

    try {
      setDownloadingFiles(prev => new Set(prev).add(downloadKey));
      await jobApi.downloadFile(jobId, fileType);
    } catch (err: any) {
      console.error('Download failed:', err);
      alert(`Download failed: ${err.response?.data?.detail || err.message || 'Unknown error'}`);
    } finally {
      setDownloadingFiles(prev => {
        const next = new Set(prev);
        next.delete(downloadKey);
        return next;
      });
    }
  };

  const handleCancel = async () => {
    if (!confirm('Are you sure you want to cancel this job?')) return;

    try {
      setProcessing(true);
      await jobApi.cancelJob(jobId);
      await fetchJob();
      alert('Job cancelled successfully');
    } catch (err: any) {
      alert(`Failed to cancel job: ${err.response?.data?.detail || err.message}`);
    } finally {
      setProcessing(false);
    }
  };

  const handleRerun = async () => {
    if (!confirm('Are you sure you want to rerun this job?')) return;

    try {
      setProcessing(true);
      await jobApi.rerunJob(jobId);
      await fetchJob();
      alert('Job rerun initiated successfully');
    } catch (err: any) {
      alert(`Failed to rerun job: ${err.response?.data?.detail || err.message}`);
    } finally {
      setProcessing(false);
    }
  };

  const handleResume = async () => {
    if (!confirm('Resume this failed job?')) return;

    try {
      setProcessing(true);
      await jobApi.resumeJob(jobId);
      await fetchJob();
      alert('Job resumed successfully');
    } catch (err: any) {
      alert(`Failed to resume job: ${err.response?.data?.detail || err.message}`);
    } finally {
      setProcessing(false);
    }
  };

  const getStatusColor = (status: Job['status']) => {
    switch (status) {
      case 'pending':
        return 'bg-yellow-500/10 text-yellow-500 border-yellow-500/20';
      case 'running':
        return 'bg-blue-500/10 text-blue-500 border-blue-500/20';
      case 'completed':
        return 'bg-green-500/10 text-green-500 border-green-500/20';
      case 'failed':
        return 'bg-red-500/10 text-red-500 border-red-500/20';
    }
  };

  const getStatusIcon = (status: Job['status']) => {
    switch (status) {
      case 'pending':
        return <Clock className="h-5 w-5" />;
      case 'running':
        return <Loader className="h-5 w-5 animate-spin" />;
      case 'completed':
        return <CheckCircle className="h-5 w-5" />;
      case 'failed':
        return <XCircle className="h-5 w-5" />;
    }
  };

  if (loading) {
    return (
      <div className="flex items-center justify-center h-96">
        <Loader className="h-8 w-8 animate-spin text-primary" />
      </div>
    );
  }

  if (error || !job) {
    return (
      <Card>
        <CardHeader>
          <CardTitle className="flex items-center gap-2 text-red-500">
            <AlertCircle className="h-5 w-5" />
            Error Loading Job
          </CardTitle>
        </CardHeader>
        <CardContent>
          <p className="text-muted-foreground mb-4">{error || 'Job not found'}</p>
          <Button onClick={onBack} variant="outline">
            <ArrowLeft className="h-4 w-4 mr-2" />
            Back to Jobs
          </Button>
        </CardContent>
      </Card>
    );
  }

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-start justify-between">
        <div className="space-y-1">
          <Button onClick={onBack} variant="ghost" size="sm" className="mb-2">
            <ArrowLeft className="h-4 w-4 mr-2" />
            Back to Jobs
          </Button>
          <div className="flex items-center gap-3">
            <h1 className="text-3xl font-bold tracking-tight">{job.sample_name}</h1>
            <Badge className={`${getStatusColor(job.status)} flex items-center gap-1.5 px-3 py-1`}>
              {getStatusIcon(job.status)}
              <span className="capitalize">{job.status}</span>
            </Badge>
          </div>
          <p className="text-sm text-muted-foreground">Job ID: {formatJobId(job.job_id)}</p>
        </div>

        {/* Action Buttons */}
        <div className="flex gap-2">
          {job.status === 'running' && (
            <Button
              onClick={handleCancel}
              variant="destructive"
              disabled={processing}
            >
              <XOctagon className="h-4 w-4 mr-2" />
              Cancel Job
            </Button>
          )}
          {job.status === 'failed' && (
            <>
              <Button
                onClick={handleResume}
                variant="outline"
                disabled={processing}
              >
                <PlayCircle className="h-4 w-4 mr-2" />
                Resume
              </Button>
              <Button
                onClick={handleRerun}
                variant="outline"
                disabled={processing}
              >
                <RotateCcw className="h-4 w-4 mr-2" />
                Rerun
              </Button>
            </>
          )}
          {job.status === 'completed' && (
            <Button
              onClick={handleRerun}
              variant="outline"
              disabled={processing}
            >
              <RotateCcw className="h-4 w-4 mr-2" />
              Rerun
            </Button>
          )}
          <Button onClick={fetchJob} variant="outline" size="icon">
            <RefreshCw className="h-4 w-4" />
          </Button>
        </div>
      </div>

      {/* Overview Cards */}
      <div className="grid gap-4 md:grid-cols-4">
        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Status</CardTitle>
            <Activity className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold capitalize">{job.status}</div>
            <p className="text-xs text-muted-foreground mt-1">{job.current_step || 'Waiting...'}</p>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Created</CardTitle>
            <Calendar className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold">{format(new Date(job.created_at), 'MMM d')}</div>
            <p className="text-xs text-muted-foreground mt-1">
              {format(new Date(job.created_at), 'HH:mm:ss')}
            </p>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Duration</CardTitle>
            <Clock className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold">
              {job.completed_at
                ? Math.round(
                    (new Date(job.completed_at).getTime() - new Date(job.created_at).getTime()) /
                      1000 /
                      60
                  ) + 'm'
                : '-'}
            </div>
            <p className="text-xs text-muted-foreground mt-1">
              {job.status === 'running' ? 'In progress' : job.status === 'completed' ? 'Completed' : 'Not started'}
            </p>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Sample</CardTitle>
            <Dna className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold truncate">{job.sample_name}</div>
            <p className="text-xs text-muted-foreground mt-1">Whole Exome Seq</p>
          </CardContent>
        </Card>
      </div>

      {/* Error Message */}
      {job.status === 'failed' && job.error_message && (
        <Card className="border-red-500/50 bg-red-500/5">
          <CardHeader>
            <CardTitle className="flex items-center gap-2 text-red-500">
              <AlertCircle className="h-5 w-5" />
              Error Details
            </CardTitle>
          </CardHeader>
          <CardContent>
            <p className="text-sm font-mono bg-background p-3 rounded border">{job.error_message}</p>
          </CardContent>
        </Card>
      )}

      {/* Main Content Tabs */}
      <Tabs defaultValue="pipeline" className="space-y-4">
        <TabsList className="grid w-full grid-cols-5">
          <TabsTrigger value="pipeline">Pipeline Progress</TabsTrigger>
          <TabsTrigger value="files">Files & Downloads</TabsTrigger>
          <TabsTrigger value="filtering" disabled={job?.status !== 'completed'}>
            Variant Filtering
          </TabsTrigger>
          <TabsTrigger value="analysis">Analysis Tools</TabsTrigger>
          <TabsTrigger value="details">Details</TabsTrigger>
        </TabsList>

        {/* Pipeline Tab */}
        <TabsContent value="pipeline" className="space-y-4">
          <Card>
            <CardHeader>
              <CardTitle>Pipeline Execution</CardTitle>
              <CardDescription>Real-time progress of the WES analysis pipeline</CardDescription>
            </CardHeader>
            <CardContent>
              <PipelineProgress job={job} />
            </CardContent>
          </Card>
        </TabsContent>

        {/* Files Tab */}
        <TabsContent value="files" className="space-y-4">
          <Card>
            <CardHeader>
              <CardTitle>Output Files</CardTitle>
              <CardDescription>Download results from completed pipeline stages</CardDescription>
            </CardHeader>
            <CardContent>
              {job.status !== 'completed' ? (
                <div className="text-center py-8 text-muted-foreground">
                  <Database className="h-12 w-12 mx-auto mb-3 opacity-50" />
                  <p>Files will be available when the job completes</p>
                </div>
              ) : (
                <div className="grid gap-3">
                  {[
                    { type: 'bam' as const, label: 'BAM File', desc: 'Aligned and processed reads', icon: Database },
                    { type: 'raw_vcf' as const, label: 'Raw VCF', desc: 'Raw variant calls', icon: FileText },
                    {
                      type: 'annotated_vcf' as const,
                      label: 'Annotated VCF',
                      desc: 'ANNOVAR annotated variants',
                      icon: FileText,
                    },
                    { type: 'filtered_tsv' as const, label: 'Final TSV', desc: 'Processed variant table', icon: FileText },
                  ].map(({ type, label, desc, icon: Icon }) => (
                    <div
                      key={type}
                      className="flex items-center justify-between p-4 border rounded-lg hover:bg-accent transition-colors"
                    >
                      <div className="flex items-center gap-3">
                        <div className="p-2 rounded-lg bg-primary/10">
                          <Icon className="h-5 w-5 text-primary" />
                        </div>
                        <div>
                          <p className="font-medium">{label}</p>
                          <p className="text-sm text-muted-foreground">{desc}</p>
                        </div>
                      </div>
                      <Button
                        onClick={() => handleDownload(type)}
                        disabled={downloadingFiles.has(`${jobId}-${type}`)}
                        variant="outline"
                        size="sm"
                      >
                        {downloadingFiles.has(`${jobId}-${type}`) ? (
                          <Loader className="h-4 w-4 animate-spin" />
                        ) : (
                          <>
                            <Download className="h-4 w-4 mr-2" />
                            Download
                          </>
                        )}
                      </Button>
                    </div>
                  ))}
                </div>
              )}
            </CardContent>
          </Card>
        </TabsContent>

        {/* Variant Filtering Tab */}
        <TabsContent value="filtering" className="space-y-4">
          {job && (
            <VariantFilterPanel
              jobId={jobId}
              onClassifyResults={(results) => {
                setAcmgResults(results);
                setShowACMGModal(true);
              }}
              onVisualizeResults={(metrics) => {
                setFilteredMetrics(metrics);
                // Optionally switch to a visualization view or update charts
              }}
            />
          )}
        </TabsContent>

        {/* Analysis Tab */}
        <TabsContent value="analysis" className="space-y-4">
          <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-4">
            <Card className="cursor-pointer hover:border-primary transition-colors" onClick={() => onVariantsClick?.(jobId)}>
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <BarChart3 className="h-5 w-5 text-primary" />
                  Variant Visualization
                </CardTitle>
                <CardDescription>Explore variant distributions and metrics</CardDescription>
              </CardHeader>
              <CardContent>
                <Button className="w-full" disabled={job.status !== 'completed'}>
                  <Eye className="h-4 w-4 mr-2" />
                  View Visualizations
                </Button>
              </CardContent>
            </Card>

            <Card
              className="cursor-pointer hover:border-primary transition-colors"
              onClick={() => onClassifyClick?.(jobId, job.sample_name)}
            >
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <Dna className="h-5 w-5 text-primary" />
                  ACMG Classification
                </CardTitle>
                <CardDescription>Classify variants using ACMG 2015 guidelines</CardDescription>
              </CardHeader>
              <CardContent>
                <Button className="w-full" disabled={job.status !== 'completed'}>
                  <Eye className="h-4 w-4 mr-2" />
                  Classify Variants
                </Button>
              </CardContent>
            </Card>

            <Card
              className="cursor-pointer hover:border-primary transition-colors"
              onClick={() => onIGVClick?.(jobId, job.sample_name)}
            >
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <Activity className="h-5 w-5 text-primary" />
                  IGV Browser
                </CardTitle>
                <CardDescription>Visualize alignments and variants in genome browser</CardDescription>
              </CardHeader>
              <CardContent>
                <Button className="w-full" disabled={job.status !== 'completed'}>
                  <Eye className="h-4 w-4 mr-2" />
                  Open IGV
                </Button>
              </CardContent>
            </Card>

            <Card
              className="cursor-pointer hover:border-primary transition-colors"
              onClick={() => onGenePanelClick?.(jobId)}
            >
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <Filter className="h-5 w-5 text-primary" />
                  Gene Panel Filter
                </CardTitle>
                <CardDescription>Filter variants by gene panels (ACMG SF, PanelApp)</CardDescription>
              </CardHeader>
              <CardContent>
                <Button className="w-full" disabled={job.status !== 'completed'}>
                  <Eye className="h-4 w-4 mr-2" />
                  Apply Gene Panel
                </Button>
              </CardContent>
            </Card>
          </div>
        </TabsContent>

        {/* Details Tab */}
        <TabsContent value="details" className="space-y-4">
          <Card>
            <CardHeader>
              <CardTitle>Job Information</CardTitle>
            </CardHeader>
            <CardContent className="space-y-4">
              <div className="grid grid-cols-2 gap-4">
                <div>
                  <p className="text-sm font-medium text-muted-foreground">Job ID</p>
                  <p className="text-sm font-mono mt-1">{formatJobId(job.job_id)}</p>
                </div>
                <div>
                  <p className="text-sm font-medium text-muted-foreground">Sample Name</p>
                  <p className="text-sm mt-1">{job.sample_name}</p>
                </div>
                <div>
                  <p className="text-sm font-medium text-muted-foreground">Created At</p>
                  <p className="text-sm mt-1">{format(new Date(job.created_at), 'PPpp')}</p>
                </div>
                {job.started_at && (
                  <div>
                    <p className="text-sm font-medium text-muted-foreground">Started At</p>
                    <p className="text-sm mt-1">{format(new Date(job.started_at), 'PPpp')}</p>
                  </div>
                )}
                {job.completed_at && (
                  <div>
                    <p className="text-sm font-medium text-muted-foreground">Completed At</p>
                    <p className="text-sm mt-1">{format(new Date(job.completed_at), 'PPpp')}</p>
                  </div>
                )}
                <div>
                  <p className="text-sm font-medium text-muted-foreground">Current Step</p>
                  <p className="text-sm mt-1">{job.current_step || 'N/A'}</p>
                </div>
              </div>
            </CardContent>
          </Card>
        </TabsContent>
      </Tabs>

      {/* ACMG Results Modal */}
      <ACMGResultsModal
        open={showACMGModal}
        onClose={() => setShowACMGModal(false)}
        results={acmgResults}
      />
    </div>
  );
}
