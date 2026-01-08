'use client';

import { useEffect, useState } from 'react';
import { jobApi, Job } from '@/lib/api';
import { format } from 'date-fns';
import { Download, RefreshCw, Clock, CheckCircle, XCircle, Loader, Dna, FlaskConical, Microscope, XOctagon, RotateCcw, PlayCircle, BarChart3, X } from 'lucide-react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from './ui/card';
import { Button } from './ui/button';
import { Badge } from './ui/badge';
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from './ui/table';
import PipelineProgress from './PipelineProgress';

interface JobListProps {
  onJobClick?: (jobId: string) => void;
  onClassifyClick?: (jobId: string, sampleName: string) => void;
  onIGVClick?: (jobId: string, sampleName: string) => void;
  onVariantsClick?: (jobId: string) => void;
}

export default function JobList({ onJobClick, onClassifyClick, onIGVClick, onVariantsClick }: JobListProps = {}) {
  const [jobs, setJobs] = useState<Job[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState('');
  const [downloadingFiles, setDownloadingFiles] = useState<Set<string>>(new Set());
  const [processingActions, setProcessingActions] = useState<Set<string>>(new Set());
  const [selectedJob, setSelectedJob] = useState<Job | null>(null);
  const [sideMenuOpen, setSideMenuOpen] = useState(false);

  useEffect(() => {
    fetchJobs();
    const interval = setInterval(fetchJobs, 10000);
    return () => clearInterval(interval);
  }, []);

  const fetchJobs = async () => {
    try {
      const data = await jobApi.getJobs();
      setJobs(data);
      setError('');
    } catch (err: any) {
      setError('Failed to fetch jobs');
    } finally {
      setLoading(false);
    }
  };

  const getStatusIcon = (status: Job['status']) => {
    switch (status) {
      case 'pending':
        return <Clock className="h-5 w-5 text-yellow-500" />;
      case 'running':
        return <Loader className="h-5 w-5 text-blue-500 animate-spin" />;
      case 'completed':
        return <CheckCircle className="h-5 w-5 text-green-500" />;
      case 'failed':
        return <XCircle className="h-5 w-5 text-red-500" />;
    }
  };

  const getStatusBadgeVariant = (status: Job['status']) => {
    const variants: Record<Job['status'], 'pending' | 'success' | 'destructive' | 'secondary'> = {
      pending: 'pending',
      running: 'secondary',
      completed: 'success',
      failed: 'destructive',
    };
    return variants[status];
  };

  const handleDownload = async (jobId: string, fileType: 'bam' | 'raw_vcf' | 'annotated_vcf' | 'filtered_tsv') => {
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

  const handleCancel = async (jobId: string) => {
    if (processingActions.has(jobId)) return;
    if (!confirm('Are you sure you want to cancel this job?')) return;

    try {
      setProcessingActions(prev => new Set(prev).add(jobId));
      await jobApi.cancelJob(jobId);
      await fetchJobs(); // Refresh jobs list
    } catch (err: any) {
      console.error('Cancel failed:', err);
      alert(`Cancel failed: ${err.response?.data?.detail || err.message || 'Unknown error'}`);
    } finally {
      setProcessingActions(prev => {
        const next = new Set(prev);
        next.delete(jobId);
        return next;
      });
    }
  };

  const handleRerun = async (jobId: string) => {
    if (processingActions.has(jobId)) return;
    if (!confirm('This will delete previous results and restart the pipeline from scratch. Continue?')) return;

    try {
      setProcessingActions(prev => new Set(prev).add(jobId));
      await jobApi.rerunJob(jobId);
      await fetchJobs(); // Refresh jobs list
    } catch (err: any) {
      console.error('Rerun failed:', err);
      alert(`Rerun failed: ${err.response?.data?.detail || err.message || 'Unknown error'}`);
    } finally {
      setProcessingActions(prev => {
        const next = new Set(prev);
        next.delete(jobId);
        return next;
      });
    }
  };

  const handleResume = async (jobId: string) => {
    if (processingActions.has(jobId)) return;

    try {
      setProcessingActions(prev => new Set(prev).add(jobId));
      await jobApi.resumeJob(jobId);
      await fetchJobs(); // Refresh jobs list
    } catch (err: any) {
      console.error('Resume failed:', err);
      alert(`Resume failed: ${err.response?.data?.detail || err.message || 'Unknown error'}`);
    } finally {
      setProcessingActions(prev => {
        const next = new Set(prev);
        next.delete(jobId);
        return next;
      });
    }
  };

  const openSideMenu = (job: Job) => {
    setSelectedJob(job);
    setSideMenuOpen(true);
  };

  const closeSideMenu = () => {
    setSideMenuOpen(false);
    setTimeout(() => setSelectedJob(null), 300); // Delay clearing to allow animation
  };

  if (loading) {
    return (
      <div className="flex items-center justify-center py-12">
        <Loader className="h-8 w-8 animate-spin text-primary" />
      </div>
    );
  }

  return (
    <div className="space-y-6">
      <div className="flex flex-col sm:flex-row sm:items-center sm:justify-between gap-4">
        <div>
          <h2 className="text-2xl sm:text-3xl font-bold tracking-tight">My Jobs</h2>
          <p className="text-sm text-muted-foreground">
            Track and manage your Analysis  jobs
          </p>
        </div>
        <Button onClick={fetchJobs} variant="outline" className="w-full sm:w-auto">
          <RefreshCw className="mr-2 h-4 w-4" />
          Refresh
        </Button>
      </div>

      {error && (
        <Card className="border-destructive/50 bg-destructive/5">
          <CardContent className="pt-6">
            <p className="text-sm text-destructive">{error}</p>
          </CardContent>
        </Card>
      )}

      <Card>
        <CardHeader>
          <CardTitle>All Jobs</CardTitle>
          <CardDescription>
            View all your submitted analysis jobs and download results
          </CardDescription>
        </CardHeader>
        <CardContent>
          {jobs.length === 0 ? (
            <div className="flex flex-col items-center justify-center py-12 text-center">
              <Dna className="mb-4 h-12 w-12 text-muted-foreground" />
              <p className="text-muted-foreground">
                No jobs submitted yet. Upload files to get started!
              </p>
            </div>
          ) : (
            <>
              {/* Mobile Card View */}
              <div className="block md:hidden space-y-4">
                {jobs.map((job) => (
                  <Card
                    key={job.job_id}
                    className="overflow-hidden cursor-pointer hover:border-primary transition-colors"
                    onClick={() => onJobClick?.(job.job_id)}
                  >
                    <CardHeader className="pb-3" onClick={(e) => e.stopPropagation()}>
                      <div className="flex items-start justify-between gap-2">
                        <div className="flex-1 min-w-0">
                          <CardTitle className="text-base truncate">{job.sample_name}</CardTitle>
                          <CardDescription className="text-xs mt-1">
                            {format(new Date(job.created_at), 'PPp')}
                          </CardDescription>
                        </div>
                        <div className="flex items-center gap-2 flex-shrink-0">
                          {getStatusIcon(job.status)}
                          <Badge variant={getStatusBadgeVariant(job.status)} className="text-xs">
                            {job.status}
                          </Badge>
                        </div>
                      </div>
                    </CardHeader>
                    <CardContent className="space-y-3">
                      <div>
                        <PipelineProgress job={job} />
                      </div>

                      {/* Actions for Running Jobs */}
                      {job.status === 'running' && (
                        <div className="flex flex-wrap gap-2">
                          <Button
                            variant="destructive"
                            onClick={() => handleCancel(job.job_id)}
                            disabled={processingActions.has(job.job_id)}
                            className="flex-1 min-w-[120px]"
                            size="sm"
                          >
                            {processingActions.has(job.job_id) ? (
                              <Loader className="mr-2 h-4 w-4 animate-spin" />
                            ) : (
                              <XOctagon className="mr-2 h-4 w-4" />
                            )}
                            Cancel
                          </Button>
                        </div>
                      )}

                      {/* Actions for Failed Jobs */}
                      {job.status === 'failed' && (
                        <div className="flex flex-wrap gap-2">
                          <Button
                            variant="outline"
                            onClick={() => handleResume(job.job_id)}
                            disabled={processingActions.has(job.job_id)}
                            className="flex-1 min-w-[120px]"
                            size="sm"
                          >
                            {processingActions.has(job.job_id) ? (
                              <Loader className="mr-2 h-4 w-4 animate-spin" />
                            ) : (
                              <PlayCircle className="mr-2 h-4 w-4" />
                            )}
                            Resume
                          </Button>
                          <Button
                            variant="outline"
                            onClick={() => handleRerun(job.job_id)}
                            disabled={processingActions.has(job.job_id)}
                            className="flex-1 min-w-[120px]"
                            size="sm"
                          >
                            {processingActions.has(job.job_id) ? (
                              <Loader className="mr-2 h-4 w-4 animate-spin" />
                            ) : (
                              <RotateCcw className="mr-2 h-4 w-4" />
                            )}
                            Rerun
                          </Button>
                        </div>
                      )}

                      {/* Actions for Completed Jobs */}
                      {job.status === 'completed' && (
                        <div className="space-y-2">
                          <Button
                            variant="outline"
                            onClick={() => openSideMenu(job)}
                            className="w-full"
                            size="sm"
                          >
                            <BarChart3 className="mr-2 h-4 w-4" />
                            Analysis Tools
                          </Button>
                          <div className="grid grid-cols-2 gap-2">
                            <Button
                              variant="outline"
                              onClick={() => handleDownload(job.job_id, 'bam')}
                              disabled={downloadingFiles.has(`${job.job_id}-bam`)}
                              size="sm"
                            >
                              {downloadingFiles.has(`${job.job_id}-bam`) ? (
                                <Loader className="mr-2 h-3 w-3 animate-spin" />
                              ) : (
                                <Download className="mr-2 h-3 w-3" />
                              )}
                              BAM
                            </Button>
                            <Button
                              variant="outline"
                              onClick={() => handleDownload(job.job_id, 'raw_vcf')}
                              disabled={downloadingFiles.has(`${job.job_id}-raw_vcf`)}
                              size="sm"
                            >
                              {downloadingFiles.has(`${job.job_id}-raw_vcf`) ? (
                                <Loader className="mr-2 h-3 w-3 animate-spin" />
                              ) : (
                                <Download className="mr-2 h-3 w-3" />
                              )}
                              VCF
                            </Button>
                            <Button
                              variant="outline"
                              onClick={() => handleDownload(job.job_id, 'annotated_vcf')}
                              disabled={downloadingFiles.has(`${job.job_id}-annotated_vcf`)}
                              size="sm"
                            >
                              {downloadingFiles.has(`${job.job_id}-annotated_vcf`) ? (
                                <Loader className="mr-2 h-3 w-3 animate-spin" />
                              ) : (
                                <Download className="mr-2 h-3 w-3" />
                              )}
                              Annotated
                            </Button>
                            <Button
                              variant="outline"
                              onClick={() => handleDownload(job.job_id, 'filtered_tsv')}
                              disabled={downloadingFiles.has(`${job.job_id}-filtered_tsv`)}
                              size="sm"
                            >
                              {downloadingFiles.has(`${job.job_id}-filtered_tsv`) ? (
                                <Loader className="mr-2 h-3 w-3 animate-spin" />
                              ) : (
                                <Download className="mr-2 h-3 w-3" />
                              )}
                              TSV
                            </Button>
                          </div>
                          <Button
                            variant="ghost"
                            onClick={() => handleRerun(job.job_id)}
                            disabled={processingActions.has(job.job_id)}
                            className="w-full"
                            size="sm"
                          >
                            {processingActions.has(job.job_id) ? (
                              <Loader className="mr-2 h-4 w-4 animate-spin" />
                            ) : (
                              <RotateCcw className="mr-2 h-4 w-4" />
                            )}
                            Rerun Pipeline
                          </Button>
                        </div>
                      )}

                      {/* Pending Jobs */}
                      {job.status === 'pending' && (
                        <p className="text-xs text-muted-foreground">Waiting to start...</p>
                      )}

                      {/* Error Message for Failed Jobs */}
                      {job.error_message && (
                        <p className="text-xs text-destructive">{job.error_message}</p>
                      )}
                    </CardContent>
                  </Card>
                ))}
              </div>

              {/* Desktop Table View */}
              <div className="hidden md:block">
                <Table>
              <TableHeader>
                <TableRow>
                  <TableHead>Status</TableHead>
                  <TableHead>Sample Name</TableHead>
                  <TableHead>Pipeline Progress</TableHead>
                  <TableHead>Submitted</TableHead>
                  <TableHead className="text-right">Actions</TableHead>
                </TableRow>
              </TableHeader>
              <TableBody>
                {jobs.map((job) => (
                  <TableRow
                    key={job.job_id}
                    className="cursor-pointer hover:bg-accent/50"
                    onClick={() => onJobClick?.(job.job_id)}
                  >
                    <TableCell>
                      <div className="flex items-center gap-2">
                        {getStatusIcon(job.status)}
                        <Badge variant={getStatusBadgeVariant(job.status)}>
                          {job.status}
                        </Badge>
                      </div>
                    </TableCell>
                    <TableCell className="font-medium">{job.sample_name}</TableCell>
                    <TableCell>
                      <PipelineProgress job={job} />
                    </TableCell>
                    <TableCell className="text-sm">
                      {format(new Date(job.created_at), 'PPp')}
                    </TableCell>
                    <TableCell className="text-right" onClick={(e) => e.stopPropagation()}>
                      {job.status === 'running' ? (
                        <div className="flex justify-end gap-2">
                          <Button
                            variant="ghost"
                            onClick={() => handleCancel(job.job_id)}
                            disabled={processingActions.has(job.job_id)}
                            className="h-8 px-2 text-destructive hover:text-destructive"
                            title="Cancel Pipeline"
                          >
                            {processingActions.has(job.job_id) ? (
                              <Loader className="mr-1 h-3 w-3 animate-spin" />
                            ) : (
                              <XOctagon className="mr-1 h-3 w-3" />
                            )}
                            Cancel
                          </Button>
                        </div>
                      ) : job.status === 'failed' ? (
                        <div className="flex justify-end gap-2">
                          <Button
                            variant="ghost"
                            onClick={() => handleResume(job.job_id)}
                            disabled={processingActions.has(job.job_id)}
                            className="h-8 px-2"
                            title="Resume Pipeline"
                          >
                            {processingActions.has(job.job_id) ? (
                              <Loader className="mr-1 h-3 w-3 animate-spin" />
                            ) : (
                              <PlayCircle className="mr-1 h-3 w-3" />
                            )}
                            Resume
                          </Button>
                          <Button
                            variant="ghost"
                            onClick={() => handleRerun(job.job_id)}
                            disabled={processingActions.has(job.job_id)}
                            className="h-8 px-2"
                            title="Rerun from Scratch"
                          >
                            {processingActions.has(job.job_id) ? (
                              <Loader className="mr-1 h-3 w-3 animate-spin" />
                            ) : (
                              <RotateCcw className="mr-1 h-3 w-3" />
                            )}
                            Rerun
                          </Button>
                        </div>
                      ) : job.status === 'completed' ? (
                        <div className="flex justify-end gap-2">
                          <Button
                            variant="outline"
                            onClick={() => openSideMenu(job)}
                            className="h-8 px-3"
                            title="Analysis Tools"
                          >
                            <BarChart3 className="mr-1 h-4 w-4" />
                            Analysis
                          </Button>
                          <Button
                            variant="ghost"
                            onClick={() => handleDownload(job.job_id, 'bam')}
                            disabled={downloadingFiles.has(`${job.job_id}-bam`)}
                            className="h-8 px-2"
                          >
                            {downloadingFiles.has(`${job.job_id}-bam`) ? (
                              <Loader className="mr-1 h-3 w-3 animate-spin" />
                            ) : (
                              <Download className="mr-1 h-3 w-3" />
                            )}
                            BAM
                          </Button>
                          <Button
                            variant="ghost"
                            onClick={() => handleDownload(job.job_id, 'raw_vcf')}
                            disabled={downloadingFiles.has(`${job.job_id}-raw_vcf`)}
                            className="h-8 px-2"
                          >
                            {downloadingFiles.has(`${job.job_id}-raw_vcf`) ? (
                              <Loader className="mr-1 h-3 w-3 animate-spin" />
                            ) : (
                              <Download className="mr-1 h-3 w-3" />
                            )}
                            VCF
                          </Button>
                          <Button
                            variant="ghost"
                            onClick={() => handleDownload(job.job_id, 'annotated_vcf')}
                            disabled={downloadingFiles.has(`${job.job_id}-annotated_vcf`)}
                            className="h-8 px-2"
                          >
                            {downloadingFiles.has(`${job.job_id}-annotated_vcf`) ? (
                              <Loader className="mr-1 h-3 w-3 animate-spin" />
                            ) : (
                              <Download className="mr-1 h-3 w-3" />
                            )}
                            Annotated
                          </Button>
                          <Button
                            variant="ghost"
                            onClick={() => handleDownload(job.job_id, 'filtered_tsv')}
                            disabled={downloadingFiles.has(`${job.job_id}-filtered_tsv`)}
                            className="h-8 px-2"
                          >
                            {downloadingFiles.has(`${job.job_id}-filtered_tsv`) ? (
                              <Loader className="mr-1 h-3 w-3 animate-spin" />
                            ) : (
                              <Download className="mr-1 h-3 w-3" />
                            )}
                            TSV
                          </Button>
                          <Button
                            variant="ghost"
                            onClick={() => handleRerun(job.job_id)}
                            disabled={processingActions.has(job.job_id)}
                            className="h-8 px-2"
                            title="Rerun Pipeline"
                          >
                            {processingActions.has(job.job_id) ? (
                              <Loader className="mr-1 h-3 w-3 animate-spin" />
                            ) : (
                              <RotateCcw className="mr-1 h-3 w-3" />
                            )}
                            Rerun
                          </Button>
                        </div>
                      ) : job.error_message ? (
                        <p className="text-xs text-destructive">{job.error_message}</p>
                      ) : (
                        <span className="text-xs text-muted-foreground">Processing...</span>
                      )}
                    </TableCell>
                  </TableRow>
                ))}
              </TableBody>
            </Table>
              </div>
            </>
          )}
        </CardContent>
      </Card>

      {/* Side Menu for Analysis Tools */}
      {sideMenuOpen && selectedJob && (
        <>
          {/* Backdrop */}
          <div
            className="fixed inset-0 bg-black/50 z-40 transition-opacity"
            onClick={closeSideMenu}
          />

          {/* Side Panel */}
          <div className="fixed right-0 top-0 h-full w-full sm:w-96 bg-white shadow-2xl z-50 transform transition-transform duration-300 ease-in-out overflow-y-auto">
            {/* Header */}
            <div className="sticky top-0 bg-gradient-to-r from-indigo-600 to-purple-600 text-white p-6 z-10">
              <div className="flex items-start justify-between">
                <div className="flex-1">
                  <h3 className="text-lg font-bold mb-1">Analysis Tools</h3>
                  <p className="text-sm opacity-90">{selectedJob.sample_name}</p>
                </div>
                <button
                  onClick={closeSideMenu}
                  className="p-2 hover:bg-white/20 rounded-lg transition-colors"
                >
                  <X className="h-5 w-5" />
                </button>
              </div>
            </div>

            {/* Content */}
            <div className="p-6 space-y-4">
              {/* Job Info Card */}
              <Card>
                <CardHeader className="pb-3">
                  <CardTitle className="text-sm">Job Information</CardTitle>
                </CardHeader>
                <CardContent className="space-y-2 text-sm">
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">Status:</span>
                    <Badge variant={getStatusBadgeVariant(selectedJob.status)}>
                      {selectedJob.status}
                    </Badge>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">Submitted:</span>
                    <span className="font-medium">
                      {format(new Date(selectedJob.created_at), 'PPp')}
                    </span>
                  </div>
                </CardContent>
              </Card>

              {/* Analysis Actions */}
              <div className="space-y-3">
                <h4 className="text-sm font-semibold text-muted-foreground uppercase tracking-wide">
                  Analysis Options
                </h4>

                <Button
                  variant="outline"
                  className="w-full justify-start h-auto py-4"
                  onClick={() => {
                    onVariantsClick?.(selectedJob.job_id);
                    closeSideMenu();
                  }}
                >
                  <div className="flex items-start gap-3 text-left">
                    <div className="p-2 bg-blue-100 rounded-lg">
                      <BarChart3 className="h-5 w-5 text-blue-600" />
                    </div>
                    <div className="flex-1">
                      <div className="font-semibold">Visualize Variants</div>
                      <div className="text-xs text-muted-foreground mt-1">
                        Interactive variant analysis and visualization
                      </div>
                    </div>
                  </div>
                </Button>

                <Button
                  variant="outline"
                  className="w-full justify-start h-auto py-4"
                  onClick={() => {
                    onIGVClick?.(selectedJob.job_id, selectedJob.sample_name);
                    closeSideMenu();
                  }}
                >
                  <div className="flex items-start gap-3 text-left">
                    <div className="p-2 bg-green-100 rounded-lg">
                      <Microscope className="h-5 w-5 text-green-600" />
                    </div>
                    <div className="flex-1">
                      <div className="font-semibold">IGV Browser</div>
                      <div className="text-xs text-muted-foreground mt-1">
                        View alignments in Integrative Genomics Viewer
                      </div>
                    </div>
                  </div>
                </Button>

                <Button
                  variant="outline"
                  className="w-full justify-start h-auto py-4"
                  onClick={() => {
                    onClassifyClick?.(selectedJob.job_id, selectedJob.sample_name);
                    closeSideMenu();
                  }}
                >
                  <div className="flex items-start gap-3 text-left">
                    <div className="p-2 bg-purple-100 rounded-lg">
                      <FlaskConical className="h-5 w-5 text-purple-600" />
                    </div>
                    <div className="flex-1">
                      <div className="font-semibold">ACMG Classification</div>
                      <div className="text-xs text-muted-foreground mt-1">
                        Classify variants using ACMG guidelines
                      </div>
                    </div>
                  </div>
                </Button>
              </div>

              {/* Download Section */}
              <div className="space-y-3 pt-4 border-t">
                <h4 className="text-sm font-semibold text-muted-foreground uppercase tracking-wide">
                  Download Files
                </h4>

                <div className="grid grid-cols-2 gap-2">
                  <Button
                    variant="outline"
                    onClick={() => handleDownload(selectedJob.job_id, 'bam')}
                    disabled={downloadingFiles.has(`${selectedJob.job_id}-bam`)}
                    size="sm"
                  >
                    {downloadingFiles.has(`${selectedJob.job_id}-bam`) ? (
                      <Loader className="mr-2 h-3 w-3 animate-spin" />
                    ) : (
                      <Download className="mr-2 h-3 w-3" />
                    )}
                    BAM
                  </Button>
                  <Button
                    variant="outline"
                    onClick={() => handleDownload(selectedJob.job_id, 'raw_vcf')}
                    disabled={downloadingFiles.has(`${selectedJob.job_id}-raw_vcf`)}
                    size="sm"
                  >
                    {downloadingFiles.has(`${selectedJob.job_id}-raw_vcf`) ? (
                      <Loader className="mr-2 h-3 w-3 animate-spin" />
                    ) : (
                      <Download className="mr-2 h-3 w-3" />
                    )}
                    VCF
                  </Button>
                  <Button
                    variant="outline"
                    onClick={() => handleDownload(selectedJob.job_id, 'annotated_vcf')}
                    disabled={downloadingFiles.has(`${selectedJob.job_id}-annotated_vcf`)}
                    size="sm"
                  >
                    {downloadingFiles.has(`${selectedJob.job_id}-annotated_vcf`) ? (
                      <Loader className="mr-2 h-3 w-3 animate-spin" />
                    ) : (
                      <Download className="mr-2 h-3 w-3" />
                    )}
                    Annotated
                  </Button>
                  <Button
                    variant="outline"
                    onClick={() => handleDownload(selectedJob.job_id, 'filtered_tsv')}
                    disabled={downloadingFiles.has(`${selectedJob.job_id}-filtered_tsv`)}
                    size="sm"
                  >
                    {downloadingFiles.has(`${selectedJob.job_id}-filtered_tsv`) ? (
                      <Loader className="mr-2 h-3 w-3 animate-spin" />
                    ) : (
                      <Download className="mr-2 h-3 w-3" />
                    )}
                    TSV
                  </Button>
                </div>
              </div>

              {/* Rerun Section */}
              <div className="pt-4 border-t">
                <Button
                  variant="ghost"
                  onClick={() => {
                    handleRerun(selectedJob.job_id);
                    closeSideMenu();
                  }}
                  disabled={processingActions.has(selectedJob.job_id)}
                  className="w-full"
                >
                  {processingActions.has(selectedJob.job_id) ? (
                    <Loader className="mr-2 h-4 w-4 animate-spin" />
                  ) : (
                    <RotateCcw className="mr-2 h-4 w-4" />
                  )}
                  Rerun Pipeline
                </Button>
              </div>
            </div>
          </div>
        </>
      )}
    </div>
  );
}
