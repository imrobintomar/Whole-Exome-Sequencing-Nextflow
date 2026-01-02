'use client';

import { useEffect, useState } from 'react';
import { jobApi, Job } from '@/lib/api';
import { format } from 'date-fns';
import { Download, RefreshCw, Clock, CheckCircle, XCircle, Loader, Dna } from 'lucide-react';
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

export default function JobList() {
  const [jobs, setJobs] = useState<Job[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState('');

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

  if (loading) {
    return (
      <div className="flex items-center justify-center py-12">
        <Loader className="h-8 w-8 animate-spin text-primary" />
      </div>
    );
  }

  return (
    <div className="space-y-6">
      <div className="flex items-center justify-between">
        <div>
          <h2 className="text-3xl font-bold tracking-tight">My Jobs</h2>
          <p className="text-muted-foreground">
            Track and manage your sequencing pipeline jobs
          </p>
        </div>
        <Button onClick={fetchJobs} variant="outline">
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
            View all your submitted pipeline jobs and download results
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
                  <TableRow key={job.job_id}>
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
                    <TableCell className="text-right">
                      {job.status === 'completed' ? (
                        <div className="flex justify-end gap-2">
                          <Button
                            variant="ghost"
                            onClick={() => jobApi.downloadFile(job.job_id, 'bam')}
                            className="h-8 px-2"
                          >
                            <Download className="mr-1 h-3 w-3" />
                            BAM
                          </Button>
                          <Button
                            variant="ghost"
                            onClick={() => jobApi.downloadFile(job.job_id, 'raw_vcf')}
                            className="h-8 px-2"
                          >
                            <Download className="mr-1 h-3 w-3" />
                            VCF
                          </Button>
                          <Button
                            variant="ghost"
                            onClick={() => jobApi.downloadFile(job.job_id, 'annotated_vcf')}
                            className="h-8 px-2"
                          >
                            <Download className="mr-1 h-3 w-3" />
                            Annotated
                          </Button>
                          <Button
                            variant="ghost"
                            onClick={() => jobApi.downloadFile(job.job_id, 'filtered_tsv')}
                            className="h-8 px-2"
                          >
                            <Download className="mr-1 h-3 w-3" />
                            TSV
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
          )}
        </CardContent>
      </Card>
    </div>
  );
}
