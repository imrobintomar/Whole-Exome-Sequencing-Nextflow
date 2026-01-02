'use client';

import { useEffect, useState } from 'react';
import { jobApi, Job } from '@/lib/api';
import { format } from 'date-fns';
import { Download, RefreshCw, Clock, CheckCircle, XCircle, Loader } from 'lucide-react';

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

  const getStatusBadge = (status: Job['status']) => {
    const badges = {
      pending: 'bg-yellow-100 text-yellow-800',
      running: 'bg-blue-100 text-blue-800',
      completed: 'bg-green-100 text-green-800',
      failed: 'bg-red-100 text-red-800',
    };
    return badges[status];
  };

  if (loading) {
    return (
      <div className="bg-white rounded-lg shadow-xl p-8 text-center">
        <Loader className="h-8 w-8 text-blue-500 animate-spin mx-auto mb-4" />
        <p className="text-gray-600">Loading jobs...</p>
      </div>
    );
  }

  return (
    <div className="bg-white rounded-lg shadow-xl p-8">
      <div className="flex justify-between items-center mb-6">
        <h2 className="text-2xl font-bold text-gray-800">My Jobs</h2>
        <button
          onClick={fetchJobs}
          className="flex items-center space-x-2 px-4 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700 transition"
        >
          <RefreshCw className="h-4 w-4" />
          <span>Refresh</span>
        </button>
      </div>

      {error && (
        <div className="bg-red-50 text-red-600 p-3 rounded-lg text-sm mb-4">
          {error}
        </div>
      )}

      {jobs.length === 0 ? (
        <p className="text-center text-gray-500 py-8">
          No jobs submitted yet. Upload files to get started!
        </p>
      ) : (
        <div className="space-y-4">
          {jobs.map((job) => (
            <div
              key={job.job_id}
              className="border border-gray-200 rounded-lg p-4 hover:shadow-md transition"
            >
              <div className="flex items-start justify-between">
                <div className="flex-1">
                  <div className="flex items-center space-x-3 mb-2">
                    {getStatusIcon(job.status)}
                    <h3 className="text-lg font-semibold text-gray-800">
                      {job.sample_name}
                    </h3>
                    <span
                      className={`px-2 py-1 text-xs font-medium rounded-full ${getStatusBadge(
                        job.status
                      )}`}
                    >
                      {job.status.toUpperCase()}
                    </span>
                  </div>

                  <p className="text-sm text-gray-500 mb-1">
                    Job ID: {job.job_id}
                  </p>
                  <p className="text-sm text-gray-500">
                    Submitted: {format(new Date(job.created_at), 'PPpp')}
                  </p>
                  {job.completed_at && (
                    <p className="text-sm text-gray-500">
                      Completed: {format(new Date(job.completed_at), 'PPpp')}
                    </p>
                  )}
                  {job.error_message && (
                    <p className="text-sm text-red-600 mt-2">
                      Error: {job.error_message}
                    </p>
                  )}
                </div>

                {job.status === 'completed' && (
                  <div className="flex flex-col space-y-2">
                    <button
                      onClick={() => jobApi.downloadFile(job.job_id, 'bam')}
                      className="flex items-center space-x-2 px-3 py-1 bg-green-600 text-white rounded text-sm hover:bg-green-700 transition"
                    >
                      <Download className="h-4 w-4" />
                      <span>BAM</span>
                    </button>
                    <button
                      onClick={() => jobApi.downloadFile(job.job_id, 'raw_vcf')}
                      className="flex items-center space-x-2 px-3 py-1 bg-green-600 text-white rounded text-sm hover:bg-green-700 transition"
                    >
                      <Download className="h-4 w-4" />
                      <span>Raw VCF</span>
                    </button>
                    <button
                      onClick={() => jobApi.downloadFile(job.job_id, 'annotated_vcf')}
                      className="flex items-center space-x-2 px-3 py-1 bg-green-600 text-white rounded text-sm hover:bg-green-700 transition"
                    >
                      <Download className="h-4 w-4" />
                      <span>Annotated VCF</span>
                    </button>
                    <button
                      onClick={() => jobApi.downloadFile(job.job_id, 'filtered_tsv')}
                      className="flex items-center space-x-2 px-3 py-1 bg-green-600 text-white rounded text-sm hover:bg-green-700 transition"
                    >
                      <Download className="h-4 w-4" />
                      <span>Filtered TSV</span>
                    </button>
                  </div>
                )}
              </div>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
