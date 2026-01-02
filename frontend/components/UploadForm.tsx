'use client';

import { useState } from 'react';
import { useDropzone } from 'react-dropzone';
import { jobApi } from '@/lib/api';
import { Upload, FileText, Loader } from 'lucide-react';

interface UploadFormProps {
  onJobSubmitted: () => void;
}

export default function UploadForm({ onJobSubmitted }: UploadFormProps) {
  const [sampleName, setSampleName] = useState('');
  const [fastqR1, setFastqR1] = useState<File | null>(null);
  const [fastqR2, setFastqR2] = useState<File | null>(null);
  const [uploading, setUploading] = useState(false);
  const [uploadProgress, setUploadProgress] = useState(0);
  const [error, setError] = useState('');
  const [success, setSuccess] = useState('');

  const onDropR1 = (acceptedFiles: File[]) => {
    if (acceptedFiles.length > 0) {
      setFastqR1(acceptedFiles[0]);
      setError('');
    }
  };

  const onDropR2 = (acceptedFiles: File[]) => {
    if (acceptedFiles.length > 0) {
      setFastqR2(acceptedFiles[0]);
      setError('');
    }
  };

  const dropzoneR1 = useDropzone({
    onDrop: onDropR1,
    accept: { 'application/gzip': ['.gz', '.fastq.gz', '.fq.gz'] },
    maxFiles: 1,
  });

  const dropzoneR2 = useDropzone({
    onDrop: onDropR2,
    accept: { 'application/gzip': ['.gz', '.fastq.gz', '.fq.gz'] },
    maxFiles: 1,
  });

  const formatFileSize = (bytes: number) => {
    if (bytes === 0) return '0 Bytes';
    const k = 1024;
    const sizes = ['Bytes', 'KB', 'MB', 'GB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return Math.round(bytes / Math.pow(k, i) * 100) / 100 + ' ' + sizes[i];
  };

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setError('');
    setSuccess('');
    setUploadProgress(0);

    if (!sampleName || !fastqR1 || !fastqR2) {
      setError('Please provide sample name and both FASTQ files');
      return;
    }

    setUploading(true);

    try {
      await jobApi.submitJob(sampleName, fastqR1, fastqR2, (progress) => {
        setUploadProgress(progress);
      });
      setSuccess('Job submitted successfully! Check the "My Jobs" tab to monitor progress.');
      setSampleName('');
      setFastqR1(null);
      setFastqR2(null);
      setUploadProgress(0);
      setTimeout(() => {
        onJobSubmitted();
      }, 1500);
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to submit job');
      setUploadProgress(0);
    } finally {
      setUploading(false);
    }
  };

  return (
    <div className="bg-white rounded-lg shadow-xl p-8">
      <h2 className="text-2xl font-bold mb-6 text-gray-800">
        Submit New Analysis
      </h2>

      <form onSubmit={handleSubmit} className="space-y-6">
        <div>
          <label className="block text-sm font-medium text-gray-700 mb-2">
            Sample Name
          </label>
          <input
            type="text"
            value={sampleName}
            onChange={(e) => setSampleName(e.target.value)}
            placeholder="e.g., Sample001"
            className="w-full px-4 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-blue-500 focus:border-transparent"
            required
            disabled={uploading}
          />
        </div>

        <div>
          <label className="block text-sm font-medium text-gray-700 mb-2">
            FASTQ R1 File (paired-end read 1)
          </label>
          <div
            {...dropzoneR1.getRootProps()}
            className={`border-2 border-dashed rounded-lg p-8 text-center cursor-pointer transition ${
              uploading
                ? 'border-gray-200 bg-gray-50 cursor-not-allowed'
                : dropzoneR1.isDragActive
                ? 'border-blue-500 bg-blue-50'
                : 'border-gray-300 hover:border-blue-400'
            }`}
          >
            <input {...dropzoneR1.getInputProps()} disabled={uploading} />
            <Upload className="mx-auto h-12 w-12 text-gray-400 mb-3" />
            {fastqR1 ? (
              <div className="space-y-1">
                <div className="flex items-center justify-center space-x-2">
                  <FileText className="h-5 w-5 text-green-500" />
                  <p className="text-sm text-gray-700">{fastqR1.name}</p>
                </div>
                <p className="text-xs text-gray-500">{formatFileSize(fastqR1.size)}</p>
              </div>
            ) : (
              <p className="text-sm text-gray-600">
                Drop R1 FASTQ file here or click to browse
              </p>
            )}
          </div>
        </div>

        <div>
          <label className="block text-sm font-medium text-gray-700 mb-2">
            FASTQ R2 File (paired-end read 2)
          </label>
          <div
            {...dropzoneR2.getRootProps()}
            className={`border-2 border-dashed rounded-lg p-8 text-center cursor-pointer transition ${
              uploading
                ? 'border-gray-200 bg-gray-50 cursor-not-allowed'
                : dropzoneR2.isDragActive
                ? 'border-blue-500 bg-blue-50'
                : 'border-gray-300 hover:border-blue-400'
            }`}
          >
            <input {...dropzoneR2.getInputProps()} disabled={uploading} />
            <Upload className="mx-auto h-12 w-12 text-gray-400 mb-3" />
            {fastqR2 ? (
              <div className="space-y-1">
                <div className="flex items-center justify-center space-x-2">
                  <FileText className="h-5 w-5 text-green-500" />
                  <p className="text-sm text-gray-700">{fastqR2.name}</p>
                </div>
                <p className="text-xs text-gray-500">{formatFileSize(fastqR2.size)}</p>
              </div>
            ) : (
              <p className="text-sm text-gray-600">
                Drop R2 FASTQ file here or click to browse
              </p>
            )}
          </div>
        </div>

        {uploading && (
          <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
            <div className="flex items-center justify-between mb-2">
              <div className="flex items-center space-x-2">
                <Loader className="h-5 w-5 text-blue-600 animate-spin" />
                <span className="text-sm font-medium text-blue-800">
                  Uploading files...
                </span>
              </div>
              <span className="text-sm font-semibold text-blue-600">
                {uploadProgress}%
              </span>
            </div>
            <div className="w-full bg-blue-200 rounded-full h-2.5">
              <div
                className="bg-blue-600 h-2.5 rounded-full transition-all duration-300"
                style={{ width: `${uploadProgress}%` }}
              ></div>
            </div>
            <p className="text-xs text-blue-600 mt-2">
              Please wait while your files are being uploaded...
            </p>
          </div>
        )}

        {error && (
          <div className="bg-red-50 border border-red-200 text-red-600 p-3 rounded-lg text-sm">
            {error}
          </div>
        )}

        {success && (
          <div className="bg-green-50 border border-green-200 text-green-600 p-3 rounded-lg text-sm">
            {success}
          </div>
        )}

        <button
          type="submit"
          disabled={uploading}
          className="w-full bg-blue-600 text-white py-3 rounded-lg font-semibold hover:bg-blue-700 disabled:bg-blue-300 disabled:cursor-not-allowed transition flex items-center justify-center space-x-2"
        >
          {uploading ? (
            <>
              <Loader className="h-5 w-5 animate-spin" />
              <span>Uploading {uploadProgress}%...</span>
            </>
          ) : (
            <span>Submit Analysis</span>
          )}
        </button>
      </form>
    </div>
  );
}
