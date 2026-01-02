'use client';

import { useState } from 'react';
import { useDropzone } from 'react-dropzone';
import { jobApi } from '@/lib/api';
import { Upload, FileText, Loader, CheckCircle2 } from 'lucide-react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from './ui/card';
import { Input } from './ui/input';
import { Label } from './ui/label';
import { Button } from './ui/button';

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
    <div className="space-y-6">
      <div>
        <h2 className="text-3xl font-bold tracking-tight">Submit New Analysis</h2>
        <p className="text-muted-foreground">
          Upload your paired-end FASTQ files for whole exome sequencing analysis
        </p>
      </div>

      <Card className="max-w-3xl">
        <CardHeader>
          <CardTitle>Sample Information</CardTitle>
          <CardDescription>
            Provide a unique name for your sample and upload the paired FASTQ files
          </CardDescription>
        </CardHeader>
        <CardContent>
          <form onSubmit={handleSubmit} className="space-y-6">
            <div className="space-y-2">
              <Label htmlFor="sampleName">Sample Name</Label>
              <Input
                id="sampleName"
                type="text"
                value={sampleName}
                onChange={(e) => setSampleName(e.target.value)}
                placeholder="e.g., Sample001"
                required
                disabled={uploading}
              />
            </div>

            <div className="space-y-2">
              <Label>FASTQ R1 File (paired-end read 1)</Label>
              <div
                {...dropzoneR1.getRootProps()}
                className={`border-2 border-dashed rounded-lg p-8 text-center cursor-pointer transition-colors ${
                  uploading
                    ? 'border-muted bg-muted/10 cursor-not-allowed'
                    : dropzoneR1.isDragActive
                    ? 'border-primary bg-primary/5'
                    : 'border-muted-foreground/25 hover:border-primary/50 hover:bg-accent/50'
                }`}
              >
                <input {...dropzoneR1.getInputProps()} disabled={uploading} />
                {fastqR1 ? (
                  <div className="space-y-2">
                    <div className="flex items-center justify-center gap-2">
                      <CheckCircle2 className="h-5 w-5 text-green-500" />
                      <FileText className="h-5 w-5 text-primary" />
                    </div>
                    <p className="font-medium">{fastqR1.name}</p>
                    <p className="text-sm text-muted-foreground">{formatFileSize(fastqR1.size)}</p>
                  </div>
                ) : (
                  <>
                    <Upload className="mx-auto h-12 w-12 text-muted-foreground mb-3" />
                    <p className="text-sm text-muted-foreground">
                      Drop R1 FASTQ file here or click to browse
                    </p>
                    <p className="text-xs text-muted-foreground mt-1">
                      Accepts .gz, .fastq.gz, .fq.gz files
                    </p>
                  </>
                )}
              </div>
            </div>

            <div className="space-y-2">
              <Label>FASTQ R2 File (paired-end read 2)</Label>
              <div
                {...dropzoneR2.getRootProps()}
                className={`border-2 border-dashed rounded-lg p-8 text-center cursor-pointer transition-colors ${
                  uploading
                    ? 'border-muted bg-muted/10 cursor-not-allowed'
                    : dropzoneR2.isDragActive
                    ? 'border-primary bg-primary/5'
                    : 'border-muted-foreground/25 hover:border-primary/50 hover:bg-accent/50'
                }`}
              >
                <input {...dropzoneR2.getInputProps()} disabled={uploading} />
                {fastqR2 ? (
                  <div className="space-y-2">
                    <div className="flex items-center justify-center gap-2">
                      <CheckCircle2 className="h-5 w-5 text-green-500" />
                      <FileText className="h-5 w-5 text-primary" />
                    </div>
                    <p className="font-medium">{fastqR2.name}</p>
                    <p className="text-sm text-muted-foreground">{formatFileSize(fastqR2.size)}</p>
                  </div>
                ) : (
                  <>
                    <Upload className="mx-auto h-12 w-12 text-muted-foreground mb-3" />
                    <p className="text-sm text-muted-foreground">
                      Drop R2 FASTQ file here or click to browse
                    </p>
                    <p className="text-xs text-muted-foreground mt-1">
                      Accepts .gz, .fastq.gz, .fq.gz files
                    </p>
                  </>
                )}
              </div>
            </div>

            {uploading && (
              <Card className="border-primary/50 bg-primary/5">
                <CardContent className="pt-6">
                  <div className="flex items-center justify-between mb-2">
                    <div className="flex items-center gap-2">
                      <Loader className="h-5 w-5 text-primary animate-spin" />
                      <span className="text-sm font-medium">
                        Uploading files...
                      </span>
                    </div>
                    <span className="text-sm font-semibold text-primary">
                      {uploadProgress}%
                    </span>
                  </div>
                  <div className="w-full bg-muted rounded-full h-2.5">
                    <div
                      className="bg-primary h-2.5 rounded-full transition-all duration-300"
                      style={{ width: `${uploadProgress}%` }}
                    ></div>
                  </div>
                  <p className="text-xs text-muted-foreground mt-2">
                    Please wait while your files are being uploaded...
                  </p>
                </CardContent>
              </Card>
            )}

            {error && (
              <Card className="border-destructive/50 bg-destructive/5">
                <CardContent className="pt-6">
                  <p className="text-sm text-destructive">{error}</p>
                </CardContent>
              </Card>
            )}

            {success && (
              <Card className="border-green-500/50 bg-green-500/5">
                <CardContent className="pt-6">
                  <div className="flex items-center gap-2">
                    <CheckCircle2 className="h-5 w-5 text-green-600" />
                    <p className="text-sm text-green-600">{success}</p>
                  </div>
                </CardContent>
              </Card>
            )}

            <Button
              type="submit"
              disabled={uploading}
              className="w-full h-11"
            >
              {uploading ? (
                <>
                  <Loader className="mr-2 h-5 w-5 animate-spin" />
                  Uploading {uploadProgress}%...
                </>
              ) : (
                <>
                  <Upload className="mr-2 h-5 w-5" />
                  Submit Analysis
                </>
              )}
            </Button>
          </form>
        </CardContent>
      </Card>
    </div>
  );
}
