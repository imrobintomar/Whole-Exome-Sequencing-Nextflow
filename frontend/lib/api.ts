import axios from 'axios';
import { auth } from './firebase';

const API_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

const api = axios.create({
  baseURL: API_URL,
});

// Add Firebase token to requests
api.interceptors.request.use(async (config) => {
  const user = auth.currentUser;
  if (user) {
    const token = await user.getIdToken();
    config.headers.Authorization = `Bearer ${token}`;
  }
  return config;
});

export interface User {
  uid: string;
  email: string;
  displayName?: string;
}

export interface Job {
  id: number;
  job_id: string;
  sample_name: string;
  status: 'pending' | 'running' | 'completed' | 'failed';
  created_at: string;
  started_at?: string;
  completed_at?: string;
  error_message?: string;
}

export interface JobDetail extends Job {
  bam_path?: string;
  raw_vcf_path?: string;
  annotated_vcf_path?: string;
  filtered_tsv_path?: string;
}

export const authApi = {
  getCurrentUser: () => {
    const user = auth.currentUser;
    if (!user) return null;
    return {
      uid: user.uid,
      email: user.email || '',
      displayName: user.displayName || user.email || '',
    };
  },
};

export const jobApi = {
  submitJob: async (
    sampleName: string,
    fastqR1: File,
    fastqR2: File,
    onProgress?: (progress: number) => void
  ) => {
    const formData = new FormData();
    formData.append('sample_name', sampleName);
    formData.append('fastq_r1', fastqR1);
    formData.append('fastq_r2', fastqR2);

    const response = await api.post<Job>('/jobs/submit', formData, {
      headers: {
        'Content-Type': 'multipart/form-data',
      },
      onUploadProgress: (progressEvent) => {
        if (progressEvent.total) {
          const percentCompleted = Math.round((progressEvent.loaded * 100) / progressEvent.total);
          onProgress?.(percentCompleted);
        }
      },
    });
    return response.data;
  },

  getJobs: async () => {
    const response = await api.get<Job[]>('/jobs');
    return response.data;
  },

  getJob: async (jobId: string) => {
    const response = await api.get<JobDetail>(`/jobs/${jobId}`);
    return response.data;
  },

  downloadFile: (jobId: string, fileType: 'bam' | 'raw_vcf' | 'annotated_vcf' | 'filtered_tsv') => {
    window.open(`${API_URL}/jobs/${jobId}/download/${fileType}`, '_blank');
  },
};

export default api;
