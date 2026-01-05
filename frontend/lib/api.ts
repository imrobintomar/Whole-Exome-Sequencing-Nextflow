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
  current_step?: string;
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

  downloadFile: async (jobId: string, fileType: 'bam' | 'raw_vcf' | 'annotated_vcf' | 'filtered_tsv') => {
    try {
      const response = await api.get(`/jobs/${jobId}/download/${fileType}`, {
        responseType: 'blob',
      });

      // Extract filename from Content-Disposition header
      const contentDisposition = response.headers['content-disposition'];
      let filename = `download_${fileType}`;

      if (contentDisposition) {
        // Try multiple regex patterns for filename extraction
        const patterns = [
          /filename\*?=["']?(?:UTF-\d['"])?([^;\r\n"']+)["']?/i,
          /filename=["']([^"']+)["']/i,
          /filename=([^;\r\n]+)/i
        ];

        for (const pattern of patterns) {
          const match = contentDisposition.match(pattern);
          if (match && match[1]) {
            filename = match[1].trim();
            break;
          }
        }
      }

      // Create blob URL and trigger download
      const blob = new Blob([response.data]);
      const url = window.URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = url;
      link.download = filename;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      window.URL.revokeObjectURL(url);
    } catch (error) {
      console.error('Download failed:', error);
      throw error;
    }
  },

  // Pipeline control actions
  cancelJob: async (jobId: string) => {
    const response = await api.post(`/jobs/${jobId}/cancel`);
    return response.data;
  },

  rerunJob: async (jobId: string) => {
    const response = await api.post(`/jobs/${jobId}/rerun`);
    return response.data;
  },

  resumeJob: async (jobId: string) => {
    const response = await api.post(`/jobs/${jobId}/resume`);
    return response.data;
  },
};

// Convenience exports for direct access
export const getJobs = jobApi.getJobs;
export const getJob = jobApi.getJob;
export const submitJob = jobApi.submitJob;
export const downloadFile = jobApi.downloadFile;
export const cancelJob = jobApi.cancelJob;
export const rerunJob = jobApi.rerunJob;
export const resumeJob = jobApi.resumeJob;

// ACMG Classification API
export interface ACMGClassification {
  classification: 'Pathogenic' | 'Likely Pathogenic' | 'Uncertain Significance' | 'Likely Benign' | 'Benign';
  evidence: Record<string, {
    met: boolean;
    strength: string;
    description: string;
    auto_applied: boolean;
  }>;
  evidence_summary: {
    PVS: number;
    PS: number;
    PM: number;
    PP: number;
    BA: number;
    BS: number;
    BP: number;
  };
  met_criteria: string[];
}

export interface JobClassificationResult {
  job_id: string;
  sample_name: string;
  total_variants: number;
  classifications: Array<{
    position: string;
    gene: string;
    consequence: string;
    classification: string;
    evidence_summary: Record<string, number>;
    met_criteria: string[];
  }>;
  summary: {
    pathogenic: number;
    likely_pathogenic: number;
    vus: number;
    likely_benign: number;
    benign: number;
  };
}

export const acmgApi = {
  classifySingleVariant: async (variant: any) => {
    const response = await api.post<ACMGClassification>('/classify/acmg', variant);
    return response.data;
  },

  classifyJobVariants: async (jobId: string) => {
    const response = await api.post<JobClassificationResult>(`/jobs/${jobId}/classify`);
    return response.data;
  }
};

// Gene Panel API
export interface GenePanel {
  id: number;
  name: string;
  disease_group: string;
  disease_sub_group: string;
  relevant_disorders: string[];
  version: string;
  genes_count: number;
}

export interface GenePanelGenes {
  panel_id: number;
  genes: string[];
  count: number;
}

export const panelApi = {
  searchPanels: async (query: string) => {
    const response = await api.get<{ results: GenePanel[] }>('/panels/search', {
      params: { query }
    });
    return response.data.results;
  },

  getPanelGenes: async (panelId: number, confidenceLevel: string = '3') => {
    const response = await api.get<GenePanelGenes>(`/panels/${panelId}/genes`, {
      params: { confidence_level: confidenceLevel }
    });
    return response.data;
  },

  getAcmgSecondaryFindings: async () => {
    const response = await api.get<{ genes: string[]; count: number; version: string }>('/panels/acmg-sf');
    return response.data;
  },

  downloadFilteredVariants: async (jobId: string, genes: string[]) => {
    try {
      const response = await api.post(`/jobs/${jobId}/download/filtered`, genes, {
        responseType: 'blob',
      });

      // Extract filename from Content-Disposition header
      const contentDisposition = response.headers['content-disposition'];
      let filename = `filtered_variants_${genes.length}genes.tsv`;

      if (contentDisposition) {
        const patterns = [
          /filename\*?=["']?(?:UTF-\d['"])?([^;\r\n"']+)["']?/i,
          /filename=["']([^"']+)["']/i,
          /filename=([^;\r\n]+)/i
        ];

        for (const pattern of patterns) {
          const match = contentDisposition.match(pattern);
          if (match && match[1]) {
            filename = match[1].trim();
            break;
          }
        }
      }

      // Create blob URL and trigger download
      const blob = new Blob([response.data]);
      const url = window.URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = url;
      link.download = filename;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      window.URL.revokeObjectURL(url);
    } catch (error) {
      console.error('Filtered download failed:', error);
      throw error;
    }
  }
};

// Variant Analysis API
export interface VariantMetrics {
  job_id: string;
  sample_name: string;
  metrics: {
    headline: {
      total_variants: number;
      total_genes: number;
      snvs: number;
      indels: number;
      high_confidence: number;
    };
    chromosome_distribution: Array<{
      chromosome: string;
      count: number;
      type: 'autosome' | 'sex';
    }>;
    chromosome_distribution_normalized: Array<{
      chromosome: string;
      count: number;
      type: 'autosome' | 'sex';
    }>;
    gene_distribution: Array<{
      gene: string;
      count: number;
    }>;
    gene_distribution_protein_altering: Array<{
      gene: string;
      count: number;
    }>;
    functional_impact: {
      categories: Record<string, number>;
      exonic_subcategories: Record<string, number>;
    };
  };
}

export const variantApi = {
  getVariantMetrics: async (jobId: string) => {
    const response = await api.get<VariantMetrics>(`/jobs/${jobId}/variant-metrics`);
    return response.data;
  }
};

export default api;
