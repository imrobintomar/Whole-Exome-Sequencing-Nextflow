import axios from 'axios';
import { auth } from './firebase';

const API_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

const api = axios.create({
  baseURL: API_URL,
  headers: {
    'ngrok-skip-browser-warning': 'true',
  },
});

// Add Firebase token to requests
api.interceptors.request.use(async (config) => {
  try {
    const user = auth.currentUser;

    if (user) {
      const token = await user.getIdToken(true);
      config.headers.Authorization = `Bearer ${token}`;
    }
  } catch (error) {
    console.error('Error getting auth token:', error);
  }

  config.headers['ngrok-skip-browser-warning'] = 'true';
  return config;
}, (error) => {
  return Promise.reject(error);
});

export interface PhenotypeStatusResponse {
  has_phenotype_analysis: boolean;
  hpo_terms: string[];
  analysis_date?: string;
  file_path?: string;
}

export const phenotypeApi = {
  getPhenotypeStatus: async (jobId: string): Promise<PhenotypeStatusResponse> => {
    const response = await api.get<PhenotypeStatusResponse>(`/jobs/${jobId}/phenotype/status`);
    return response.data;
  },

  runPhenotypeAnalysis: async (jobId: string, hpoTerms: string[]) => {
    const response = await api.post(`/jobs/${jobId}/phenotype/analyze`, {
      hpo_terms: hpoTerms,
    });
    return response.data;
  },

  downloadPhenotypeFile: async (jobId: string) => {
    try {
      const response = await api.get(`/jobs/${jobId}/phenotype/download`, {
        responseType: 'blob',
      });

      // Extract filename from Content-Disposition header
      const contentDisposition = response.headers['content-disposition'];
      let filename = `phenotype_analysis_${jobId}.tsv`;

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
      console.error('Phenotype file download failed:', error);
      throw error;
    }
  },
};

export default api;
