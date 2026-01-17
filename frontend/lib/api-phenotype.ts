import axios from 'axios';
import { auth } from './firebase';

const API_BASE_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

const getAuthHeaders = async () => {
  const user = auth.currentUser;
  if (!user) throw new Error('User not authenticated');
  const token = await user.getIdToken();
  return { Authorization: `Bearer ${token}` };
};

export interface HPOPhenotypeRequest {
  hpo_terms: string[];
}

export interface PhenotypeStatusResponse {
  job_id: string;
  has_phenotype_analysis: boolean;
  hpo_terms: string[];
  augmented_file_path: string | null;
}

export interface PhenotypeVisualizationsResponse {
  job_id: string;
  visualizations: {
    gene_rank: {
      available: boolean;
      path: string | null;
    };
    score_distribution: {
      available: boolean;
      path: string | null;
    };
  };
}

export const phenotypeApi = {
  async runPhenotypeAnalysis(jobId: string, hpoTerms: string[]) {
    const headers = await getAuthHeaders();
    const response = await axios.post(
      `${API_BASE_URL}/jobs/${jobId}/phenotype/analyze`,
      { hpo_terms: hpoTerms },
      { headers }
    );
    return response.data;
  },

  async getPhenotypeStatus(jobId: string): Promise<PhenotypeStatusResponse> {
    const headers = await getAuthHeaders();
    const response = await axios.get(
      `${API_BASE_URL}/jobs/${jobId}/phenotype/status`,
      { headers }
    );
    return response.data;
  },

  async downloadPhenotypeFile(jobId: string) {
    const headers = await getAuthHeaders();
    const response = await axios.get(
      `${API_BASE_URL}/jobs/${jobId}/phenotype/download`,
      { headers, responseType: 'blob' }
    );
    const blob = new Blob([response.data], { type: 'text/tab-separated-values' });
    const url = window.URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = `phenotype_analysis_${jobId}.txt`;
    link.click();
    window.URL.revokeObjectURL(url);
  },

  async getVisualizations(jobId: string): Promise<PhenotypeVisualizationsResponse> {
    const headers = await getAuthHeaders();
    const response = await axios.get(
      `${API_BASE_URL}/jobs/${jobId}/phenotype/visualizations`,
      { headers }
    );
    return response.data;
  },
};
