import axios from 'axios';
import { auth } from './firebase';

const API_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

const api = axios.create({
  baseURL: API_URL,
  headers: {
    // Add ngrok-skip-browser-warning header for all requests
    'ngrok-skip-browser-warning': 'true',
  },
});

// Add Firebase token to requests
api.interceptors.request.use(async (config) => {
  try {
    // Wait a bit for auth state to be ready after login
    const user = auth.currentUser;

    if (user) {
      // Force refresh to ensure we have a valid token
      const token = await user.getIdToken(true);
      config.headers.Authorization = `Bearer ${token}`;
    }
  } catch (error) {
    console.error('Error getting auth token:', error);
  }

  // Ensure ngrok header is always present
  config.headers['ngrok-skip-browser-warning'] = 'true';
  return config;
}, (error) => {
  return Promise.reject(error);
});

export interface User {
  uid: string
  email: string
  displayName?: string
  photoURL?: string 
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

  deleteJob: async (jobId: string) => {
    const response = await api.delete(`/jobs/${jobId}`);
    return response.data;
  },

  checkFileAvailability: async (jobId: string, fileType: 'bam' | 'vcf') => {
    try {
      const response = await api.get(`/jobs/${jobId}/files/${fileType}`);
      return response.data;
    } catch (error) {
      return null;
    }
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
export const deleteJob = jobApi.deleteJob;

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

export interface GenePanelResult {
  job_id: string;
  sample_name: string;
  applied_genes: string[];
  gene_count: number;
  total_variants: number;
  showing_variants: number;
  variants: Record<string, any>[];
  statistics: {
    chromosome_distribution: Record<string, number>;
    gene_distribution: Record<string, number>;
    functional_impact: Record<string, number>;
    exonic_subcategories: Record<string, number>;
    clinical_significance: Record<string, number>;
  };
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

  applyGenePanel: async (jobId: string, genes: string[]) => {
    const response = await api.post<GenePanelResult>(`/jobs/${jobId}/apply-panel`, {
      genes: genes
    });
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

// SaaS Extension - Admin API
export interface AdminDashboardStats {
  users: {
    total: number;
    active_subscriptions: number;
  };
  jobs: {
    total: number;
    running: number;
  };
  revenue: {
    monthly_recurring_cents: number;
    monthly_recurring_usd: number;
  };
  system: {
    cpu: { usage_percent: number; count: number };
    memory: { total_gb: number; used_gb: number; percent: number };
    disk: { total_gb: number; used_gb: number; percent: number };
  };
  nextflow_processes: number;
  unread_chats: number;
}

export interface AdminUser {
  uid: string;
  email: string;
  created_at: string;
  subscription: {
    plan: string;
    status: string;
  } | null;
  usage: {
    jobs_executed: number;
    jobs_limit: number;
  } | null;
  is_active?: boolean;
  is_banned?: boolean;
  ban_reason?: string;
  banned_at?: string;
  banned_by?: string;
}

export interface AdminJob extends Job {
  user_id: string;
}

export const adminApi = {
  getDashboardStats: async () => {
    const response = await api.get<AdminDashboardStats>('/admin/dashboard/stats');
    return response.data;
  },

  getAllJobs: async (params?: { status?: string; user_id?: string; search?: string; date_from?: string; date_to?: string; limit?: number; offset?: number }) => {
    const response = await api.get<{ jobs: AdminJob[]; total: number; limit: number; offset: number }>('/admin/jobs', { params });
    return response.data;
  },

  bulkJobAction: async (action: 'cancel' | 'delete', jobIds: string[]) => {
    const response = await api.post<{ success: string[]; failed: Array<{ job_id: string; error: string }> }>('/admin/jobs/bulk-action', {
      action,
      job_ids: jobIds
    });
    return response.data;
  },

  exportUsers: async () => {
    const response = await api.get('/admin/users/export', {
      responseType: 'blob'
    });

    // Download the CSV file
    const url = window.URL.createObjectURL(new Blob([response.data]));
    const link = document.createElement('a');
    link.href = url;
    link.download = `users_export_${new Date().toISOString().split('T')[0]}.csv`;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    window.URL.revokeObjectURL(url);
  },

  getJobLogs: async (jobId: string) => {
    const response = await api.get<{ logs: string | null; path: string; error?: string }>(`/admin/jobs/${jobId}/logs`);
    return response.data;
  },

  getAllUsers: async () => {
    const response = await api.get<{ users: AdminUser[]; total: number }>('/admin/users');
    return response.data;
  },

  getSystemMetrics: async () => {
    const response = await api.get<{
      system: AdminDashboardStats['system'];
      nextflow_processes: number;
      storage: { total_gb: number; used_gb: number; percent: number };
      timestamp: string;
    }>('/admin/system/metrics');
    return response.data;
  },

  updateUserUsage: async (userUid: string, jobsLimit?: number, resetCount?: boolean) => {
    const response = await api.post<{
      success: boolean;
      usage: {
        user_id: string;
        jobs_executed: number;
        jobs_limit: number;
        month: number;
      };
    }>(`/admin/users/${userUid}/usage`, null, {
      params: { jobs_limit: jobsLimit, reset_count: resetCount }
    });
    return response.data;
  },

  updateUserSubscription: async (userUid: string, planName: string) => {
    const response = await api.post<{
      success: boolean;
      subscription: {
        id: number;
        user_id: string;
        plan_name: string;
        plan_id: number;
        status: string;
        monthly_jobs_limit: number;
        price_usd: number;
        current_period_start: string;
        current_period_end: string;
      };
    }>(`/admin/users/${userUid}/subscription`, null, {
      params: { plan_name: planName }
    });
    return response.data;
  },

  banUser: async (userUid: string, reason: string) => {
    const response = await api.post<{ success: boolean; message: string }>(`/admin/users/${userUid}/ban`, null, {
      params: { reason }
    });
    return response.data;
  },

  unbanUser: async (userUid: string) => {
    const response = await api.post<{ success: boolean; message: string }>(`/admin/users/${userUid}/unban`);
    return response.data;
  },

  suspendUser: async (userUid: string) => {
    const response = await api.post<{ success: boolean; message: string }>(`/admin/users/${userUid}/suspend`);
    return response.data;
  },

  activateUser: async (userUid: string) => {
    const response = await api.post<{ success: boolean; message: string }>(`/admin/users/${userUid}/activate`);
    return response.data;
  },

  // User Details endpoints
  getUserDetails: async (userUid: string) => {
    const response = await api.get<UserDetailsResponse>(`/admin/users/${userUid}/details`);
    return response.data;
  },

  addUserNote: async (userUid: string, noteText: string) => {
    const response = await api.post<{ success: boolean; note: UserNote }>(`/admin/users/${userUid}/notes`, null, {
      params: { note_text: noteText }
    });
    return response.data;
  },

  addUserTag: async (userUid: string, tagName: string, color: string = 'blue') => {
    const response = await api.post<{ success: boolean; tag: UserTag }>(`/admin/users/${userUid}/tags`, null, {
      params: { tag_name: tagName, color }
    });
    return response.data;
  },

  removeUserTag: async (userUid: string, tagId: number) => {
    const response = await api.delete<{ success: boolean; message: string }>(`/admin/users/${userUid}/tags/${tagId}`);
    return response.data;
  },

  // Email Notification endpoints
  testEmailConnection: async () => {
    const response = await api.post<{ status: string; message: string; smtp_host: string; smtp_port: number; smtp_user: string; admin_email: string }>('/admin/email/test');
    return response.data;
  },

  sendCustomEmail: async (userEmail: string, subject: string, message: string, userName?: string) => {
    const response = await api.post<{ success: boolean; message: string }>('/admin/email/send-custom', null, {
      params: { user_email: userEmail, subject, message, user_name: userName }
    });
    return response.data;
  },

  sendPaymentReminder: async (userUid: string, amount: number, dueDate: string, invoiceUrl?: string) => {
    const response = await api.post<{ success: boolean; message: string }>('/admin/email/payment-reminder', null, {
      params: { user_uid: userUid, amount, due_date: dueDate, invoice_url: invoiceUrl }
    });
    return response.data;
  },

  sendSubscriptionExpiry: async (userUid: string) => {
    const response = await api.post<{ success: boolean; message: string }>('/admin/email/subscription-expiry', null, {
      params: { user_uid: userUid }
    });
    return response.data;
  },

  getHealthAlerts: async (limit: number = 50, severity?: string) => {
    const response = await api.get<{ alerts: any[]; count: number }>('/admin/email/health-alerts', {
      params: { limit, severity }
    });
    return response.data;
  },

  resolveHealthAlert: async (alertId: number) => {
    const response = await api.post<{ success: boolean; message: string }>(`/admin/email/health-alerts/${alertId}/resolve`);
    return response.data;
  },

  sendTestHealthAlert: async () => {
    const response = await api.post<{ success: boolean; message: string }>('/admin/email/health-alerts/test');
    return response.data;
  },

  getEmailSettings: async () => {
    const response = await api.get<{
      smtp_host: string;
      smtp_port: number;
      smtp_user: string;
      admin_email: string;
      health_alerts_enabled: boolean;
      health_thresholds: {
        cpu: number;
        memory: number;
        disk: number;
      };
    }>('/admin/email/settings');
    return response.data;
  },

  // Analytics endpoints
  getAnalyticsSummary: async (period: 'day' | 'week' | 'month' | 'year' = 'month') => {
    const response = await api.get<{
      users: {
        total: number;
        new_this_period: number;
        active_subscriptions: number;
      };
      revenue: {
        mrr: number;
        total_subscriptions: number;
      };
      jobs: {
        total: number;
        this_period: number;
        completed: number;
        failed: number;
        success_rate: number;
      };
      period: string;
      start_date: string;
      end_date: string;
    }>('/admin/analytics/summary', {
      params: { period }
    });
    return response.data;
  },

  getRevenueAnalytics: async (fromDate?: string, toDate?: string) => {
    const response = await api.get<{
      revenue_over_time: Array<{ month: string; revenue: number }>;
      revenue_by_plan: Array<{ plan: string; revenue: number }>;
      current_mrr: number;
      arpu: number;
      total_revenue: number;
      from_date: string;
      to_date: string;
    }>('/admin/analytics/revenue', {
      params: { from_date: fromDate, to_date: toDate }
    });
    return response.data;
  },

  getUserAnalytics: async (fromDate?: string, toDate?: string) => {
    const response = await api.get<{
      user_growth: Array<{ month: string; new_users: number; total_users: number }>;
      users_by_plan: Array<{ plan: string; count: number }>;
      total_users: number;
      from_date: string;
      to_date: string;
    }>('/admin/analytics/users', {
      params: { from_date: fromDate, to_date: toDate }
    });
    return response.data;
  },

  getJobAnalytics: async (fromDate?: string, toDate?: string) => {
    const response = await api.get<{
      jobs_over_time: Array<{ date: string; completed: number; failed: number; running: number; pending: number }>;
      jobs_by_status: Array<{ status: string; count: number }>;
      success_rate: number;
      avg_duration_hours: number;
      total_jobs: number;
      from_date: string;
      to_date: string;
    }>('/admin/analytics/jobs', {
      params: { from_date: fromDate, to_date: toDate }
    });
    return response.data;
  },
};

// User Details interfaces
export interface UserDetailsResponse {
  user: {
    uid: string;
    email: string;
    username: string;
    created_at: string;
    is_active: boolean;
    is_banned: boolean;
    ban_reason: string | null;
    banned_at: string | null;
  };
  subscription: {
    plan: string;
    status: string;
    stripe_customer_id: string | null;
    current_period_start: string | null;
    current_period_end: string | null;
    price_cents: number;
    monthly_jobs_limit: number;
  };
  usage: {
    jobs_executed: number;
    jobs_limit: number;
    month: number;
  };
  jobs: Array<{
    job_id: string;
    sample_name: string;
    status: string;
    current_step: string | null;
    error_message: string | null;
    created_at: string;
    updated_at: string;
    completed_at: string | null;
  }>;
  payment_history: Array<{
    id: number;
    event_type: string;
    created_at: string;
    processed: boolean;
    amount: number | null;
    currency: string;
  }>;
  support_tickets: Array<{
    id: number;
    subject: string;
    status: string;
    job_id: string | null;
    message_count: number;
    last_message_at: string;
    created_at: string;
  }>;
  activity_timeline: Array<{
    id: number;
    action: string;
    resource_type: string | null;
    resource_id: string | null;
    created_at: string;
    metadata: string | null;
  }>;
  notes: UserNote[];
  tags: UserTag[];
}

export interface UserNote {
  id: number;
  note_text: string;
  admin_id: string;
  created_at: string;
  updated_at: string;
}

export interface UserTag {
  id: number;
  tag_name: string;
  color: string;
  created_at: string;
  created_by: string;
}

// SaaS Extension - Billing API
export interface SubscriptionPlan {
  id: number;
  name: string;
  price_cents: number;
  price_usd: number;
  monthly_jobs_limit: number;
  chat_support: boolean;
  features: string[];
}

export interface UserSubscription {
  subscription: {
    status: string;
    current_period_start: string;
    current_period_end: string;
    cancel_at_period_end: boolean;
    stripe_customer_id: string;
    stripe_subscription_id: string;
  } | null;
  plan: {
    name: string;
    monthly_jobs_limit: number;
    chat_support: boolean;
    price_cents: number;
    features?: string[];
  };
  usage: {
    jobs_executed: number;
    jobs_limit: number;
    month: number;
  };
}

export interface UsageStats {
  month: number;
  jobs_executed: number;
  jobs_limit: number;
  jobs_remaining: number;
  usage_percent: number;
  plan_name: string;
  chat_support_enabled: boolean;
}

export const billingApi = {
  getPlans: async () => {
    const response = await api.get<{ plans: SubscriptionPlan[] }>('/billing/plans');
    return response.data.plans;
  },

  getSubscription: async () => {
    const response = await api.get<UserSubscription>('/billing/subscription');
    return response.data;
  },

  createCheckout: async (planId: number) => {
    const response = await api.post<{ checkout_url: string; plan: { name: string; price_cents: number } }>('/billing/checkout', null, {
      params: { plan_id: planId }
    });
    return response.data;
  },

  createPortalSession: async () => {
    const response = await api.post<{ portal_url: string }>('/billing/portal');
    return response.data;
  },

  getUsage: async () => {
    const response = await api.get<UsageStats>('/billing/usage');
    return response.data;
  },
};

// SaaS Extension - Chat API
export interface ChatConversation {
  id: number;
  subject: string;
  status: 'open' | 'resolved' | 'closed';
  created_at: string;
  last_message_at: string;
  admin_id: string | null;
}

export interface ChatMessage {
  id: number;
  sender_role: 'user' | 'admin';
  sender_id?: string;
  message: string;
  created_at: string;
}

export interface ConversationDetail {
  conversation: ChatConversation;
  messages: ChatMessage[];
}

export const chatApi = {
  createConversation: async (subject: string, initialMessage: string) => {
    const response = await api.post<{ conversation_id: number; subject: string; status: string; created_at: string }>('/chat/conversations', {
      subject,
      initial_message: initialMessage,
    });
    return response.data;
  },

  getConversations: async () => {
    const response = await api.get<{ conversations: ChatConversation[] }>('/chat/conversations');
    return response.data.conversations;
  },

  getConversationMessages: async (conversationId: number) => {
    const response = await api.get<ConversationDetail>(`/chat/conversations/${conversationId}/messages`);
    return response.data;
  },

  sendMessage: async (conversationId: number, message: string) => {
    const response = await api.post<{ message_id: number; created_at: string }>(`/chat/conversations/${conversationId}/messages`, {
      message,
    });
    return response.data;
  },

  closeConversation: async (conversationId: number) => {
    const response = await api.patch<{ status: string }>(`/chat/conversations/${conversationId}/close`);
    return response.data;
  },

  // Admin endpoints
  adminGetConversations: async (params?: { status?: string; limit?: number; offset?: number }) => {
    const response = await api.get<{
      conversations: Array<ChatConversation & { user_email: string }>;
      total: number;
      limit: number;
      offset: number;
    }>('/chat/admin/conversations', { params });
    return response.data;
  },

  adminGetConversationMessages: async (conversationId: number) => {
    const response = await api.get<ConversationDetail & { conversation: ChatConversation & { user_email: string } }>(`/chat/admin/conversations/${conversationId}/messages`);
    return response.data;
  },

  adminSendMessage: async (conversationId: number, message: string) => {
    const response = await api.post<{ message_id: number; created_at: string }>(`/chat/admin/conversations/${conversationId}/messages`, {
      message,
    });
    return response.data;
  },

  adminUpdateStatus: async (conversationId: number, status: 'open' | 'resolved' | 'closed') => {
    const response = await api.patch<{ status: string }>(`/chat/admin/conversations/${conversationId}/status`, null, {
      params: { status }
    });
    return response.data;
  },
};

export default api;
