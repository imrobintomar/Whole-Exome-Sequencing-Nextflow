'use client';

import { useState, useEffect } from 'react';
import { adminApi, AdminDashboardStats, AdminUser, AdminJob } from '../lib/api';
import { auth } from '../lib/firebase';
import { onAuthStateChanged, signOut } from 'firebase/auth';
import AdminLogin from './AdminLogin';
import AdminChatPanel from './AdminChatPanel';

export default function AdminDashboard() {
  const [stats, setStats] = useState<AdminDashboardStats | null>(null);
  const [users, setUsers] = useState<AdminUser[]>([]);
  const [jobs, setJobs] = useState<AdminJob[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [activeTab, setActiveTab] = useState<'overview' | 'users' | 'jobs' | 'system' | 'chat'>('overview');
  const [authReady, setAuthReady] = useState(false);
  const [isAuthenticated, setIsAuthenticated] = useState(false);

  // Wait for Firebase auth to be ready
  useEffect(() => {
    const unsubscribe = onAuthStateChanged(auth, (user) => {
      if (user) {
        setAuthReady(true);
        setIsAuthenticated(true);
      } else {
        setIsAuthenticated(false);
        setLoading(false);
      }
    });
    return () => unsubscribe();
  }, []);

  useEffect(() => {
    if (!authReady) return;

    loadDashboardData();
    // Refresh every 30 seconds
    const interval = setInterval(loadDashboardData, 30000);
    return () => clearInterval(interval);
  }, [authReady]);

  const loadDashboardData = async () => {
    try {
      const [dashboardStats, usersData, jobsData] = await Promise.all([
        adminApi.getDashboardStats(),
        adminApi.getAllUsers(),
        adminApi.getAllJobs({ limit: 20 }),
      ]);

      setStats(dashboardStats);
      setUsers(usersData.users);
      setJobs(jobsData.jobs);
      setError(null);
    } catch (err: any) {
      setError(err.response?.status === 403 ? 'Access denied. Admin privileges required.' : 'Failed to load dashboard data');
      console.error('Dashboard error:', err);
    } finally {
      setLoading(false);
    }
  };

  const handleUpdateUsage = async (userUid: string, currentLimit: number, action: 'increase' | 'decrease' | 'reset') => {
    try {
      let newLimit = currentLimit;
      let resetCount = false;

      if (action === 'increase') {
        newLimit = currentLimit + 1;
      } else if (action === 'decrease') {
        newLimit = Math.max(0, currentLimit - 1);
      } else if (action === 'reset') {
        resetCount = true;
      }

      await adminApi.updateUserUsage(userUid, action === 'reset' ? undefined : newLimit, resetCount);

      // Reload users data
      const usersData = await adminApi.getAllUsers();
      setUsers(usersData.users);
    } catch (err: any) {
      console.error('Error updating usage:', err);
      alert('Failed to update user usage');
    }
  };

  const handleLogout = async () => {
    try {
      await signOut(auth);
    } catch (err) {
      console.error('Logout error:', err);
    }
  };

  // Show login page if not authenticated
  if (!isAuthenticated && !loading) {
    return <AdminLogin />;
  }

  if (loading) {
    return (
      <div className="flex items-center justify-center min-h-screen">
        <div className="text-center">
          <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-600 mx-auto mb-4"></div>
          <p className="text-gray-600">Loading admin dashboard...</p>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="flex items-center justify-center min-h-screen">
        <div className="bg-red-50 border border-red-200 rounded-lg p-6 max-w-md">
          <h3 className="text-red-800 font-semibold mb-2">Error</h3>
          <p className="text-red-600">{error}</p>
          <button
            onClick={loadDashboardData}
            className="mt-4 px-4 py-2 bg-red-600 text-white rounded hover:bg-red-700"
          >
            Retry
          </button>
        </div>
      </div>
    );
  }

  if (!stats) return null;

  return (
    <div className="min-h-screen bg-gray-50">
      {/* Header */}
      <div className="bg-white border-b">
        <div className="max-w-7xl mx-auto px-4 py-6">
          <div className="flex justify-between items-center">
            <div>
              <h1 className="text-3xl font-bold text-gray-900">Admin Dashboard</h1>
              <p className="text-gray-600 mt-1">Platform management and monitoring</p>
            </div>
            <button
              onClick={handleLogout}
              className="px-4 py-2 bg-gray-600 text-white rounded-lg hover:bg-gray-700 transition-colors"
            >
              Logout
            </button>
          </div>
        </div>
      </div>

      <div className="max-w-7xl mx-auto px-4 py-6">
        {/* Tabs */}
        <div className="border-b border-gray-200 mb-6">
          <nav className="-mb-px flex space-x-8">
            {[
              { id: 'overview', label: 'Overview' },
              { id: 'users', label: 'Users' },
              { id: 'jobs', label: 'Jobs' },
              { id: 'system', label: 'System' },
              { id: 'chat', label: 'Support Chat' },
            ].map((tab) => (
              <button
                key={tab.id}
                onClick={() => setActiveTab(tab.id as any)}
                className={`${
                  activeTab === tab.id
                    ? 'border-blue-500 text-blue-600'
                    : 'border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300'
                } whitespace-nowrap py-4 px-1 border-b-2 font-medium text-sm`}
              >
                {tab.label}
              </button>
            ))}
          </nav>
        </div>

        {/* Overview Tab */}
        {activeTab === 'overview' && (
          <div className="space-y-6">
            {/* Stats Cards */}
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4">
              <StatCard
                title="Total Users"
                value={stats.users.total}
                subtitle={`${stats.users.active_subscriptions} active subscriptions`}
                icon="👥"
                color="blue"
              />
              <StatCard
                title="Total Jobs"
                value={stats.jobs.total}
                subtitle={`${stats.jobs.running} currently running`}
                icon="🔬"
                color="green"
              />
              <StatCard
                title="Monthly Revenue"
                value={`$${stats.revenue.monthly_recurring_usd.toFixed(2)}`}
                subtitle="MRR"
                icon="💰"
                color="purple"
              />
              <StatCard
                title="Unread Chats"
                value={stats.unread_chats}
                subtitle="Support requests"
                icon="💬"
                color="orange"
              />
            </div>

            {/* System Status */}
            <div className="bg-white rounded-lg shadow p-6">
              <h3 className="text-lg font-semibold mb-4">System Health</h3>
              <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
                <SystemMetric
                  label="CPU Usage"
                  value={stats.system.cpu.usage_percent}
                  suffix="%"
                  max={100}
                  color={stats.system.cpu.usage_percent > 80 ? 'red' : stats.system.cpu.usage_percent > 60 ? 'yellow' : 'green'}
                />
                <SystemMetric
                  label="Memory Usage"
                  value={stats.system.memory.percent}
                  suffix="%"
                  max={100}
                  color={stats.system.memory.percent > 80 ? 'red' : stats.system.memory.percent > 60 ? 'yellow' : 'green'}
                />
                <SystemMetric
                  label="Disk Usage"
                  value={stats.system.disk.percent}
                  suffix="%"
                  max={100}
                  color={stats.system.disk.percent > 80 ? 'red' : stats.system.disk.percent > 60 ? 'yellow' : 'green'}
                />
              </div>
              <div className="mt-4 pt-4 border-t">
                <p className="text-sm text-gray-600">
                  Nextflow Processes: <span className="font-semibold text-gray-900">{stats.nextflow_processes}</span>
                </p>
              </div>
            </div>
          </div>
        )}

        {/* Users Tab */}
        {activeTab === 'users' && (
          <div className="bg-white rounded-lg shadow overflow-hidden">
            <div className="px-6 py-4 border-b border-gray-200">
              <h3 className="text-lg font-semibold">All Users ({users.length})</h3>
            </div>
            <div className="overflow-x-auto">
              <table className="min-w-full divide-y divide-gray-200">
                <thead className="bg-gray-50">
                  <tr>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Email</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Plan</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Status</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Usage</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Created</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Actions</th>
                  </tr>
                </thead>
                <tbody className="bg-white divide-y divide-gray-200">
                  {users.map((user) => (
                    <tr key={user.uid} className="hover:bg-gray-50">
                      <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-900">{user.email}</td>
                      <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-900">
                        {user.subscription?.plan || 'Free'}
                      </td>
                      <td className="px-6 py-4 whitespace-nowrap">
                        <span
                          className={`px-2 py-1 text-xs font-medium rounded-full ${
                            user.subscription?.status === 'active'
                              ? 'bg-green-100 text-green-800'
                              : 'bg-gray-100 text-gray-800'
                          }`}
                        >
                          {user.subscription?.status || 'none'}
                        </span>
                      </td>
                      <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-600">
                        {user.usage ? `${user.usage.jobs_executed}/${user.usage.jobs_limit}` : '0/2'}
                      </td>
                      <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-600">
                        {new Date(user.created_at).toLocaleDateString()}
                      </td>
                      <td className="px-6 py-4 whitespace-nowrap text-sm font-medium">
                        <div className="flex space-x-2">
                          <button
                            onClick={() => handleUpdateUsage(user.uid, user.usage?.jobs_limit || 2, 'increase')}
                            className="text-green-600 hover:text-green-900"
                            title="Increase limit"
                          >
                            +
                          </button>
                          <button
                            onClick={() => handleUpdateUsage(user.uid, user.usage?.jobs_limit || 2, 'decrease')}
                            className="text-red-600 hover:text-red-900"
                            title="Decrease limit"
                          >
                            -
                          </button>
                          <button
                            onClick={() => handleUpdateUsage(user.uid, user.usage?.jobs_limit || 2, 'reset')}
                            className="text-blue-600 hover:text-blue-900 text-xs"
                            title="Reset usage count"
                          >
                            Reset
                          </button>
                        </div>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        )}

        {/* Jobs Tab */}
        {activeTab === 'jobs' && (
          <div className="bg-white rounded-lg shadow overflow-hidden">
            <div className="px-6 py-4 border-b border-gray-200">
              <h3 className="text-lg font-semibold">Recent Jobs ({jobs.length})</h3>
            </div>
            <div className="overflow-x-auto">
              <table className="min-w-full divide-y divide-gray-200">
                <thead className="bg-gray-50">
                  <tr>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Sample</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">User ID</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Status</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Step</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Created</th>
                  </tr>
                </thead>
                <tbody className="bg-white divide-y divide-gray-200">
                  {jobs.map((job) => (
                    <tr key={job.job_id} className="hover:bg-gray-50">
                      <td className="px-6 py-4 whitespace-nowrap text-sm font-medium text-gray-900">{job.sample_name}</td>
                      <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-600 font-mono text-xs">
                        User #{job.user_id}
                      </td>
                      <td className="px-6 py-4 whitespace-nowrap">
                        <span
                          className={`px-2 py-1 text-xs font-medium rounded-full ${
                            job.status === 'completed'
                              ? 'bg-green-100 text-green-800'
                              : job.status === 'running'
                              ? 'bg-blue-100 text-blue-800'
                              : job.status === 'failed'
                              ? 'bg-red-100 text-red-800'
                              : 'bg-gray-100 text-gray-800'
                          }`}
                        >
                          {job.status}
                        </span>
                      </td>
                      <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-600">{job.current_step || '-'}</td>
                      <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-600">
                        {new Date(job.created_at).toLocaleDateString()}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        )}

        {/* System Tab */}
        {activeTab === 'system' && (
          <div className="space-y-6">
            <div className="bg-white rounded-lg shadow p-6">
              <h3 className="text-lg font-semibold mb-4">System Resources</h3>
              <div className="space-y-4">
                <div>
                  <div className="flex justify-between text-sm mb-1">
                    <span className="text-gray-600">CPU</span>
                    <span className="font-medium">{stats.system.cpu.usage_percent.toFixed(1)}% ({stats.system.cpu.count} cores)</span>
                  </div>
                  <div className="w-full bg-gray-200 rounded-full h-2">
                    <div
                      className={`h-2 rounded-full ${
                        stats.system.cpu.usage_percent > 80 ? 'bg-red-500' : stats.system.cpu.usage_percent > 60 ? 'bg-yellow-500' : 'bg-green-500'
                      }`}
                      style={{ width: `${stats.system.cpu.usage_percent}%` }}
                    ></div>
                  </div>
                </div>

                <div>
                  <div className="flex justify-between text-sm mb-1">
                    <span className="text-gray-600">Memory</span>
                    <span className="font-medium">
                      {stats.system.memory.used_gb.toFixed(1)} GB / {stats.system.memory.total_gb.toFixed(1)} GB ({stats.system.memory.usage_percent.toFixed(1)}%)
                    </span>
                  </div>
                  <div className="w-full bg-gray-200 rounded-full h-2">
                    <div
                      className={`h-2 rounded-full ${
                        stats.system.memory.usage_percent > 80 ? 'bg-red-500' : stats.system.memory.usage_percent > 60 ? 'bg-yellow-500' : 'bg-green-500'
                      }`}
                      style={{ width: `${stats.system.memory.usage_percent}%` }}
                    ></div>
                  </div>
                </div>

                <div>
                  <div className="flex justify-between text-sm mb-1">
                    <span className="text-gray-600">Disk</span>
                    <span className="font-medium">
                      {stats.system.disk.used_gb.toFixed(1)} GB / {stats.system.disk.total_gb.toFixed(1)} GB ({stats.system.disk.usage_percent.toFixed(1)}%)
                    </span>
                  </div>
                  <div className="w-full bg-gray-200 rounded-full h-2">
                    <div
                      className={`h-2 rounded-full ${
                        stats.system.disk.usage_percent > 80 ? 'bg-red-500' : stats.system.disk.usage_percent > 60 ? 'bg-yellow-500' : 'bg-green-500'
                      }`}
                      style={{ width: `${stats.system.disk.usage_percent}%` }}
                    ></div>
                  </div>
                </div>
              </div>
            </div>

            <div className="bg-white rounded-lg shadow p-6">
              <h3 className="text-lg font-semibold mb-4">Pipeline Status</h3>
              <div className="grid grid-cols-2 gap-4">
                <div className="border border-gray-200 rounded-lg p-4">
                  <div className="text-3xl font-bold text-blue-600">{stats.nextflow_processes}</div>
                  <div className="text-sm text-gray-600 mt-1">Active Nextflow Processes</div>
                </div>
                <div className="border border-gray-200 rounded-lg p-4">
                  <div className="text-3xl font-bold text-green-600">{stats.jobs.running}</div>
                  <div className="text-sm text-gray-600 mt-1">Running Jobs</div>
                </div>
              </div>
            </div>
          </div>
        )}

        {/* Chat Tab */}
        {activeTab === 'chat' && (
          <AdminChatPanel />
        )}
      </div>
    </div>
  );
}

// Helper Components
function StatCard({ title, value, subtitle, icon, color }: any) {
  const colorClasses = {
    blue: 'bg-blue-50 text-blue-600',
    green: 'bg-green-50 text-green-600',
    purple: 'bg-purple-50 text-purple-600',
    orange: 'bg-orange-50 text-orange-600',
  };

  return (
    <div className="bg-white rounded-lg shadow p-6">
      <div className="flex items-center justify-between">
        <div>
          <p className="text-sm text-gray-600 mb-1">{title}</p>
          <p className="text-2xl font-bold text-gray-900">{value}</p>
          <p className="text-xs text-gray-500 mt-1">{subtitle}</p>
        </div>
        <div className={`text-3xl ${colorClasses[color as keyof typeof colorClasses]} p-3 rounded-full`}>{icon}</div>
      </div>
    </div>
  );
}

function SystemMetric({ label, value, suffix, max, color }: any) {
  const safeValue = value ?? 0;
  const percentage = (safeValue / max) * 100;
  const colorClasses = {
    green: 'bg-green-500',
    yellow: 'bg-yellow-500',
    red: 'bg-red-500',
  };

  return (
    <div>
      <div className="flex justify-between text-sm mb-2">
        <span className="text-gray-600">{label}</span>
        <span className="font-semibold text-gray-900">
          {safeValue.toFixed(1)}
          {suffix}
        </span>
      </div>
      <div className="w-full bg-gray-200 rounded-full h-2">
        <div className={`h-2 rounded-full ${colorClasses[color as keyof typeof colorClasses]}`} style={{ width: `${percentage}%` }}></div>
      </div>
    </div>
  );
}
