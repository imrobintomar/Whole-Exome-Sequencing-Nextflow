'use client';

import { useState, useEffect } from 'react';
import { adminApi, AdminDashboardStats, AdminUser, AdminJob } from '../lib/api';
import { auth } from '../lib/firebase';
import { onAuthStateChanged, signOut } from 'firebase/auth';
import AdminLogin from './AdminLogin';
import AdminChatPanel from './AdminChatPanel';
import { useRouter } from 'next/navigation';

export default function AdminDashboard() {
  const router = useRouter();
  const [stats, setStats] = useState<AdminDashboardStats | null>(null);
  const [users, setUsers] = useState<AdminUser[]>([]);
  const [jobs, setJobs] = useState<AdminJob[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [activeTab, setActiveTab] = useState<'overview' | 'users' | 'jobs' | 'system' | 'chat' | 'email' | 'analytics'>('overview');
  const [authReady, setAuthReady] = useState(false);
  const [isAuthenticated, setIsAuthenticated] = useState(false);

  // New features state
  const [autoRefresh, setAutoRefresh] = useState(true);
  const [refreshInterval, setRefreshInterval] = useState(30000);
  const [selectedJobs, setSelectedJobs] = useState<string[]>([]);
  const [filters, setFilters] = useState({
    status: '',
    search: '',
    dateFrom: '',
    dateTo: ''
  });
  const [userSearchQuery, setUserSearchQuery] = useState('');

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
    if (!authReady || !autoRefresh) return;

    loadDashboardData();
    const interval = setInterval(loadDashboardData, refreshInterval);
    return () => clearInterval(interval);
  }, [authReady, autoRefresh, refreshInterval]);

  // Reload when filters change
  useEffect(() => {
    if (!authReady) return;
    loadDashboardData();
  }, [filters]);

  const loadDashboardData = async () => {
    try {
      const jobParams: any = { limit: 20 };
      if (filters.status) jobParams.status = filters.status;
      if (filters.search) jobParams.search = filters.search;
      if (filters.dateFrom) jobParams.date_from = filters.dateFrom;
      if (filters.dateTo) jobParams.date_to = filters.dateTo;

      const [dashboardStats, usersData, jobsData] = await Promise.all([
        adminApi.getDashboardStats(),
        adminApi.getAllUsers(),
        adminApi.getAllJobs(jobParams),
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

  const handleChangePlan = async (userUid: string, newPlan: string) => {
    if (!newPlan) return;

    try {
      await adminApi.updateUserSubscription(userUid, newPlan);
      alert(`Successfully updated plan to ${newPlan}`);
      // Reload users to show updated data
      const usersData = await adminApi.getAllUsers();
      setUsers(usersData.users);
    } catch (err: any) {
      console.error('Error updating subscription:', err);
      alert(`Failed to update subscription plan: ${err.response?.data?.detail || err.message}`);
    }
  };

  const handleBanUser = async (userUid: string, email: string) => {
    const reason = prompt(`Enter reason for banning ${email}:`);
    if (!reason) return;

    try {
      await adminApi.banUser(userUid, reason);
      alert(`User ${email} has been banned`);
      const usersData = await adminApi.getAllUsers();
      setUsers(usersData.users);
    } catch (err: any) {
      console.error('Error banning user:', err);
      alert(`Failed to ban user: ${err.response?.data?.detail || err.message}`);
    }
  };

  const handleUnbanUser = async (userUid: string, email: string) => {
    if (!confirm(`Unban user ${email}?`)) return;

    try {
      await adminApi.unbanUser(userUid);
      alert(`User ${email} has been unbanned`);
      const usersData = await adminApi.getAllUsers();
      setUsers(usersData.users);
    } catch (err: any) {
      console.error('Error unbanning user:', err);
      alert(`Failed to unban user: ${err.response?.data?.detail || err.message}`);
    }
  };

  const handleSuspendUser = async (userUid: string, email: string) => {
    if (!confirm(`Suspend user ${email}? They will not be able to use the service until activated.`)) return;

    try {
      await adminApi.suspendUser(userUid);
      alert(`User ${email} has been suspended`);
      const usersData = await adminApi.getAllUsers();
      setUsers(usersData.users);
    } catch (err: any) {
      console.error('Error suspending user:', err);
      alert(`Failed to suspend user: ${err.response?.data?.detail || err.message}`);
    }
  };

  const handleActivateUser = async (userUid: string, email: string) => {
    if (!confirm(`Activate user ${email}?`)) return;

    try {
      await adminApi.activateUser(userUid);
      alert(`User ${email} has been activated`);
      const usersData = await adminApi.getAllUsers();
      setUsers(usersData.users);
    } catch (err: any) {
      console.error('Error activating user:', err);
      alert(`Failed to activate user: ${err.response?.data?.detail || err.message}`);
    }
  };

  const handleBulkAction = async (action: 'cancel' | 'delete') => {
    if (selectedJobs.length === 0) {
      alert('Please select jobs first');
      return;
    }

    if (!confirm(`${action.toUpperCase()} ${selectedJobs.length} selected job(s)?`)) {
      return;
    }

    try {
      const result = await adminApi.bulkJobAction(action, selectedJobs);
      alert(`Success: ${result.success.length}, Failed: ${result.failed.length}`);

      if (result.failed.length > 0) {
        console.error('Failed jobs:', result.failed);
      }

      setSelectedJobs([]);
      loadDashboardData();
    } catch (err: any) {
      console.error('Bulk action error:', err);
      alert('Failed to perform bulk action');
    }
  };

  const handleExportUsers = async () => {
    try {
      await adminApi.exportUsers();
    } catch (err: any) {
      console.error('Export error:', err);
      alert('Failed to export users');
    }
  };

  const toggleJobSelection = (jobId: string) => {
    setSelectedJobs(prev =>
      prev.includes(jobId)
        ? prev.filter(id => id !== jobId)
        : [...prev, jobId]
    );
  };

  const toggleSelectAll = () => {
    if (selectedJobs.length === jobs.length) {
      setSelectedJobs([]);
    } else {
      setSelectedJobs(jobs.map(j => j.job_id));
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
              <p className="text-gray-600 mt-1">ATGCFLOW</p>
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
              { id: 'analytics', label: 'Analytics' },
              { id: 'chat', label: 'Support Chat' },
              { id: 'email', label: 'Email Alerts' },
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
                trend="+8.3%"
                trendDirection="up"
              />
              <StatCard
                title="Total Jobs"
                value={stats.jobs.total}
                subtitle={`${stats.jobs.running} currently running`}
                icon="🔬"
                color="green"
                trend="+15.2%"
                trendDirection="up"
              />
              <StatCard
                title="Monthly Revenue"
                value={`$${stats.revenue.monthly_recurring_usd.toFixed(2)}`}
                subtitle="MRR"
                icon="💰"
                color="purple"
                trend="+12.5%"
                trendDirection="up"
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
              <div className="flex justify-between items-center mb-4">
                <h3 className="text-lg font-semibold">All Users ({users.filter(u =>
                  u.email.toLowerCase().includes(userSearchQuery.toLowerCase())
                ).length})</h3>
                <button
                  onClick={handleExportUsers}
                  className="px-4 py-2 bg-blue-500 text-white rounded hover:bg-blue-600 text-sm flex items-center gap-2"
                >
                  📥 Export to CSV
                </button>
              </div>
              {/* User Search */}
              <div className="relative">
                <input
                  type="text"
                  placeholder="🔍 Search users by email..."
                  value={userSearchQuery}
                  onChange={(e) => setUserSearchQuery(e.target.value)}
                  className="w-full px-4 py-2 pl-10 border border-gray-300 rounded-lg focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                />
                <div className="absolute left-3 top-2.5 text-gray-400">
                  🔍
                </div>
                {userSearchQuery && (
                  <button
                    onClick={() => setUserSearchQuery('')}
                    className="absolute right-3 top-2.5 text-gray-400 hover:text-gray-600"
                  >
                    ✕
                  </button>
                )}
              </div>
            </div>
            <div className="overflow-x-auto">
              <table className="min-w-full divide-y divide-gray-200">
                <thead className="bg-gray-50">
                  <tr>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Email</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Plan</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Subscription</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Account Status</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Usage</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Created</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Actions</th>
                  </tr>
                </thead>
                <tbody className="bg-white divide-y divide-gray-200">
                  {users.filter(u =>
                    u.email.toLowerCase().includes(userSearchQuery.toLowerCase())
                  ).map((user) => {
                    const accountStatus = user.is_banned ? 'banned' : !user.is_active ? 'suspended' : 'active';
                    return (
                      <tr key={user.uid} className="hover:bg-gray-50">
                        <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-900">
                          <div>
                            <button
                              onClick={() => router.push(`/admin/user/${user.uid}`)}
                              className="text-blue-600 hover:text-blue-800 hover:underline"
                            >
                              {user.email}
                            </button>
                            {user.is_banned && user.ban_reason && (
                              <div className="text-xs text-red-600 mt-1">Reason: {user.ban_reason}</div>
                            )}
                          </div>
                        </td>
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
                        <td className="px-6 py-4 whitespace-nowrap">
                          <span
                            className={`px-2 py-1 text-xs font-medium rounded-full ${
                              accountStatus === 'active'
                                ? 'bg-green-100 text-green-800'
                                : accountStatus === 'banned'
                                ? 'bg-red-100 text-red-800'
                                : 'bg-yellow-100 text-yellow-800'
                            }`}
                          >
                            {accountStatus}
                          </span>
                        </td>
                        <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-600">
                          {user.usage ? `${user.usage.jobs_executed}/${user.usage.jobs_limit}` : '0/2'}
                        </td>
                        <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-600">
                          {new Date(user.created_at).toLocaleDateString()}
                        </td>
                        <td className="px-6 py-4 whitespace-nowrap text-sm font-medium">
                          <div className="flex flex-col space-y-2">
                            {/* Usage Controls */}
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
                            {/* Plan Selector */}
                            <select
                              value={user.subscription?.plan || 'Free'}
                              onChange={(e) => handleChangePlan(user.uid, e.target.value)}
                              className="px-2 py-1 text-xs border border-gray-300 rounded focus:outline-none focus:ring-2 focus:ring-purple-500 focus:border-transparent"
                              title="Change subscription plan"
                            >
                              <option value="Free">Free (2 jobs/mo)</option>
                              <option value="Basic">Basic (10 jobs/mo, $29)</option>
                              <option value="Pro">Pro (50 jobs/mo, $99)</option>
                            </select>
                            {/* Ban/Suspend Controls */}
                            <div className="flex flex-col space-y-1">
                              {user.is_banned ? (
                                <button
                                  onClick={() => handleUnbanUser(user.uid, user.email)}
                                  className="px-2 py-1 bg-green-500 text-white rounded text-xs hover:bg-green-600"
                                >
                                  Unban
                                </button>
                              ) : (
                                <button
                                  onClick={() => handleBanUser(user.uid, user.email)}
                                  className="px-2 py-1 bg-red-500 text-white rounded text-xs hover:bg-red-600"
                                >
                                  Ban
                                </button>
                              )}
                              {!user.is_active ? (
                                <button
                                  onClick={() => handleActivateUser(user.uid, user.email)}
                                  className="px-2 py-1 bg-blue-500 text-white rounded text-xs hover:bg-blue-600"
                                >
                                  Activate
                                </button>
                              ) : (
                                <button
                                  onClick={() => handleSuspendUser(user.uid, user.email)}
                                  className="px-2 py-1 bg-yellow-500 text-white rounded text-xs hover:bg-yellow-600"
                                >
                                  Suspend
                                </button>
                              )}
                            </div>
                          </div>
                        </td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </div>
          </div>
        )}

        {/* Jobs Tab */}
        {activeTab === 'jobs' && (
          <div className="space-y-4">
            {/* Job Status Summary Cards */}
            <div className="grid grid-cols-1 md:grid-cols-4 gap-4">
              <div className="bg-white rounded-lg shadow p-6 border-l-4 border-blue-500">
                <div className="flex items-center justify-between">
                  <div>
                    <p className="text-sm text-gray-600 mb-1">Running</p>
                    <p className="text-3xl font-bold text-blue-600">{jobs.filter(j => j.status === 'running').length}</p>
                    <p className="text-xs text-gray-500 mt-1">
                      {jobs.length > 0 ? ((jobs.filter(j => j.status === 'running').length / jobs.length) * 100).toFixed(1) : 0}% of total
                    </p>
                  </div>
                  <div className="text-4xl">🔵</div>
                </div>
              </div>
              <div className="bg-white rounded-lg shadow p-6 border-l-4 border-yellow-500">
                <div className="flex items-center justify-between">
                  <div>
                    <p className="text-sm text-gray-600 mb-1">Pending</p>
                    <p className="text-3xl font-bold text-yellow-600">{jobs.filter(j => j.status === 'pending').length}</p>
                    <p className="text-xs text-gray-500 mt-1">
                      {jobs.length > 0 ? ((jobs.filter(j => j.status === 'pending').length / jobs.length) * 100).toFixed(1) : 0}% of total
                    </p>
                  </div>
                  <div className="text-4xl">🟡</div>
                </div>
              </div>
              <div className="bg-white rounded-lg shadow p-6 border-l-4 border-green-500">
                <div className="flex items-center justify-between">
                  <div>
                    <p className="text-sm text-gray-600 mb-1">Completed</p>
                    <p className="text-3xl font-bold text-green-600">{jobs.filter(j => j.status === 'completed').length}</p>
                    <p className="text-xs text-gray-500 mt-1">
                      {jobs.length > 0 ? ((jobs.filter(j => j.status === 'completed').length / jobs.length) * 100).toFixed(1) : 0}% of total
                    </p>
                  </div>
                  <div className="text-4xl">🟢</div>
                </div>
              </div>
              <div className="bg-white rounded-lg shadow p-6 border-l-4 border-red-500">
                <div className="flex items-center justify-between">
                  <div>
                    <p className="text-sm text-gray-600 mb-1">Failed</p>
                    <p className="text-3xl font-bold text-red-600">{jobs.filter(j => j.status === 'failed').length}</p>
                    <p className="text-xs text-gray-500 mt-1">
                      {jobs.length > 0 ? ((jobs.filter(j => j.status === 'failed').length / jobs.length) * 100).toFixed(1) : 0}% of total
                    </p>
                  </div>
                  <div className="text-4xl">🔴</div>
                </div>
              </div>
            </div>

            {/* Jobs Table */}
            <div className="bg-white rounded-lg shadow overflow-hidden">
              <div className="px-6 py-4 border-b border-gray-200">
                <div className="flex justify-between items-center mb-4">
                  <h3 className="text-lg font-semibold">Jobs Management ({jobs.length})</h3>
                <div className="flex gap-2">
                  <button
                    onClick={() => setAutoRefresh(!autoRefresh)}
                    className={`px-3 py-1 text-sm rounded ${autoRefresh ? 'bg-green-100 text-green-700' : 'bg-gray-100 text-gray-700'}`}
                  >
                    {autoRefresh ? '⏸️ Pause' : '▶️ Resume'} Auto-refresh
                  </button>
                  <select
                    value={refreshInterval}
                    onChange={(e) => setRefreshInterval(Number(e.target.value))}
                    className="px-3 py-1 text-sm border rounded"
                  >
                    <option value={10000}>10s</option>
                    <option value={30000}>30s</option>
                    <option value={60000}>1m</option>
                  </select>
                </div>
              </div>

              {/* Filters */}
              <div className="grid grid-cols-1 md:grid-cols-4 gap-3 mb-4">
                <input
                  type="text"
                  placeholder="Search by job ID or sample name..."
                  value={filters.search}
                  onChange={(e) => setFilters({...filters, search: e.target.value})}
                  className="px-3 py-2 border rounded text-sm"
                />
                <select
                  value={filters.status}
                  onChange={(e) => setFilters({...filters, status: e.target.value})}
                  className="px-3 py-2 border rounded text-sm"
                >
                  <option value="">All Status</option>
                  <option value="pending">Pending</option>
                  <option value="running">Running</option>
                  <option value="completed">Completed</option>
                  <option value="failed">Failed</option>
                </select>
                <input
                  type="date"
                  value={filters.dateFrom}
                  onChange={(e) => setFilters({...filters, dateFrom: e.target.value})}
                  className="px-3 py-2 border rounded text-sm"
                  placeholder="From"
                />
                <input
                  type="date"
                  value={filters.dateTo}
                  onChange={(e) => setFilters({...filters, dateTo: e.target.value})}
                  className="px-3 py-2 border rounded text-sm"
                  placeholder="To"
                />
              </div>

              {/* Bulk Actions */}
              {selectedJobs.length > 0 && (
                <div className="bg-blue-50 border border-blue-200 rounded p-3 mb-4 flex items-center gap-3">
                  <span className="text-sm font-medium">{selectedJobs.length} selected</span>
                  <button
                    onClick={() => handleBulkAction('cancel')}
                    className="px-3 py-1 bg-yellow-500 text-white rounded text-sm hover:bg-yellow-600"
                  >
                    Cancel Selected
                  </button>
                  <button
                    onClick={() => handleBulkAction('delete')}
                    className="px-3 py-1 bg-red-500 text-white rounded text-sm hover:bg-red-600"
                  >
                    Delete Selected
                  </button>
                  <button
                    onClick={() => setSelectedJobs([])}
                    className="px-3 py-1 bg-gray-500 text-white rounded text-sm hover:bg-gray-600"
                  >
                    Clear Selection
                  </button>
                </div>
              )}
            </div>

            <div className="overflow-x-auto">
              <table className="min-w-full divide-y divide-gray-200">
                <thead className="bg-gray-50">
                  <tr>
                    <th className="px-6 py-3 text-left">
                      <input
                        type="checkbox"
                        checked={selectedJobs.length === jobs.length && jobs.length > 0}
                        onChange={toggleSelectAll}
                        className="rounded"
                      />
                    </th>
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
                      <td className="px-6 py-4 whitespace-nowrap">
                        <input
                          type="checkbox"
                          checked={selectedJobs.includes(job.job_id)}
                          onChange={() => toggleJobSelection(job.job_id)}
                          className="rounded"
                        />
                      </td>
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
                      {stats.system.memory.used_gb.toFixed(1)} GB / {stats.system.memory.total_gb.toFixed(1)} GB ({stats.system.memory.percent.toFixed(1)}%)
                    </span>
                  </div>
                  <div className="w-full bg-gray-200 rounded-full h-2">
                    <div
                      className={`h-2 rounded-full ${
                        stats.system.memory.percent > 80 ? 'bg-red-500' : stats.system.memory.percent > 60 ? 'bg-yellow-500' : 'bg-green-500'
                      }`}
                      style={{ width: `${stats.system.memory.percent}%` }}
                    ></div>
                  </div>
                </div>

                <div>
                  <div className="flex justify-between text-sm mb-1">
                    <span className="text-gray-600">Disk</span>
                    <span className="font-medium">
                      {stats.system.disk.used_gb.toFixed(1)} GB / {stats.system.disk.total_gb.toFixed(1)} GB ({stats.system.disk.percent.toFixed(1)}%)
                    </span>
                  </div>
                  <div className="w-full bg-gray-200 rounded-full h-2">
                    <div
                      className={`h-2 rounded-full ${
                        stats.system.disk.percent > 80 ? 'bg-red-500' : stats.system.disk.percent > 60 ? 'bg-yellow-500' : 'bg-green-500'
                      }`}
                      style={{ width: `${stats.system.disk.percent}%` }}
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
                  <div className="text-sm text-gray-600 mt-1">Active Processes</div>
                </div>
                <div className="border border-gray-200 rounded-lg p-4">
                  <div className="text-3xl font-bold text-green-600">{stats.jobs.running}</div>
                  <div className="text-sm text-gray-600 mt-1">Running Jobs</div>
                </div>
              </div>
            </div>
          </div>
        )}

        {/* Analytics Tab */}
        {activeTab === 'analytics' && (
          <AnalyticsPanel />
        )}

        {/* Chat Tab */}
        {activeTab === 'chat' && (
          <AdminChatPanel />
        )}

        {/* Email Alerts Tab */}
        {activeTab === 'email' && (
          <EmailAlertsPanel />
        )}
      </div>
    </div>
  );
}

// Helper Components
function StatCard({ title, value, subtitle, icon, color, trend, trendDirection }: any) {
  const colorClasses = {
    blue: 'bg-blue-50 text-blue-600',
    green: 'bg-green-50 text-green-600',
    purple: 'bg-purple-50 text-purple-600',
    orange: 'bg-orange-50 text-orange-600',
  };

  return (
    <div className="bg-white rounded-lg shadow p-6">
      <div className="flex items-center justify-between">
        <div className="flex-1">
          <p className="text-sm text-gray-600 mb-1">{title}</p>
          <div className="flex items-baseline gap-2">
            <p className="text-2xl font-bold text-gray-900">{value}</p>
            {trend && (
              <span className={`text-sm font-semibold ${
                trendDirection === 'up' ? 'text-green-600' : trendDirection === 'down' ? 'text-red-600' : 'text-gray-600'
              }`}>
                {trendDirection === 'up' ? '↑' : trendDirection === 'down' ? '↓' : '→'} {trend}
              </span>
            )}
          </div>
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

function AnalyticsPanel() {
  const [loading, setLoading] = useState(true);
  const [period, setPeriod] = useState<'day' | 'week' | 'month' | 'year'>('month');
  const [summary, setSummary] = useState<any>(null);
  const [revenueData, setRevenueData] = useState<any>(null);
  const [userData, setUserData] = useState<any>(null);
  const [jobData, setJobData] = useState<any>(null);

  useEffect(() => {
    loadAnalyticsData();
  }, [period]);

  const loadAnalyticsData = async () => {
    try {
      setLoading(true);
      const [summaryData, revenue, users, jobs] = await Promise.all([
        adminApi.getAnalyticsSummary(period),
        adminApi.getRevenueAnalytics(),
        adminApi.getUserAnalytics(),
        adminApi.getJobAnalytics()
      ]);

      setSummary(summaryData);
      setRevenueData(revenue);
      setUserData(users);
      setJobData(jobs);
    } catch (error) {
      console.error('Failed to load analytics:', error);
    } finally {
      setLoading(false);
    }
  };

  if (loading) {
    return (
      <div className="flex items-center justify-center h-64">
        <div className="text-center">
          <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-600 mx-auto mb-4"></div>
          <p className="text-gray-600">Loading analytics...</p>
        </div>
      </div>
    );
  }

  if (!summary || !revenueData || !userData || !jobData) {
    return (
      <div className="bg-yellow-50 border border-yellow-200 rounded-lg p-6">
        <p className="text-yellow-800">No analytics data available</p>
      </div>
    );
  }

  // Import Recharts components dynamically
  const {
    LineChart, Line, BarChart, Bar, PieChart, Pie, Cell, AreaChart, Area,
    XAxis, YAxis, Tooltip, Legend, ResponsiveContainer
  } = require('recharts');

  const COLORS = ['#8b5cf6', '#10b981', '#3b82f6', '#f59e0b', '#ef4444', '#6b7280'];

  return (
    <div className="space-y-6">
      {/* Header with Period Selector */}
      <div className="flex justify-between items-center">
        <div>
          <h2 className="text-2xl font-bold text-gray-900">Analytics Dashboard</h2>
          <p className="text-gray-600 mt-1">Insights and trends for your platform</p>
        </div>
        <select
          value={period}
          onChange={(e) => setPeriod(e.target.value as any)}
          className="px-4 py-2 border border-gray-300 rounded-lg focus:ring-2 focus:ring-blue-500 focus:border-transparent"
        >
          <option value="day">Last 24 Hours</option>
          <option value="week">Last 7 Days</option>
          <option value="month">Last 30 Days</option>
          <option value="year">Last 12 Months</option>
        </select>
      </div>

      {/* Key Metrics Cards */}
      <div className="grid grid-cols-1 md:grid-cols-4 gap-4">
        <div className="bg-white rounded-lg shadow p-6">
          <div className="flex items-center justify-between">
            <div>
              <p className="text-sm text-gray-600 mb-1">MRR</p>
              <p className="text-2xl font-bold text-gray-900">${revenueData.current_mrr.toFixed(2)}</p>
              <p className="text-xs text-gray-500 mt-1">ARPU: ${revenueData.arpu.toFixed(2)}</p>
            </div>
            <div className="text-3xl">💰</div>
          </div>
        </div>

        <div className="bg-white rounded-lg shadow p-6">
          <div className="flex items-center justify-between">
            <div>
              <p className="text-sm text-gray-600 mb-1">New Users</p>
              <p className="text-2xl font-bold text-gray-900">{summary.users.new_this_period}</p>
              <p className="text-xs text-gray-500 mt-1">Total: {summary.users.total}</p>
            </div>
            <div className="text-3xl">👥</div>
          </div>
        </div>

        <div className="bg-white rounded-lg shadow p-6">
          <div className="flex items-center justify-between">
            <div>
              <p className="text-sm text-gray-600 mb-1">Success Rate</p>
              <p className="text-2xl font-bold text-gray-900">{summary.jobs.success_rate.toFixed(1)}%</p>
              <p className="text-xs text-gray-500 mt-1">{summary.jobs.completed}/{summary.jobs.this_period} jobs</p>
            </div>
            <div className="text-3xl">✅</div>
          </div>
        </div>

        <div className="bg-white rounded-lg shadow p-6">
          <div className="flex items-center justify-between">
            <div>
              <p className="text-sm text-gray-600 mb-1">Avg Duration</p>
              <p className="text-2xl font-bold text-gray-900">{jobData.avg_duration_hours.toFixed(1)}h</p>
              <p className="text-xs text-gray-500 mt-1">{jobData.total_jobs} jobs analyzed</p>
            </div>
            <div className="text-3xl">⏱️</div>
          </div>
        </div>
      </div>

      {/* Revenue Over Time */}
      <div className="bg-white rounded-lg shadow p-6">
        <h3 className="text-lg font-semibold mb-4">Revenue Over Time</h3>
        <ResponsiveContainer width="100%" height={300}>
          <LineChart data={revenueData.revenue_over_time}>
            <XAxis dataKey="month" />
            <YAxis />
            <Tooltip formatter={(value: any) => `$${value.toFixed(2)}`} />
            <Legend />
            <Line
              type="monotone"
              dataKey="revenue"
              stroke="#8b5cf6"
              strokeWidth={2}
              dot={{ r: 4 }}
              name="Revenue"
            />
          </LineChart>
        </ResponsiveContainer>
      </div>

      {/* Revenue by Plan & User Growth */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
        {/* Revenue by Plan */}
        <div className="bg-white rounded-lg shadow p-6">
          <h3 className="text-lg font-semibold mb-4">Revenue by Plan</h3>
          <ResponsiveContainer width="100%" height={300}>
            <PieChart>
              <Pie
                data={revenueData.revenue_by_plan}
                dataKey="revenue"
                nameKey="plan"
                cx="50%"
                cy="50%"
                outerRadius={100}
                label={({ name, value }: { name: string; value: number }) => `${name}: $${value.toFixed(0)}`}
              >
                {revenueData.revenue_by_plan.map((entry: any, index: number) => (
                  <Cell key={`cell-${index}`} fill={COLORS[index % COLORS.length]} />
                ))}
              </Pie>
              <Tooltip formatter={(value: any) => `$${value.toFixed(2)}`} />
            </PieChart>
          </ResponsiveContainer>
        </div>

        {/* Users by Plan */}
        <div className="bg-white rounded-lg shadow p-6">
          <h3 className="text-lg font-semibold mb-4">Users by Plan</h3>
          <ResponsiveContainer width="100%" height={300}>
            <PieChart>
              <Pie
                data={userData.users_by_plan}
                dataKey="count"
                nameKey="plan"
                cx="50%"
                cy="50%"
                outerRadius={100}
                label={({ name, value }: { name: string; value: number }) => `${name}: ${value}`}
              >
                {userData.users_by_plan.map((entry: any, index: number) => (
                  <Cell key={`cell-${index}`} fill={COLORS[index % COLORS.length]} />
                ))}
              </Pie>
              <Tooltip />
            </PieChart>
          </ResponsiveContainer>
        </div>
      </div>

      {/* User Growth Over Time */}
      <div className="bg-white rounded-lg shadow p-6">
        <h3 className="text-lg font-semibold mb-4">User Growth</h3>
        <ResponsiveContainer width="100%" height={300}>
          <AreaChart data={userData.user_growth}>
            <XAxis dataKey="month" />
            <YAxis />
            <Tooltip />
            <Legend />
            <Area
              type="monotone"
              dataKey="total_users"
              stroke="#10b981"
              fill="#10b981"
              fillOpacity={0.6}
              name="Total Users"
            />
            <Area
              type="monotone"
              dataKey="new_users"
              stroke="#3b82f6"
              fill="#3b82f6"
              fillOpacity={0.3}
              name="New Users"
            />
          </AreaChart>
        </ResponsiveContainer>
      </div>

      {/* Jobs Over Time */}
      <div className="bg-white rounded-lg shadow p-6">
        <h3 className="text-lg font-semibold mb-4">Jobs Over Time</h3>
        <ResponsiveContainer width="100%" height={300}>
          <BarChart data={jobData.jobs_over_time}>
            <XAxis dataKey="date" />
            <YAxis />
            <Tooltip />
            <Legend />
            <Bar dataKey="completed" stackId="a" fill="#10b981" name="Completed" />
            <Bar dataKey="failed" stackId="a" fill="#ef4444" name="Failed" />
            <Bar dataKey="running" stackId="a" fill="#3b82f6" name="Running" />
            <Bar dataKey="pending" stackId="a" fill="#f59e0b" name="Pending" />
          </BarChart>
        </ResponsiveContainer>
      </div>

      {/* Jobs by Status */}
      <div className="bg-white rounded-lg shadow p-6">
        <h3 className="text-lg font-semibold mb-4">Jobs by Status (Current Period)</h3>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
          <ResponsiveContainer width="100%" height={300}>
            <PieChart>
              <Pie
                data={jobData.jobs_by_status}
                dataKey="count"
                nameKey="status"
                cx="50%"
                cy="50%"
                outerRadius={100}
                label={({ name, value }: { name: string; value: number }) => `${name}: ${value}`}
              >
                {jobData.jobs_by_status.map((entry: any, index: number) => {
                  const statusColors: any = {
                    completed: '#10b981',
                    failed: '#ef4444',
                    running: '#3b82f6',
                    pending: '#f59e0b'
                  };
                  return (
                    <Cell key={`cell-${index}`} fill={statusColors[entry.status] || COLORS[index % COLORS.length]} />
                  );
                })}
              </Pie>
              <Tooltip />
            </PieChart>
          </ResponsiveContainer>

          <div className="flex flex-col justify-center space-y-4">
            {jobData.jobs_by_status.map((item: any) => {
              const statusColors: any = {
                completed: 'bg-green-100 text-green-800 border-green-200',
                failed: 'bg-red-100 text-red-800 border-red-200',
                running: 'bg-blue-100 text-blue-800 border-blue-200',
                pending: 'bg-yellow-100 text-yellow-800 border-yellow-200'
              };
              return (
                <div
                  key={item.status}
                  className={`p-4 rounded-lg border-2 ${statusColors[item.status] || 'bg-gray-100 text-gray-800 border-gray-200'}`}
                >
                  <div className="flex justify-between items-center">
                    <span className="font-semibold capitalize">{item.status}</span>
                    <span className="text-2xl font-bold">{item.count}</span>
                  </div>
                  <div className="text-sm mt-1">
                    {((item.count / jobData.total_jobs) * 100).toFixed(1)}% of total
                  </div>
                </div>
              );
            })}
          </div>
        </div>
      </div>

      {/* Summary Stats */}
      <div className="bg-gradient-to-r from-blue-50 to-purple-50 rounded-lg p-6 border border-blue-100">
        <h3 className="text-lg font-semibold mb-4">Summary Statistics</h3>
        <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
          <div>
            <p className="text-sm text-gray-600">Total Revenue (Period)</p>
            <p className="text-xl font-bold text-gray-900">${revenueData.total_revenue.toFixed(2)}</p>
          </div>
          <div>
            <p className="text-sm text-gray-600">Total Users</p>
            <p className="text-xl font-bold text-gray-900">{userData.total_users}</p>
          </div>
          <div>
            <p className="text-sm text-gray-600">Success Rate</p>
            <p className="text-xl font-bold text-gray-900">{jobData.success_rate.toFixed(1)}%</p>
          </div>
        </div>
      </div>
    </div>
  );
}

function EmailAlertsPanel() {
  const [emailSettings, setEmailSettings] = useState<any>(null);
  const [healthAlerts, setHealthAlerts] = useState<any[]>([]);
  const [loading, setLoading] = useState(true);
  const [testingConnection, setTestingConnection] = useState(false);
  const [sendingEmail, setSendingEmail] = useState(false);
  const [message, setMessage] = useState<{ type: 'success' | 'error'; text: string } | null>(null);

  // Custom email form state
  const [customEmail, setCustomEmail] = useState({
    userEmail: '',
    subject: '',
    message: '',
    userName: ''
  });

  // Payment reminder form state
  const [paymentReminder, setPaymentReminder] = useState({
    userUid: '',
    amount: '',
    dueDate: '',
    invoiceUrl: ''
  });

  useEffect(() => {
    loadEmailData();
  }, []);

  const loadEmailData = async () => {
    try {
      setLoading(true);
      const [settings, alerts] = await Promise.all([
        adminApi.getEmailSettings(),
        adminApi.getHealthAlerts(50)
      ]);
      setEmailSettings(settings);
      setHealthAlerts(alerts.alerts);
    } catch (error) {
      console.error('Failed to load email data:', error);
      setMessage({ type: 'error', text: 'Failed to load email data' });
    } finally {
      setLoading(false);
    }
  };

  const testEmailConnection = async () => {
    try {
      setTestingConnection(true);
      setMessage(null);
      const result = await adminApi.testEmailConnection();
      setMessage({
        type: result.status === 'success' ? 'success' : 'error',
        text: result.message
      });
    } catch (error: any) {
      setMessage({ type: 'error', text: error.response?.data?.detail || 'Failed to test connection' });
    } finally {
      setTestingConnection(false);
    }
  };

  const sendTestHealthAlert = async () => {
    try {
      setSendingEmail(true);
      setMessage(null);
      await adminApi.sendTestHealthAlert();
      setMessage({ type: 'success', text: 'Test health alert email sent successfully!' });
    } catch (error: any) {
      setMessage({ type: 'error', text: error.response?.data?.detail || 'Failed to send test email' });
    } finally {
      setSendingEmail(false);
    }
  };

  const sendCustomEmail = async (e: React.FormEvent) => {
    e.preventDefault();
    try {
      setSendingEmail(true);
      setMessage(null);
      const result = await adminApi.sendCustomEmail(
        customEmail.userEmail,
        customEmail.subject,
        customEmail.message,
        customEmail.userName || undefined
      );
      setMessage({
        type: 'success',
        text: ` Email sent successfully to ${customEmail.userEmail}! Check their inbox (and spam folder).`
      });
      setCustomEmail({ userEmail: '', subject: '', message: '', userName: '' });
    } catch (error: any) {
      console.error('Email send error:', error);
      const errorMsg = error.response?.data?.detail || error.message || 'Failed to send email';
      setMessage({ type: 'error', text: ` ${errorMsg}` });
    } finally {
      setSendingEmail(false);
    }
  };

  const sendPaymentReminder = async (e: React.FormEvent) => {
    e.preventDefault();
    try {
      setSendingEmail(true);
      setMessage(null);
      const result = await adminApi.sendPaymentReminder(
        paymentReminder.userUid,
        parseFloat(paymentReminder.amount),
        paymentReminder.dueDate,
        paymentReminder.invoiceUrl || undefined
      );
      setMessage({
        type: 'success',
        text: ` Payment reminder sent successfully! Amount: $${paymentReminder.amount}, Due: ${paymentReminder.dueDate}`
      });
      setPaymentReminder({ userUid: '', amount: '', dueDate: '', invoiceUrl: '' });
    } catch (error: any) {
      console.error('Payment reminder error:', error);
      const errorMsg = error.response?.data?.detail || error.message || 'Failed to send payment reminder';
      setMessage({ type: 'error', text: ` ${errorMsg}` });
    } finally {
      setSendingEmail(false);
    }
  };

  const resolveAlert = async (alertId: number) => {
    try {
      await adminApi.resolveHealthAlert(alertId);
      setMessage({ type: 'success', text: 'Alert resolved successfully!' });
      loadEmailData();
    } catch (error: any) {
      setMessage({ type: 'error', text: error.response?.data?.detail || 'Failed to resolve alert' });
    }
  };

  if (loading) {
    return (
      <div className="flex items-center justify-center h-64">
        <div className="text-gray-500">Loading email settings...</div>
      </div>
    );
  }

  return (
    <div className="space-y-6">
      {/* Message Banner */}
      {message && (
        <div className={`p-4 rounded-lg ${message.type === 'success' ? 'bg-green-50 text-green-800' : 'bg-red-50 text-red-800'}`}>
          {message.text}
        </div>
      )}

      {/* Email Settings */}
      <div className="bg-white rounded-lg shadow p-6">
        <h2 className="text-xl font-bold mb-4">Email Configuration</h2>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <div>
            <p className="text-sm text-gray-600">SMTP Host</p>
            <p className="font-semibold">{emailSettings?.smtp_host}</p>
          </div>
          <div>
            <p className="text-sm text-gray-600">SMTP Port</p>
            <p className="font-semibold">{emailSettings?.smtp_port}</p>
          </div>
          <div>
            <p className="text-sm text-gray-600">SMTP User</p>
            <p className="font-semibold">{emailSettings?.smtp_user}</p>
          </div>
          <div>
            <p className="text-sm text-gray-600">Admin Email</p>
            <p className="font-semibold">{emailSettings?.admin_email}</p>
          </div>
          <div>
            <p className="text-sm text-gray-600">Health Alerts</p>
            <p className={`font-semibold ${emailSettings?.health_alerts_enabled ? 'text-green-600' : 'text-red-600'}`}>
              {emailSettings?.health_alerts_enabled ? 'Enabled' : 'Disabled'}
            </p>
          </div>
          <div>
            <p className="text-sm text-gray-600">Health Thresholds</p>
            <p className="font-semibold text-xs">
              CPU: {emailSettings?.health_thresholds?.cpu}% |
              Memory: {emailSettings?.health_thresholds?.memory}% |
              Disk: {emailSettings?.health_thresholds?.disk}%
            </p>
          </div>
        </div>
        <div className="mt-4 flex gap-2">
          <button
            onClick={testEmailConnection}
            disabled={testingConnection}
            className="px-4 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700 disabled:bg-gray-400"
          >
            {testingConnection ? 'Testing...' : 'Test Connection'}
          </button>
          <button
            onClick={sendTestHealthAlert}
            disabled={sendingEmail}
            className="px-4 py-2 bg-orange-600 text-white rounded-lg hover:bg-orange-700 disabled:bg-gray-400"
          >
            {sendingEmail ? 'Sending...' : 'Send Test Alert'}
          </button>
        </div>
      </div>

      {/* Send Custom Email */}
      <div className="bg-white rounded-lg shadow p-6">
        <h2 className="text-xl font-bold mb-4">Send Custom Email</h2>
        <form onSubmit={sendCustomEmail} className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div>
              <label className="block text-sm font-medium text-gray-700 mb-1">
                Recipient Email *
              </label>
              <input
                type="email"
                value={customEmail.userEmail}
                onChange={(e) => setCustomEmail({ ...customEmail, userEmail: e.target.value })}
                className="w-full px-3 py-2 border rounded-lg focus:ring-2 focus:ring-blue-500"
                required
              />
            </div>
            <div>
              <label className="block text-sm font-medium text-gray-700 mb-1">
                Recipient Name (Optional)
              </label>
              <input
                type="text"
                value={customEmail.userName}
                onChange={(e) => setCustomEmail({ ...customEmail, userName: e.target.value })}
                className="w-full px-3 py-2 border rounded-lg focus:ring-2 focus:ring-blue-500"
              />
            </div>
          </div>
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Subject *
            </label>
            <input
              type="text"
              value={customEmail.subject}
              onChange={(e) => setCustomEmail({ ...customEmail, subject: e.target.value })}
              className="w-full px-3 py-2 border rounded-lg focus:ring-2 focus:ring-blue-500"
              required
            />
          </div>
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Message *
            </label>
            <textarea
              value={customEmail.message}
              onChange={(e) => setCustomEmail({ ...customEmail, message: e.target.value })}
              className="w-full px-3 py-2 border rounded-lg focus:ring-2 focus:ring-blue-500"
              rows={4}
              required
            />
          </div>
          <button
            type="submit"
            disabled={sendingEmail}
            className="px-6 py-2 bg-green-600 text-white rounded-lg hover:bg-green-700 disabled:bg-gray-400"
          >
            {sendingEmail ? 'Sending...' : 'Send Email'}
          </button>
        </form>
      </div>

      {/* Send Payment Reminder */}
      <div className="bg-white rounded-lg shadow p-6">
        <h2 className="text-xl font-bold mb-4">Send Payment Reminder</h2>
        <form onSubmit={sendPaymentReminder} className="space-y-4">
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div>
              <label className="block text-sm font-medium text-gray-700 mb-1">
                User UID *
              </label>
              <input
                type="text"
                value={paymentReminder.userUid}
                onChange={(e) => setPaymentReminder({ ...paymentReminder, userUid: e.target.value })}
                className="w-full px-3 py-2 border rounded-lg focus:ring-2 focus:ring-blue-500"
                required
                placeholder="Firebase User UID"
              />
            </div>
            <div>
              <label className="block text-sm font-medium text-gray-700 mb-1">
                Amount (USD) *
              </label>
              <input
                type="number"
                step="0.01"
                value={paymentReminder.amount}
                onChange={(e) => setPaymentReminder({ ...paymentReminder, amount: e.target.value })}
                className="w-full px-3 py-2 border rounded-lg focus:ring-2 focus:ring-blue-500"
                required
                placeholder="29.99"
              />
            </div>
            <div>
              <label className="block text-sm font-medium text-gray-700 mb-1">
                Due Date *
              </label>
              <input
                type="date"
                value={paymentReminder.dueDate}
                onChange={(e) => setPaymentReminder({ ...paymentReminder, dueDate: e.target.value })}
                className="w-full px-3 py-2 border rounded-lg focus:ring-2 focus:ring-blue-500"
                required
              />
            </div>
            <div>
              <label className="block text-sm font-medium text-gray-700 mb-1">
                Invoice URL (Optional)
              </label>
              <input
                type="url"
                value={paymentReminder.invoiceUrl}
                onChange={(e) => setPaymentReminder({ ...paymentReminder, invoiceUrl: e.target.value })}
                className="w-full px-3 py-2 border rounded-lg focus:ring-2 focus:ring-blue-500"
                placeholder="https://..."
              />
            </div>
          </div>
          <button
            type="submit"
            disabled={sendingEmail}
            className="px-6 py-2 bg-purple-600 text-white rounded-lg hover:bg-purple-700 disabled:bg-gray-400"
          >
            {sendingEmail ? 'Sending...' : 'Send Payment Reminder'}
          </button>
        </form>
      </div>

      {/* Health Alerts History */}
      <div className="bg-white rounded-lg shadow p-6">
        <div className="flex justify-between items-center mb-4">
          <h2 className="text-xl font-bold">Health Alerts History</h2>
          <button
            onClick={loadEmailData}
            className="text-sm text-blue-600 hover:text-blue-700"
          >
            Refresh
          </button>
        </div>
        {healthAlerts.length === 0 ? (
          <p className="text-gray-500 text-center py-8">No health alerts found</p>
        ) : (
          <div className="space-y-2">
            {healthAlerts.map((alert) => (
              <div
                key={alert.id}
                className={`p-4 rounded-lg border-l-4 ${
                  alert.severity === 'critical'
                    ? 'border-red-500 bg-red-50'
                    : 'border-yellow-500 bg-yellow-50'
                }`}
              >
                <div className="flex justify-between items-start">
                  <div className="flex-1">
                    <div className="flex items-center gap-2 mb-1">
                      <span className={`text-xs font-semibold px-2 py-1 rounded ${
                        alert.severity === 'critical' ? 'bg-red-200 text-red-800' : 'bg-yellow-200 text-yellow-800'
                      }`}>
                        {alert.severity.toUpperCase()}
                      </span>
                      <span className="text-xs text-gray-500">{alert.type.toUpperCase()}</span>
                      {alert.resolved && (
                        <span className="text-xs bg-green-200 text-green-800 px-2 py-1 rounded">
                          RESOLVED
                        </span>
                      )}
                    </div>
                    <p className="text-sm font-medium text-gray-900">{alert.message}</p>
                    <p className="text-xs text-gray-600 mt-1">
                      Value: {alert.value.toFixed(1)}% | Threshold: {alert.threshold.toFixed(1)}%
                    </p>
                    <p className="text-xs text-gray-500 mt-1">
                      {new Date(alert.created_at).toLocaleString()}
                    </p>
                  </div>
                  {!alert.resolved && (
                    <button
                      onClick={() => resolveAlert(alert.id)}
                      className="ml-4 px-3 py-1 text-sm bg-gray-600 text-white rounded hover:bg-gray-700"
                    >
                      Resolve
                    </button>
                  )}
                </div>
              </div>
            ))}
          </div>
        )}
      </div>
    </div>
  );
}
