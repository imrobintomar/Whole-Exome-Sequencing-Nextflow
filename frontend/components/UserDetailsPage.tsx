'use client';

import { useState, useEffect } from 'react';
import { adminApi, UserDetailsResponse, UserNote, UserTag } from '../lib/api';
import { useRouter } from 'next/navigation';

interface UserDetailsPageProps {
  userUid: string;
}

export default function UserDetailsPage({ userUid }: UserDetailsPageProps) {
  const router = useRouter();
  const [data, setData] = useState<UserDetailsResponse | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [activeTab, setActiveTab] = useState<'profile' | 'jobs' | 'payments' | 'tickets' | 'activity' | 'notes'>('profile');

  // Notes/Tags state
  const [newNote, setNewNote] = useState('');
  const [newTagName, setNewTagName] = useState('');
  const [newTagColor, setNewTagColor] = useState('blue');
  const [addingNote, setAddingNote] = useState(false);
  const [addingTag, setAddingTag] = useState(false);

  useEffect(() => {
    loadUserDetails();
  }, [userUid]);

  const loadUserDetails = async () => {
    try {
      setLoading(true);
      const response = await adminApi.getUserDetails(userUid);
      setData(response);
      setError(null);
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to load user details');
      console.error('Error loading user details:', err);
    } finally {
      setLoading(false);
    }
  };

  const handleAddNote = async () => {
    if (!newNote.trim()) return;

    try {
      setAddingNote(true);
      await adminApi.addUserNote(userUid, newNote);
      setNewNote('');
      await loadUserDetails();
    } catch (err: any) {
      alert('Failed to add note: ' + (err.response?.data?.detail || err.message));
    } finally {
      setAddingNote(false);
    }
  };

  const handleAddTag = async () => {
    if (!newTagName.trim()) return;

    try {
      setAddingTag(true);
      await adminApi.addUserTag(userUid, newTagName, newTagColor);
      setNewTagName('');
      await loadUserDetails();
    } catch (err: any) {
      alert('Failed to add tag: ' + (err.response?.data?.detail || err.message));
    } finally {
      setAddingTag(false);
    }
  };

  const handleRemoveTag = async (tagId: number) => {
    if (!confirm('Remove this tag?')) return;

    try {
      await adminApi.removeUserTag(userUid, tagId);
      await loadUserDetails();
    } catch (err: any) {
      alert('Failed to remove tag: ' + (err.response?.data?.detail || err.message));
    }
  };

  if (loading) {
    return (
      <div className="flex items-center justify-center min-h-screen">
        <div className="text-center">
          <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-600 mx-auto mb-4"></div>
          <p className="text-gray-600">Loading user details...</p>
        </div>
      </div>
    );
  }

  if (error || !data) {
    return (
      <div className="flex items-center justify-center min-h-screen">
        <div className="bg-red-50 border border-red-200 rounded-lg p-6 max-w-md">
          <h3 className="text-red-800 font-semibold mb-2">Error</h3>
          <p className="text-red-600">{error || 'User not found'}</p>
          <button
            onClick={() => router.back()}
            className="mt-4 px-4 py-2 bg-gray-600 text-white rounded hover:bg-gray-700"
          >
            Go Back
          </button>
        </div>
      </div>
    );
  }

  const { user, subscription, usage, jobs, payment_history, support_tickets, activity_timeline, notes, tags } = data;

  const tagColors: Record<string, string> = {
    blue: 'bg-blue-100 text-blue-800',
    green: 'bg-green-100 text-green-800',
    red: 'bg-red-100 text-red-800',
    yellow: 'bg-yellow-100 text-yellow-800',
    purple: 'bg-purple-100 text-purple-800',
  };

  return (
    <div className="min-h-screen bg-gray-50">
      {/* Header */}
      <div className="bg-white border-b">
        <div className="max-w-7xl mx-auto px-4 py-6">
          <div className="flex justify-between items-start">
            <div className="flex-1">
              <div className="flex items-center gap-4 mb-2">
                <button
                  onClick={() => router.back()}
                  className="px-3 py-1 text-sm bg-gray-100 hover:bg-gray-200 rounded"
                >
                  ← Back
                </button>
                <h1 className="text-3xl font-bold text-gray-900">{user.email}</h1>
              </div>

              {/* Tags */}
              <div className="flex flex-wrap gap-2 mb-3">
                {tags.map((tag) => (
                  <span
                    key={tag.id}
                    className={`px-3 py-1 text-xs font-medium rounded-full ${tagColors[tag.color] || tagColors.blue} flex items-center gap-2`}
                  >
                    {tag.tag_name}
                    <button
                      onClick={() => handleRemoveTag(tag.id)}
                      className="hover:opacity-70"
                      title="Remove tag"
                    >
                      ×
                    </button>
                  </span>
                ))}
              </div>

              {/* Status badges */}
              <div className="flex gap-2">
                <span
                  className={`px-2 py-1 text-xs font-medium rounded-full ${
                    user.is_active ? 'bg-green-100 text-green-800' : 'bg-gray-100 text-gray-800'
                  }`}
                >
                  {user.is_active ? 'Active' : 'Inactive'}
                </span>
                {user.is_banned && (
                  <span className="px-2 py-1 text-xs font-medium rounded-full bg-red-100 text-red-800">
                    Banned
                  </span>
                )}
                <span className="px-2 py-1 text-xs font-medium rounded-full bg-blue-100 text-blue-800">
                  {subscription.plan}
                </span>
              </div>
            </div>

            {/* Quick stats */}
            <div className="grid grid-cols-3 gap-4 text-center">
              <div className="bg-blue-50 rounded-lg p-3">
                <div className="text-2xl font-bold text-blue-600">{jobs.length}</div>
                <div className="text-xs text-gray-600">Total Jobs</div>
              </div>
              <div className="bg-green-50 rounded-lg p-3">
                <div className="text-2xl font-bold text-green-600">{usage.jobs_executed}/{usage.jobs_limit}</div>
                <div className="text-xs text-gray-600">Usage</div>
              </div>
              <div className="bg-purple-50 rounded-lg p-3">
                <div className="text-2xl font-bold text-purple-600">${(subscription.price_cents / 100).toFixed(0)}</div>
                <div className="text-xs text-gray-600">Monthly</div>
              </div>
            </div>
          </div>
        </div>
      </div>

      <div className="max-w-7xl mx-auto px-4 py-6">
        {/* Tabs */}
        <div className="border-b border-gray-200 mb-6">
          <nav className="-mb-px flex space-x-8">
            {[
              { id: 'profile', label: 'Profile' },
              { id: 'jobs', label: `Jobs (${jobs.length})` },
              { id: 'payments', label: `Payments (${payment_history.length})` },
              { id: 'tickets', label: `Support (${support_tickets.length})` },
              { id: 'activity', label: `Activity (${activity_timeline.length})` },
              { id: 'notes', label: `Notes (${notes.length})` },
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

        {/* Profile Tab */}
        {activeTab === 'profile' && (
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            {/* User Info */}
            <div className="bg-white rounded-lg shadow p-6">
              <h3 className="text-lg font-semibold mb-4">User Information</h3>
              <dl className="space-y-3">
                <div>
                  <dt className="text-sm text-gray-600">User ID</dt>
                  <dd className="text-sm font-mono text-gray-900">{user.uid}</dd>
                </div>
                <div>
                  <dt className="text-sm text-gray-600">Email</dt>
                  <dd className="text-sm text-gray-900">{user.email}</dd>
                </div>
                <div>
                  <dt className="text-sm text-gray-600">Username</dt>
                  <dd className="text-sm text-gray-900">{user.username || '-'}</dd>
                </div>
                <div>
                  <dt className="text-sm text-gray-600">Created</dt>
                  <dd className="text-sm text-gray-900">{new Date(user.created_at).toLocaleString()}</dd>
                </div>
                {user.is_banned && user.ban_reason && (
                  <div>
                    <dt className="text-sm text-gray-600">Ban Reason</dt>
                    <dd className="text-sm text-red-600">{user.ban_reason}</dd>
                    <dd className="text-xs text-gray-500 mt-1">
                      Banned on {new Date(user.banned_at!).toLocaleString()}
                    </dd>
                  </div>
                )}
              </dl>
            </div>

            {/* Subscription Info */}
            <div className="bg-white rounded-lg shadow p-6">
              <h3 className="text-lg font-semibold mb-4">Subscription</h3>
              <dl className="space-y-3">
                <div>
                  <dt className="text-sm text-gray-600">Plan</dt>
                  <dd className="text-sm font-semibold text-gray-900">{subscription.plan}</dd>
                </div>
                <div>
                  <dt className="text-sm text-gray-600">Status</dt>
                  <dd>
                    <span className={`px-2 py-1 text-xs font-medium rounded-full ${
                      subscription.status === 'active' ? 'bg-green-100 text-green-800' : 'bg-gray-100 text-gray-800'
                    }`}>
                      {subscription.status}
                    </span>
                  </dd>
                </div>
                <div>
                  <dt className="text-sm text-gray-600">Price</dt>
                  <dd className="text-sm text-gray-900">${(subscription.price_cents / 100).toFixed(2)}/month</dd>
                </div>
                <div>
                  <dt className="text-sm text-gray-600">Jobs Limit</dt>
                  <dd className="text-sm text-gray-900">{subscription.monthly_jobs_limit} jobs/month</dd>
                </div>
                {subscription.current_period_start && (
                  <>
                    <div>
                      <dt className="text-sm text-gray-600">Period Start</dt>
                      <dd className="text-sm text-gray-900">
                        {new Date(subscription.current_period_start).toLocaleDateString()}
                      </dd>
                    </div>
                    <div>
                      <dt className="text-sm text-gray-600">Period End</dt>
                      <dd className="text-sm text-gray-900">
                        {new Date(subscription.current_period_end!).toLocaleDateString()}
                      </dd>
                    </div>
                  </>
                )}
                {subscription.stripe_customer_id && (
                  <div>
                    <dt className="text-sm text-gray-600">Stripe Customer ID</dt>
                    <dd className="text-sm font-mono text-gray-900">{subscription.stripe_customer_id}</dd>
                  </div>
                )}
              </dl>
            </div>

            {/* Usage Info */}
            <div className="bg-white rounded-lg shadow p-6">
              <h3 className="text-lg font-semibold mb-4">Usage This Month</h3>
              <div className="space-y-4">
                <div>
                  <div className="flex justify-between text-sm mb-2">
                    <span className="text-gray-600">Jobs Executed</span>
                    <span className="font-semibold">{usage.jobs_executed} / {usage.jobs_limit}</span>
                  </div>
                  <div className="w-full bg-gray-200 rounded-full h-2">
                    <div
                      className={`h-2 rounded-full ${
                        (usage.jobs_executed / usage.jobs_limit) > 0.9 ? 'bg-red-500' :
                        (usage.jobs_executed / usage.jobs_limit) > 0.7 ? 'bg-yellow-500' : 'bg-green-500'
                      }`}
                      style={{ width: `${Math.min((usage.jobs_executed / usage.jobs_limit) * 100, 100)}%` }}
                    ></div>
                  </div>
                </div>
                <div className="text-xs text-gray-500">
                  Month: {String(usage.month).substring(0, 4)}-{String(usage.month).substring(4)}
                </div>
              </div>
            </div>

            {/* Add Tag Form */}
            <div className="bg-white rounded-lg shadow p-6">
              <h3 className="text-lg font-semibold mb-4">Add Tag</h3>
              <div className="space-y-3">
                <input
                  type="text"
                  value={newTagName}
                  onChange={(e) => setNewTagName(e.target.value)}
                  placeholder="Tag name"
                  className="w-full px-3 py-2 border rounded text-sm"
                />
                <select
                  value={newTagColor}
                  onChange={(e) => setNewTagColor(e.target.value)}
                  className="w-full px-3 py-2 border rounded text-sm"
                >
                  <option value="blue">Blue</option>
                  <option value="green">Green</option>
                  <option value="red">Red</option>
                  <option value="yellow">Yellow</option>
                  <option value="purple">Purple</option>
                </select>
                <button
                  onClick={handleAddTag}
                  disabled={addingTag || !newTagName.trim()}
                  className="w-full px-4 py-2 bg-blue-500 text-white rounded text-sm hover:bg-blue-600 disabled:opacity-50"
                >
                  {addingTag ? 'Adding...' : 'Add Tag'}
                </button>
              </div>
            </div>
          </div>
        )}

        {/* Jobs Tab */}
        {activeTab === 'jobs' && (
          <div className="bg-white rounded-lg shadow overflow-hidden">
            <div className="overflow-x-auto">
              <table className="min-w-full divide-y divide-gray-200">
                <thead className="bg-gray-50">
                  <tr>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">Sample Name</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">Job ID</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">Status</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">Current Step</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">Created</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">Completed</th>
                  </tr>
                </thead>
                <tbody className="bg-white divide-y divide-gray-200">
                  {jobs.map((job) => (
                    <tr key={job.job_id} className="hover:bg-gray-50">
                      <td className="px-6 py-4 whitespace-nowrap text-sm font-medium text-gray-900">
                        {job.sample_name}
                      </td>
                      <td className="px-6 py-4 whitespace-nowrap text-sm font-mono text-gray-600">
                        {job.job_id.substring(0, 8)}...
                      </td>
                      <td className="px-6 py-4 whitespace-nowrap">
                        <span className={`px-2 py-1 text-xs font-medium rounded-full ${
                          job.status === 'completed' ? 'bg-green-100 text-green-800' :
                          job.status === 'running' ? 'bg-blue-100 text-blue-800' :
                          job.status === 'failed' ? 'bg-red-100 text-red-800' :
                          'bg-gray-100 text-gray-800'
                        }`}>
                          {job.status}
                        </span>
                      </td>
                      <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-600">
                        {job.current_step || '-'}
                      </td>
                      <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-600">
                        {new Date(job.created_at).toLocaleString()}
                      </td>
                      <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-600">
                        {job.completed_at ? new Date(job.completed_at).toLocaleString() : '-'}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
              {jobs.length === 0 && (
                <div className="text-center py-8 text-gray-500">No jobs found</div>
              )}
            </div>
          </div>
        )}

        {/* Payments Tab */}
        {activeTab === 'payments' && (
          <div className="bg-white rounded-lg shadow overflow-hidden">
            <div className="overflow-x-auto">
              <table className="min-w-full divide-y divide-gray-200">
                <thead className="bg-gray-50">
                  <tr>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">Date</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">Event Type</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">Amount</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">Status</th>
                  </tr>
                </thead>
                <tbody className="bg-white divide-y divide-gray-200">
                  {payment_history.map((payment) => (
                    <tr key={payment.id} className="hover:bg-gray-50">
                      <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-900">
                        {new Date(payment.created_at).toLocaleString()}
                      </td>
                      <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-600">
                        {payment.event_type}
                      </td>
                      <td className="px-6 py-4 whitespace-nowrap text-sm font-medium text-gray-900">
                        {payment.amount ? `$${payment.amount.toFixed(2)} ${payment.currency.toUpperCase()}` : '-'}
                      </td>
                      <td className="px-6 py-4 whitespace-nowrap">
                        <span className={`px-2 py-1 text-xs font-medium rounded-full ${
                          payment.processed ? 'bg-green-100 text-green-800' : 'bg-yellow-100 text-yellow-800'
                        }`}>
                          {payment.processed ? 'Processed' : 'Pending'}
                        </span>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
              {payment_history.length === 0 && (
                <div className="text-center py-8 text-gray-500">No payment history</div>
              )}
            </div>
          </div>
        )}

        {/* Support Tickets Tab */}
        {activeTab === 'tickets' && (
          <div className="bg-white rounded-lg shadow overflow-hidden">
            <div className="overflow-x-auto">
              <table className="min-w-full divide-y divide-gray-200">
                <thead className="bg-gray-50">
                  <tr>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">Subject</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">Status</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">Messages</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">Last Message</th>
                    <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">Created</th>
                  </tr>
                </thead>
                <tbody className="bg-white divide-y divide-gray-200">
                  {support_tickets.map((ticket) => (
                    <tr key={ticket.id} className="hover:bg-gray-50">
                      <td className="px-6 py-4 text-sm font-medium text-gray-900">
                        {ticket.subject}
                      </td>
                      <td className="px-6 py-4 whitespace-nowrap">
                        <span className={`px-2 py-1 text-xs font-medium rounded-full ${
                          ticket.status === 'open' ? 'bg-blue-100 text-blue-800' :
                          ticket.status === 'resolved' ? 'bg-green-100 text-green-800' :
                          'bg-gray-100 text-gray-800'
                        }`}>
                          {ticket.status}
                        </span>
                      </td>
                      <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-600">
                        {ticket.message_count}
                      </td>
                      <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-600">
                        {new Date(ticket.last_message_at).toLocaleString()}
                      </td>
                      <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-600">
                        {new Date(ticket.created_at).toLocaleString()}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
              {support_tickets.length === 0 && (
                <div className="text-center py-8 text-gray-500">No support tickets</div>
              )}
            </div>
          </div>
        )}

        {/* Activity Timeline Tab */}
        {activeTab === 'activity' && (
          <div className="bg-white rounded-lg shadow p-6">
            <h3 className="text-lg font-semibold mb-4">Activity Timeline</h3>
            <div className="space-y-4">
              {activity_timeline.map((activity) => (
                <div key={activity.id} className="border-l-2 border-blue-500 pl-4 pb-4">
                  <div className="flex items-start justify-between">
                    <div>
                      <div className="font-medium text-gray-900">{activity.action}</div>
                      {activity.resource_type && (
                        <div className="text-sm text-gray-600">
                          {activity.resource_type}: {activity.resource_id}
                        </div>
                      )}
                      {activity.metadata && (
                        <div className="text-xs text-gray-500 mt-1">
                          <pre className="whitespace-pre-wrap">{activity.metadata}</pre>
                        </div>
                      )}
                    </div>
                    <div className="text-xs text-gray-500 whitespace-nowrap ml-4">
                      {new Date(activity.created_at).toLocaleString()}
                    </div>
                  </div>
                </div>
              ))}
              {activity_timeline.length === 0 && (
                <div className="text-center py-8 text-gray-500">No activity recorded</div>
              )}
            </div>
          </div>
        )}

        {/* Notes Tab */}
        {activeTab === 'notes' && (
          <div className="space-y-6">
            {/* Add Note Form */}
            <div className="bg-white rounded-lg shadow p-6">
              <h3 className="text-lg font-semibold mb-4">Add Note</h3>
              <textarea
                value={newNote}
                onChange={(e) => setNewNote(e.target.value)}
                placeholder="Enter note..."
                rows={4}
                className="w-full px-3 py-2 border rounded text-sm"
              ></textarea>
              <button
                onClick={handleAddNote}
                disabled={addingNote || !newNote.trim()}
                className="mt-3 px-4 py-2 bg-blue-500 text-white rounded text-sm hover:bg-blue-600 disabled:opacity-50"
              >
                {addingNote ? 'Adding...' : 'Add Note'}
              </button>
            </div>

            {/* Notes List */}
            <div className="bg-white rounded-lg shadow p-6">
              <h3 className="text-lg font-semibold mb-4">Notes ({notes.length})</h3>
              <div className="space-y-4">
                {notes.map((note) => (
                  <div key={note.id} className="border-l-4 border-blue-500 pl-4 py-2">
                    <div className="text-sm text-gray-900 whitespace-pre-wrap">{note.note_text}</div>
                    <div className="text-xs text-gray-500 mt-2">
                      Added by {note.admin_id} on {new Date(note.created_at).toLocaleString()}
                      {note.updated_at !== note.created_at && (
                        <span> (edited {new Date(note.updated_at).toLocaleString()})</span>
                      )}
                    </div>
                  </div>
                ))}
                {notes.length === 0 && (
                  <div className="text-center py-8 text-gray-500">No notes yet</div>
                )}
              </div>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
