'use client';

import { useState } from 'react';
import { User } from '@/lib/api';
import { Sidebar } from './sidebar';
import { DashboardHeader } from './dashboard-header';
import UploadForm from './UploadForm';
import JobList from './JobList';
import DashboardOverview from './DashboardOverview';
import AnalyticsDashboard from './AnalyticsDashboard';

interface DashboardProps {
  user: User;
  onLogout: () => void;
}

export default function Dashboard({ user, onLogout }: DashboardProps) {
  const [currentView, setCurrentView] = useState<'overview' | 'upload' | 'jobs' | 'analytics'>('overview');
  const [refreshTrigger, setRefreshTrigger] = useState(0);

  const handleJobSubmitted = () => {
    setCurrentView('jobs');
    setRefreshTrigger(prev => prev + 1);
  };

  return (
    <div className="flex h-screen overflow-hidden">
      <Sidebar currentView={currentView} onViewChange={setCurrentView} />

      <div className="flex flex-1 flex-col overflow-hidden">
        <DashboardHeader user={user} onLogout={onLogout} />

        <main className="flex-1 overflow-y-auto bg-background p-6">
          {currentView === 'overview' && (
            <DashboardOverview />
          )}
          {currentView === 'upload' && (
            <UploadForm onJobSubmitted={handleJobSubmitted} />
          )}
          {currentView === 'jobs' && (
            <JobList key={refreshTrigger} />
          )}
          {currentView === 'analytics' && (
            <AnalyticsDashboard />
          )}
        </main>
      </div>
    </div>
  );
}
