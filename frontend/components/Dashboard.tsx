'use client';

import { useState } from 'react';
import { User } from '@/lib/api';
import { Sidebar } from './sidebar';
import { DashboardHeader } from './dashboard-header';
import UploadForm from './UploadForm';
import JobList from './JobList';
import DashboardOverview from './DashboardOverview';
import AnalyticsDashboard from './AnalyticsDashboard';
import GenePanelFilter from './GenePanelFilter';
import ACMGClassificationView from './ACMGClassificationView';
import IGVBrowser from './IGVBrowser';
import { VariantVisualization } from './VariantVisualization';

interface DashboardProps {
  user: User;
  onLogout: () => void;
}

export default function Dashboard({ user, onLogout }: DashboardProps) {
  const [currentView, setCurrentView] = useState<'overview' | 'upload' | 'jobs' | 'analytics' | 'panels' | 'acmg' | 'igv' | 'variants'>('overview');
  const [refreshTrigger, setRefreshTrigger] = useState(0);
  const [selectedJobForACMG, setSelectedJobForACMG] = useState<{ jobId: string; sampleName: string } | null>(null);
  const [selectedJobForIGV, setSelectedJobForIGV] = useState<{ jobId: string; sampleName: string } | null>(null);
  const [selectedJobForVariants, setSelectedJobForVariants] = useState<string | null>(null);

  const handleJobSubmitted = () => {
    setCurrentView('jobs');
    setRefreshTrigger(prev => prev + 1);
  };

  const handleViewChange = (view: string) => {
    setCurrentView(view as any);
  };

  const handleClassifyClick = (jobId: string, sampleName: string) => {
    setSelectedJobForACMG({ jobId, sampleName });
    setCurrentView('acmg');
  };

  const handleIGVClick = (jobId: string, sampleName: string) => {
    setSelectedJobForIGV({ jobId, sampleName });
    setCurrentView('igv');
  };

  const handleVariantsClick = (jobId: string) => {
    setSelectedJobForVariants(jobId);
    setCurrentView('variants');
  };

  const handleBackToJobs = () => {
    setCurrentView('jobs');
    setSelectedJobForVariants(null);
  };

  return (
    <div className="flex h-screen overflow-hidden">
      <Sidebar currentView={currentView} onViewChange={handleViewChange} />

      <div className="flex flex-1 flex-col overflow-hidden">
        <DashboardHeader
          user={user}
          onLogout={onLogout}
          currentView={currentView}
          onViewChange={handleViewChange}
        />

        <main className="flex-1 overflow-y-auto bg-background p-4 sm:p-6">
          {currentView === 'overview' && (
            <DashboardOverview />
          )}
          {currentView === 'upload' && (
            <UploadForm onJobSubmitted={handleJobSubmitted} />
          )}
          {currentView === 'jobs' && (
            <JobList key={refreshTrigger} onClassifyClick={handleClassifyClick} onIGVClick={handleIGVClick} onVariantsClick={handleVariantsClick} />
          )}
          {currentView === 'analytics' && (
            <AnalyticsDashboard />
          )}
          {currentView === 'panels' && (
            <div className="space-y-6">
              <GenePanelFilter />
              <JobList key={refreshTrigger} onClassifyClick={handleClassifyClick} onIGVClick={handleIGVClick} onVariantsClick={handleVariantsClick} />
            </div>
          )}
          {currentView === 'acmg' && selectedJobForACMG && (
            <ACMGClassificationView
              jobId={selectedJobForACMG.jobId}
              sampleName={selectedJobForACMG.sampleName}
            />
          )}
          {currentView === 'igv' && selectedJobForIGV && (
            <IGVBrowser
              jobId={selectedJobForIGV.jobId}
              sampleName={selectedJobForIGV.sampleName}
            />
          )}
          {currentView === 'variants' && selectedJobForVariants && (
            <VariantVisualization
              jobId={selectedJobForVariants}
              onBack={handleBackToJobs}
            />
          )}
        </main>
      </div>
    </div>
  );
}
