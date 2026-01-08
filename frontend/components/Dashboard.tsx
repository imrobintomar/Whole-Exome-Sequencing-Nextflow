'use client';

import { useState } from 'react';
import { User } from '@/lib/api';
import { Sidebar } from './sidebar';
import { DashboardHeader } from './dashboard-header';
import UploadForm from './UploadForm';
import JobList from './JobList';
import JobDetailsPage from './JobDetailsPage';
import DashboardOverview from './DashboardOverview';
import AnalyticsDashboard from './AnalyticsDashboard';
import GenePanelFilter from './GenePanelFilter';
import ACMGClassificationPage from './ACMGClassificationPage';
import IGVBrowserPage from './IGVBrowserPage';
import VariantVisualizationPage from './VariantVisualizationPage';
import GenePanelPage from './GenePanelPage';
import SupportChat from './SupportChat';

interface DashboardProps {
  user: User;
  onLogout: () => void;
}

type ViewType = 'overview' | 'upload' | 'jobs' | 'job-details' | 'analytics' | 'panels' | 'gene-panel' | 'acmg' | 'igv' | 'variants' | 'support';

export default function Dashboard({ user, onLogout }: DashboardProps) {
  const [currentView, setCurrentView] = useState<ViewType>('overview');
  const [refreshTrigger, setRefreshTrigger] = useState(0);
  const [selectedJobId, setSelectedJobId] = useState<string | null>(null);
  const [selectedJobForACMG, setSelectedJobForACMG] = useState<{ jobId: string; sampleName: string } | null>(null);
  const [selectedJobForIGV, setSelectedJobForIGV] = useState<{ jobId: string; sampleName: string } | null>(null);
  const [selectedJobForVariants, setSelectedJobForVariants] = useState<{ jobId: string; sampleName: string } | null>(null);
  const [selectedJobForGenePanel, setSelectedJobForGenePanel] = useState<string | null>(null);

  const handleJobSubmitted = () => {
    setCurrentView('jobs');
    setRefreshTrigger(prev => prev + 1);
  };

  const handleViewChange = (view: string) => {
    setCurrentView(view as ViewType);
  };

  const handleJobClick = (jobId: string) => {
    setSelectedJobId(jobId);
    setCurrentView('job-details');
  };

  const handleClassifyClick = (jobId: string, sampleName: string) => {
    setSelectedJobForACMG({ jobId, sampleName });
    setCurrentView('acmg');
  };

  const handleIGVClick = (jobId: string, sampleName: string) => {
    setSelectedJobForIGV({ jobId, sampleName });
    setCurrentView('igv');
  };

  const handleVariantsClick = (jobId: string, sampleName?: string) => {
    setSelectedJobForVariants({ jobId, sampleName: sampleName || 'Sample' });
    setCurrentView('variants');
  };

  const handleGenePanelClick = (jobId: string) => {
    setSelectedJobForGenePanel(jobId);
    setCurrentView('gene-panel');
  };

  const handleBackToJobs = () => {
    setCurrentView('jobs');
    setSelectedJobId(null);
    setSelectedJobForACMG(null);
    setSelectedJobForIGV(null);
    setSelectedJobForVariants(null);
    setSelectedJobForGenePanel(null);
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
            <JobList
              key={refreshTrigger}
              onJobClick={handleJobClick}
              onClassifyClick={handleClassifyClick}
              onIGVClick={handleIGVClick}
              onVariantsClick={handleVariantsClick}
            />
          )}
          {currentView === 'job-details' && selectedJobId && (
            <JobDetailsPage
              jobId={selectedJobId}
              onBack={handleBackToJobs}
              onClassifyClick={handleClassifyClick}
              onIGVClick={handleIGVClick}
              onVariantsClick={(jobId) => handleVariantsClick(jobId)}
              onGenePanelClick={handleGenePanelClick}
            />
          )}
          {currentView === 'analytics' && (
            <AnalyticsDashboard />
          )}
          {currentView === 'panels' && (
            <div className="space-y-6">
              <GenePanelFilter />
              <JobList
                key={refreshTrigger}
                onJobClick={handleJobClick}
                onClassifyClick={handleClassifyClick}
                onIGVClick={handleIGVClick}
                onVariantsClick={handleVariantsClick}
              />
            </div>
          )}
          {currentView === 'acmg' && selectedJobForACMG && (
            <ACMGClassificationPage
              jobId={selectedJobForACMG.jobId}
              sampleName={selectedJobForACMG.sampleName}
              onBack={handleBackToJobs}
            />
          )}
          {currentView === 'igv' && selectedJobForIGV && (
            <IGVBrowserPage
              jobId={selectedJobForIGV.jobId}
              sampleName={selectedJobForIGV.sampleName}
              onBack={handleBackToJobs}
            />
          )}
          {currentView === 'variants' && selectedJobForVariants && (
            <VariantVisualizationPage
              jobId={selectedJobForVariants.jobId}
              sampleName={selectedJobForVariants.sampleName}
              onBack={handleBackToJobs}
            />
          )}
          {currentView === 'gene-panel' && selectedJobForGenePanel && (
            <GenePanelPage
              jobId={selectedJobForGenePanel}
              onBack={handleBackToJobs}
            />
          )}
          {currentView === 'support' && (
            <SupportChat />
          )}
        </main>
      </div>
    </div>
  );
}
