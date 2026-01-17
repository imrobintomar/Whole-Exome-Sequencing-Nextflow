"use client"

import { useEffect, useState } from "react"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "./ui/card"
import { StatsCard } from "./ui/stats-card"
import { SkeletonStatsCard, Skeleton } from "./ui/skeleton"
import {
  Dna,
  Clock,
  CheckCircle2,
  XCircle,
  Activity,
  TrendingUp,
  ArrowRight,
  Sparkles,
  Upload,
  BarChart3
} from "lucide-react"
import { getJobs, Job } from "@/lib/api"
import { cn } from "@/lib/utils"
import { Button } from "./ui/button"
import { MobileStatsCarousel, DesktopStatsGrid } from "./ui/swipeable-cards"

export default function DashboardOverview() {
  const [jobs, setJobs] = useState<Job[]>([])
  const [loading, setLoading] = useState(true)

  useEffect(() => {
    loadJobs()
  }, [])

  const loadJobs = async () => {
    try {
      const data = await getJobs()
      setJobs(Array.isArray(data) ? data : [])
    } catch (error) {
      console.error("Failed to load jobs:", error)
      setJobs([])
    } finally {
      setLoading(false)
    }
  }

  const stats = {
    total: jobs.length,
    completed: jobs.filter((j) => j.status === "completed").length,
    running: jobs.filter((j) => j.status === "running").length,
    failed: jobs.filter((j) => j.status === "failed").length,
    pending: jobs.filter((j) => j.status === "pending").length,
  }

  // Calculate trends (mock for now - would need historical data)
  const trends = {
    total: { value: 12, direction: 'up' as const },
    completed: { value: 8, direction: 'up' as const },
    active: { value: 2, direction: 'neutral' as const },
    failed: { value: 50, direction: 'down' as const },
  }

  const recentJobs = jobs.slice(0, 5)

  const getStatusStyles = (status: string) => {
    switch (status) {
      case "completed":
        return {
          badge: "status-badge-completed",
          dot: "status-dot-completed",
          icon: <CheckCircle2 className="h-3.5 w-3.5" />,
        }
      case "running":
        return {
          badge: "status-badge-running",
          dot: "status-dot-running",
          icon: <Activity className="h-3.5 w-3.5 animate-pulse" />,
        }
      case "failed":
        return {
          badge: "status-badge-failed",
          dot: "status-dot-failed",
          icon: <XCircle className="h-3.5 w-3.5" />,
        }
      default:
        return {
          badge: "status-badge-pending",
          dot: "status-dot-pending",
          icon: <Clock className="h-3.5 w-3.5" />,
        }
    }
  }

  const formatDate = (dateString: string) => {
    const date = new Date(dateString)
    const now = new Date()
    const diffMs = now.getTime() - date.getTime()
    const diffMins = Math.floor(diffMs / 60000)
    const diffHours = Math.floor(diffMs / 3600000)
    const diffDays = Math.floor(diffMs / 86400000)

    if (diffMins < 60) return `${diffMins}m ago`
    if (diffHours < 24) return `${diffHours}h ago`
    if (diffDays < 7) return `${diffDays}d ago`
    return date.toLocaleDateString()
  }

  return (
    <div className="space-y-8 animate-fade-slide-up">
      {/* Welcome Section */}
      <div className="relative overflow-hidden rounded-2xl bg-gradient-to-br from-twilight-800/10 via-seagreen-500/5 to-mint-500/10 dark:from-twilight-800/20 dark:via-seagreen-500/10 dark:to-mint-500/20 p-8 border border-seagreen-500/20">
        <div className="absolute -right-20 -top-20 h-60 w-60 rounded-full bg-gradient-to-br from-seagreen-500/20 to-mint-500/20 blur-3xl" />
        <div className="absolute -left-20 -bottom-20 h-40 w-40 rounded-full bg-gradient-to-br from-mint-500/20 to-twilight-800/20 blur-3xl" />

        <div className="relative">
          <div className="flex items-center gap-2 mb-2">
            <Sparkles className="h-5 w-5 text-mint-500" />
            <span className="text-sm font-medium text-seagreen-600 dark:text-mint-400">Dashboard Overview</span>
          </div>
          <h1 className="text-3xl font-bold tracking-tight mb-2">
            Welcome to <span className="bg-gradient-to-r from-twilight-800 via-seagreen-500 to-mint-500 dark:from-seagreen-400 dark:via-mint-400 dark:to-mint-300 bg-clip-text text-transparent">ATGCFlow</span>
          </h1>
          <p className="text-muted-foreground max-w-lg">
            Track your Whole Exome Sequencing Analysis jobs at a glance. Upload samples, monitor progress, and download results.
          </p>

          <div className="flex gap-3 mt-6">
            <Button className="bg-gradient-to-r from-twilight-800 via-seagreen-500 to-mint-500 hover:opacity-90 text-white border-0 btn-press">
              <Upload className="h-4 w-4 mr-2" />
              New Analysis
            </Button>
            <Button variant="outline" className="rounded-xl border-seagreen-500/30 hover:bg-seagreen-500/10">
              <BarChart3 className="h-4 w-4 mr-2" />
              View Analytics
            </Button>
          </div>
        </div>
      </div>

      {/* Stats Cards - Desktop Grid */}
      <DesktopStatsGrid>
        {loading ? (
          <>
            <SkeletonStatsCard />
            <SkeletonStatsCard />
            <SkeletonStatsCard />
            <SkeletonStatsCard />
          </>
        ) : (
          <>
            <StatsCard
              title="Total Jobs"
              value={stats.total}
              icon={Dna}
              gradient="purple"
              trend={trends.total}
              description="All submissions"
            />
            <StatsCard
              title="Completed"
              value={stats.completed}
              icon={CheckCircle2}
              gradient="green"
              trend={trends.completed}
              description="Successfully finished"
            />
            <StatsCard
              title="Active"
              value={stats.running + stats.pending}
              icon={Activity}
              gradient="blue"
              trend={trends.active}
              description="Running or pending"
            />
            <StatsCard
              title="Failed"
              value={stats.failed}
              icon={XCircle}
              gradient="red"
              trend={trends.failed}
              description="Requires attention"
            />
          </>
        )}
      </DesktopStatsGrid>

      {/* Stats Cards - Mobile Swipeable Carousel */}
      {!loading && (
        <MobileStatsCarousel>
          {[
            <StatsCard
              key="total"
              title="Total Jobs"
              value={stats.total}
              icon={Dna}
              gradient="purple"
              trend={trends.total}
              description="All submissions"
            />,
            <StatsCard
              key="completed"
              title="Completed"
              value={stats.completed}
              icon={CheckCircle2}
              gradient="green"
              trend={trends.completed}
              description="Successfully finished"
            />,
            <StatsCard
              key="active"
              title="Active"
              value={stats.running + stats.pending}
              icon={Activity}
              gradient="blue"
              trend={trends.active}
              description="Running or pending"
            />,
            <StatsCard
              key="failed"
              title="Failed"
              value={stats.failed}
              icon={XCircle}
              gradient="red"
              trend={trends.failed}
              description="Requires attention"
            />,
          ]}
        </MobileStatsCarousel>
      )}

      {/* Recent Jobs */}
      <Card className="overflow-hidden border-0 shadow-lg bg-card/50 backdrop-blur-sm">
        <CardHeader className="border-b bg-muted/30">
          <div className="flex items-center justify-between">
            <div>
              <CardTitle className="text-lg">Recent Jobs</CardTitle>
              <CardDescription>Your latest analysis submissions</CardDescription>
            </div>
            <Button variant="ghost" size="sm" className="text-primary hover:text-primary">
              View all
              <ArrowRight className="h-4 w-4 ml-1" />
            </Button>
          </div>
        </CardHeader>
        <CardContent className="p-0">
          {loading ? (
            <div className="p-6 space-y-4">
              {[1, 2, 3, 4, 5].map((i) => (
                <div key={i} className="flex items-center gap-4">
                  <Skeleton variant="circular" width={40} height={40} />
                  <div className="flex-1 space-y-2">
                    <Skeleton variant="text" width="40%" height={16} />
                    <Skeleton variant="text" width="25%" height={12} />
                  </div>
                  <Skeleton variant="rectangular" width={80} height={24} className="rounded-full" />
                </div>
              ))}
            </div>
          ) : recentJobs.length === 0 ? (
            <div className="flex flex-col items-center justify-center py-16 text-center">
              <div className="relative mb-4">
                <div className="absolute inset-0 bg-gradient-to-br from-seagreen-500/20 to-mint-500/20 rounded-full blur-xl" />
                <div className="relative flex h-20 w-20 items-center justify-center rounded-full bg-gradient-to-br from-twilight-800/10 to-seagreen-500/10 border border-seagreen-500/20">
                  <Dna className="h-10 w-10 text-seagreen-500" />
                </div>
              </div>
              <h3 className="text-lg font-semibold mb-1">No jobs yet</h3>
              <p className="text-sm text-muted-foreground mb-4">
                Start by uploading your FASTQ files to begin analysis
              </p>
              <Button className="bg-gradient-to-r from-twilight-800 via-seagreen-500 to-mint-500 text-white border-0">
                <Upload className="h-4 w-4 mr-2" />
                Upload Files
              </Button>
            </div>
          ) : (
            <div className="divide-y divide-border/50">
              {recentJobs.map((job, index) => {
                const statusStyles = getStatusStyles(job.status)
                return (
                  <div
                    key={job.id}
                    className={cn(
                      "flex items-center gap-4 p-4 hover:bg-muted/30 transition-colors cursor-pointer",
                      "animate-fade-slide-up"
                    )}
                    style={{ animationDelay: `${index * 50}ms` }}
                  >
                    {/* Sample Icon */}
                    <div className={cn(
                      "flex h-10 w-10 items-center justify-center rounded-xl",
                      "bg-gradient-to-br from-twilight-800/10 to-seagreen-500/10",
                      "border border-seagreen-500/20"
                    )}>
                      <Dna className="h-5 w-5 text-seagreen-500" />
                    </div>

                    {/* Job Info */}
                    <div className="flex-1 min-w-0">
                      <p className="font-medium truncate">{job.sample_name}</p>
                      <div className="flex items-center gap-2 text-xs text-muted-foreground">
                        <Clock className="h-3 w-3" />
                        <span>{formatDate(job.created_at)}</span>
                        {job.current_step && (
                          <>
                            <span className="text-border">•</span>
                            <span className="truncate">{job.current_step}</span>
                          </>
                        )}
                      </div>
                    </div>

                    {/* Status Badge */}
                    <div className={cn("status-badge", statusStyles.badge)}>
                      {statusStyles.icon}
                      <span className="capitalize">{job.status}</span>
                    </div>

                    {/* Arrow */}
                    <ArrowRight className="h-4 w-4 text-muted-foreground" />
                  </div>
                )
              })}
            </div>
          )}
        </CardContent>
      </Card>

      {/* Quick Actions */}
      <div className="grid gap-4 md:grid-cols-3">
        <Card className="card-hover cursor-pointer border-dashed border-2 hover:border-twilight-800/50 transition-colors group">
          <CardContent className="flex flex-col items-center justify-center py-8">
            <div className="mb-4 flex h-14 w-14 items-center justify-center rounded-2xl bg-twilight-800/10 group-hover:bg-twilight-800/20 transition-colors">
              <Upload className="h-7 w-7 text-twilight-800 dark:text-seagreen-400" />
            </div>
            <h3 className="font-semibold mb-1">Upload Samples</h3>
            <p className="text-sm text-muted-foreground text-center">
              Start a new WES analysis
            </p>
          </CardContent>
        </Card>

        <Card className="card-hover cursor-pointer border-dashed border-2 hover:border-seagreen-500/50 transition-colors group">
          <CardContent className="flex flex-col items-center justify-center py-8">
            <div className="mb-4 flex h-14 w-14 items-center justify-center rounded-2xl bg-seagreen-500/10 group-hover:bg-seagreen-500/20 transition-colors">
              <BarChart3 className="h-7 w-7 text-seagreen-500" />
            </div>
            <h3 className="font-semibold mb-1">View Analytics</h3>
            <p className="text-sm text-muted-foreground text-center">
              Explore job statistics
            </p>
          </CardContent>
        </Card>

        <Card className="card-hover cursor-pointer border-dashed border-2 hover:border-mint-500/50 transition-colors group">
          <CardContent className="flex flex-col items-center justify-center py-8">
            <div className="mb-4 flex h-14 w-14 items-center justify-center rounded-2xl bg-mint-500/10 group-hover:bg-mint-500/20 transition-colors">
              <TrendingUp className="h-7 w-7 text-mint-600 dark:text-mint-400" />
            </div>
            <h3 className="font-semibold mb-1">Gene Panels</h3>
            <p className="text-sm text-muted-foreground text-center">
              Filter by gene panels
            </p>
          </CardContent>
        </Card>
      </div>
    </div>
  )
}
