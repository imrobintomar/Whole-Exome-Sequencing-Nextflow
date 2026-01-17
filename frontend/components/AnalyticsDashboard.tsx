"use client"

import { useEffect, useState } from "react"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "./ui/card"
import { getJobs, Job } from "@/lib/api"
import { BarChart, Bar, PieChart, Pie, Cell, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from "recharts"
import { TrendingUp, Clock, CheckCircle2, XCircle, Loader, Sparkles } from "lucide-react"
import { TimeRangeSelector, TimeRange, getDateRangeFromTimeRange } from "./ui/time-range-selector"
import { SuccessGauge } from "./ui/success-gauge"
import { cn } from "@/lib/utils"

export default function AnalyticsDashboard() {
  const [jobs, setJobs] = useState<Job[]>([])
  const [loading, setLoading] = useState(true)
  const [timeRange, setTimeRange] = useState<TimeRange>("30d")

  useEffect(() => {
    loadJobs()
  }, [])

  const loadJobs = async () => {
    try {
      const data = await getJobs()
      setJobs(Array.isArray(data) ? data : [])
    } catch (error) {
      console.error("Failed to load jobs:", error)
      setJobs([]) // Ensure jobs is always an array
    } finally {
      setLoading(false)
    }
  }

  if (loading) {
    return (
      <div className="flex items-center justify-center py-12">
        <Loader className="h-8 w-8 animate-spin text-primary" />
      </div>
    )
  }

  // Filter jobs by time range
  const { start: rangeStart } = getDateRangeFromTimeRange(timeRange)
  const filteredJobs = jobs.filter((j) => {
    const jobDate = new Date(j.created_at)
    return jobDate >= rangeStart
  })

  // Status distribution for pie chart (using filtered jobs)
  const statusCounts = {
    pending: filteredJobs.filter((j) => j.status === "pending").length,
    running: filteredJobs.filter((j) => j.status === "running").length,
    completed: filteredJobs.filter((j) => j.status === "completed").length,
    failed: filteredJobs.filter((j) => j.status === "failed").length,
  }

  // Brand colors for pie chart
  const pieData = [
    { name: "Completed", value: statusCounts.completed, color: "#00E897" }, // mint
    { name: "Running", value: statusCounts.running, color: "#00A0A0" }, // seagreen
    { name: "Pending", value: statusCounts.pending, color: "#F2D513" }, // golden
    { name: "Failed", value: statusCounts.failed, color: "#ef4444" }, // red
  ].filter((item) => item.value > 0)

  // Jobs over time - dynamic based on time range
  const getDaysForRange = (range: TimeRange): number => {
    switch (range) {
      case "7d": return 7
      case "30d": return 30
      case "90d": return 14 // Show 14 data points for 90 days
      case "1y": return 12 // Show 12 months
      case "all": return 30
      default: return 7
    }
  }

  const daysToShow = getDaysForRange(timeRange)
  const datePoints = Array.from({ length: Math.min(daysToShow, 14) }, (_, i) => {
    const date = new Date()
    date.setDate(date.getDate() - (Math.min(daysToShow, 14) - 1 - i))
    return date.toISOString().split("T")[0]
  })

  const jobsOverTime = datePoints.map((date) => ({
    date: new Date(date).toLocaleDateString("en-US", { month: "short", day: "numeric" }),
    jobs: filteredJobs.filter((j) => j.created_at.startsWith(date)).length,
  }))

  // Success rate calculation (using filtered jobs)
  const totalFinished = statusCounts.completed + statusCounts.failed
  const successRate = totalFinished > 0 ? (statusCounts.completed / totalFinished) * 100 : 0
  const successRateDisplay = successRate.toFixed(1)

  // Average processing time (for completed jobs in range)
  const completedJobs = filteredJobs.filter((j) => j.status === "completed" && j.completed_at)
  const avgProcessingTime = completedJobs.length > 0
    ? completedJobs.reduce((acc, job) => {
        const start = new Date(job.created_at).getTime()
        const end = new Date(job.completed_at!).getTime()
        return acc + (end - start)
      }, 0) / completedJobs.length
    : 0

  const avgHours = (avgProcessingTime / (1000 * 60 * 60)).toFixed(1)

  return (
    <div className="space-y-6 animate-fade-slide-up">
      {/* Header with Time Range Selector */}
      <div className="relative overflow-hidden rounded-2xl bg-gradient-to-br from-twilight-800/10 via-seagreen-500/5 to-mint-500/10 dark:from-twilight-800/20 dark:via-seagreen-500/10 dark:to-mint-500/20 p-6 border border-seagreen-500/20">
        <div className="absolute -right-20 -top-20 h-40 w-40 rounded-full bg-gradient-to-br from-seagreen-500/20 to-mint-500/20 blur-3xl" />
        <div className="relative flex flex-col sm:flex-row sm:items-center sm:justify-between gap-4">
          <div>
            <div className="flex items-center gap-2 mb-1">
              <Sparkles className="h-5 w-5 text-mint-500" />
              <span className="text-sm font-medium text-seagreen-600 dark:text-mint-400">Pipeline Analytics</span>
            </div>
            <h2 className="text-2xl font-bold tracking-tight">Performance Insights</h2>
            <p className="text-muted-foreground text-sm">
              Track your WES pipeline metrics and trends
            </p>
          </div>
          <TimeRangeSelector
            value={timeRange}
            onChange={setTimeRange}
            variant="dropdown"
          />
        </div>
      </div>

      {/* Key Metrics - with Success Gauge */}
      <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-5">
        {/* Success Rate Gauge */}
        <Card className={cn(
          "lg:col-span-2 overflow-hidden",
          "bg-gradient-to-br from-mint-500/5 via-transparent to-seagreen-500/5",
          "border-mint-500/20"
        )}>
          <CardContent className="p-6">
            <div className="flex items-center gap-6">
              <SuccessGauge value={successRate} size={120} />
              <div className="flex-1">
                <h3 className="text-lg font-semibold mb-1">Success Rate</h3>
                <p className="text-sm text-muted-foreground mb-3">
                  {statusCounts.completed} of {totalFinished} jobs successful
                </p>
                <div className="flex items-center gap-2">
                  <div className="flex items-center gap-1.5">
                    <div className="h-2 w-2 rounded-full bg-mint-500" />
                    <span className="text-xs text-muted-foreground">Completed</span>
                  </div>
                  <div className="flex items-center gap-1.5">
                    <div className="h-2 w-2 rounded-full bg-red-500" />
                    <span className="text-xs text-muted-foreground">Failed</span>
                  </div>
                </div>
              </div>
            </div>
          </CardContent>
        </Card>

        <Card className={cn(
          "overflow-hidden",
          "bg-gradient-to-br from-seagreen-500/5 via-transparent to-transparent",
          "border-seagreen-500/20"
        )}>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Avg Processing</CardTitle>
            <div className="flex h-8 w-8 items-center justify-center rounded-lg bg-seagreen-500/10">
              <Clock className="h-4 w-4 text-seagreen-500" />
            </div>
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold">{avgHours}h</div>
            <p className="text-xs text-muted-foreground">
              {completedJobs.length} completed
            </p>
          </CardContent>
        </Card>

        <Card className={cn(
          "overflow-hidden",
          "bg-gradient-to-br from-golden-500/5 via-transparent to-transparent",
          "border-golden-500/20"
        )}>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Active Jobs</CardTitle>
            <div className="flex h-8 w-8 items-center justify-center rounded-lg bg-golden-500/10">
              <Loader className="h-4 w-4 text-golden-500 animate-spin" />
            </div>
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold">{statusCounts.running + statusCounts.pending}</div>
            <p className="text-xs text-muted-foreground">
              In pipeline
            </p>
          </CardContent>
        </Card>

        <Card className={cn(
          "overflow-hidden",
          "bg-gradient-to-br from-twilight-800/5 via-transparent to-transparent",
          "border-twilight-800/20"
        )}>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Total Jobs</CardTitle>
            <div className="flex h-8 w-8 items-center justify-center rounded-lg bg-twilight-800/10">
              <CheckCircle2 className="h-4 w-4 text-twilight-800 dark:text-seagreen-400" />
            </div>
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold">{filteredJobs.length}</div>
            <p className="text-xs text-muted-foreground">
              In selected range
            </p>
          </CardContent>
        </Card>
      </div>

      {/* Charts */}
      <div className="grid gap-4 md:grid-cols-2">
        <Card>
          <CardHeader>
            <CardTitle>Job Status Distribution</CardTitle>
            <CardDescription>
              Breakdown of all jobs by current status
            </CardDescription>
          </CardHeader>
          <CardContent>
            {pieData.length > 0 ? (
              <ResponsiveContainer width="100%" height={300}>
                <PieChart>
                  <Pie
                    data={pieData}
                    cx="50%"
                    cy="50%"
                    labelLine={false}
                    label={({ name, percent }) => `${name} ${((percent ?? 0) * 100).toFixed(0)}%`}
                    outerRadius={80}
                    fill="#8884d8"
                    dataKey="value"
                  >
                    {pieData.map((entry, index) => (
                      <Cell key={`cell-${index}`} fill={entry.color} />
                    ))}
                  </Pie>
                  <Tooltip />
                </PieChart>
              </ResponsiveContainer>
            ) : (
              <div className="flex h-[300px] items-center justify-center text-muted-foreground">
                No job data available
              </div>
            )}
          </CardContent>
        </Card>

        <Card className="border-seagreen-500/20">
          <CardHeader>
            <CardTitle>Jobs Over Time</CardTitle>
            <CardDescription>
              Job submissions in selected time range
            </CardDescription>
          </CardHeader>
          <CardContent>
            <ResponsiveContainer width="100%" height={300}>
              <BarChart data={jobsOverTime}>
                <CartesianGrid strokeDasharray="3 3" className="stroke-muted" />
                <XAxis dataKey="date" className="text-xs" />
                <YAxis className="text-xs" />
                <Tooltip
                  contentStyle={{
                    backgroundColor: "hsl(var(--card))",
                    border: "1px solid hsl(var(--border))",
                    borderRadius: "12px",
                  }}
                />
                <Bar dataKey="jobs" fill="#00A0A0" radius={[4, 4, 0, 0]} />
              </BarChart>
            </ResponsiveContainer>
          </CardContent>
        </Card>
      </div>

      {/* Additional Insights */}
      <Card className="border-seagreen-500/20 overflow-hidden">
        <CardHeader className="bg-gradient-to-r from-twilight-800/5 via-seagreen-500/5 to-mint-500/5">
          <CardTitle>Pipeline Insights</CardTitle>
          <CardDescription>
            Performance summary and recommendations
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-4 pt-6">
          <div className="grid gap-4 md:grid-cols-3">
            <div className="space-y-2 p-4 rounded-xl bg-mint-500/5 border border-mint-500/20">
              <div className="flex items-center gap-2">
                <div className="flex h-8 w-8 items-center justify-center rounded-lg bg-mint-500/20">
                  <CheckCircle2 className="h-4 w-4 text-mint-600 dark:text-mint-400" />
                </div>
                <h4 className="font-semibold">Reliability</h4>
              </div>
              <p className="text-sm text-muted-foreground">
                {successRateDisplay}% of finished jobs completed successfully
              </p>
            </div>
            <div className="space-y-2 p-4 rounded-xl bg-seagreen-500/5 border border-seagreen-500/20">
              <div className="flex items-center gap-2">
                <div className="flex h-8 w-8 items-center justify-center rounded-lg bg-seagreen-500/20">
                  <Clock className="h-4 w-4 text-seagreen-600 dark:text-seagreen-400" />
                </div>
                <h4 className="font-semibold">Efficiency</h4>
              </div>
              <p className="text-sm text-muted-foreground">
                Average job completion time: {avgHours} hours
              </p>
            </div>
            <div className="space-y-2 p-4 rounded-xl bg-red-500/5 border border-red-500/20">
              <div className="flex items-center gap-2">
                <div className="flex h-8 w-8 items-center justify-center rounded-lg bg-red-500/20">
                  <XCircle className="h-4 w-4 text-red-600 dark:text-red-400" />
                </div>
                <h4 className="font-semibold">Failures</h4>
              </div>
              <p className="text-sm text-muted-foreground">
                {statusCounts.failed} jobs failed and may need attention
              </p>
            </div>
          </div>
        </CardContent>
      </Card>
    </div>
  )
}
