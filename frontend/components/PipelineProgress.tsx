"use client"

import { CheckCircle2, Circle, Loader2, XCircle } from "lucide-react"
import { Job } from "@/lib/api"

interface PipelineProgressProps {
  job: Job
}

const PIPELINE_STEPS = [
  "Initializing",
  "Quality Control",
  "Alignment",
  "Sorting",
  "Quality Statistics",
  "Deduplication",
  "Post-dedup Sorting",
  "Base Recalibration",
  "Applying BQSR",
  "Variant Calling",
  "ANNOVAR Annotation",
  "Filtering & Extraction",
  "Completed"
]

export default function PipelineProgress({ job }: PipelineProgressProps) {
  const getStepStatus = (stepName: string) => {
    if (job.status === "failed") {
      const currentStepIndex = PIPELINE_STEPS.indexOf(job.current_step || "")
      const stepIndex = PIPELINE_STEPS.indexOf(stepName)
      if (stepIndex < currentStepIndex) return "completed"
      if (stepIndex === currentStepIndex) return "failed"
      return "pending"
    }

    if (job.status === "completed") {
      return "completed"
    }

    if (job.status === "pending") {
      return "pending"
    }

    // For running jobs
    const currentStepIndex = PIPELINE_STEPS.indexOf(job.current_step || "Initializing")
    const stepIndex = PIPELINE_STEPS.indexOf(stepName)

    if (stepIndex < currentStepIndex) return "completed"
    if (stepIndex === currentStepIndex) return "running"
    return "pending"
  }

  const getStepIcon = (stepName: string) => {
    const status = getStepStatus(stepName)

    switch (status) {
      case "completed":
        return <CheckCircle2 className="h-4 w-4 text-green-500" />
      case "running":
        return <Loader2 className="h-4 w-4 text-blue-500 animate-spin" />
      case "failed":
        return <XCircle className="h-4 w-4 text-red-500" />
      default:
        return <Circle className="h-4 w-4 text-gray-300" />
    }
  }

  const getStepColor = (stepName: string) => {
    const status = getStepStatus(stepName)

    switch (status) {
      case "completed":
        return "text-green-700 dark:text-green-400"
      case "running":
        return "text-blue-700 dark:text-blue-400 font-semibold"
      case "failed":
        return "text-red-700 dark:text-red-400"
      default:
        return "text-gray-500 dark:text-gray-400"
    }
  }

  // Show simplified view for completed or failed jobs
  if (job.status === "completed" || job.status === "failed") {
    return (
      <div className="flex items-center gap-2 text-sm">
        {job.status === "completed" ? (
          <>
            <CheckCircle2 className="h-4 w-4 text-green-500" />
            <span className="text-green-700 dark:text-green-400">All steps completed</span>
          </>
        ) : (
          <>
            <XCircle className="h-4 w-4 text-red-500" />
            <span className="text-red-700 dark:text-red-400">
              Failed at: {job.current_step || "Unknown step"}
            </span>
          </>
        )}
      </div>
    )
  }

  // Show current step for pending/running jobs
  return (
    <div className="space-y-2">
      <div className="flex items-center gap-2 text-sm">
        {job.status === "pending" ? (
          <>
            <Circle className="h-4 w-4 text-gray-300" />
            <span className="text-gray-500 dark:text-gray-400">Waiting to start...</span>
          </>
        ) : (
          <>
            <Loader2 className="h-4 w-4 text-blue-500 animate-spin" />
            <span className="text-blue-700 dark:text-blue-400 font-medium">
              {job.current_step || "Processing..."}
            </span>
          </>
        )}
      </div>

      {/* Detailed step list (expandable) */}
      <details className="text-xs">
        <summary className="cursor-pointer text-muted-foreground hover:text-foreground">
          View all steps
        </summary>
        <div className="mt-2 space-y-1 pl-2 border-l-2 border-gray-200 dark:border-gray-700">
          {PIPELINE_STEPS.map((step) => (
            <div key={step} className="flex items-center gap-2 py-0.5">
              {getStepIcon(step)}
              <span className={getStepColor(step)}>{step}</span>
            </div>
          ))}
        </div>
      </details>
    </div>
  )
}
