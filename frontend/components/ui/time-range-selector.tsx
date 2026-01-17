"use client"

import * as React from "react"
import { cn } from "@/lib/utils"
import { Calendar, ChevronDown } from "lucide-react"

export type TimeRange = "7d" | "30d" | "90d" | "1y" | "all"

interface TimeRangeSelectorProps {
  value: TimeRange
  onChange: (value: TimeRange) => void
  className?: string
  variant?: "default" | "compact"
}

const timeRangeOptions: { value: TimeRange; label: string; shortLabel: string }[] = [
  { value: "7d", label: "Last 7 days", shortLabel: "7D" },
  { value: "30d", label: "Last 30 days", shortLabel: "30D" },
  { value: "90d", label: "Last 90 days", shortLabel: "90D" },
  { value: "1y", label: "Last year", shortLabel: "1Y" },
  { value: "all", label: "All time", shortLabel: "All" },
]

export function TimeRangeSelector({
  value,
  onChange,
  className,
  variant = "default"
}: TimeRangeSelectorProps) {
  const [isOpen, setIsOpen] = React.useState(false)
  const containerRef = React.useRef<HTMLDivElement>(null)

  // Close dropdown when clicking outside
  React.useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (containerRef.current && !containerRef.current.contains(event.target as Node)) {
        setIsOpen(false)
      }
    }

    document.addEventListener("mousedown", handleClickOutside)
    return () => document.removeEventListener("mousedown", handleClickOutside)
  }, [])

  const selectedOption = timeRangeOptions.find(opt => opt.value === value)

  if (variant === "compact") {
    return (
      <div className={cn("flex items-center gap-1 p-1 rounded-xl bg-muted/50", className)}>
        {timeRangeOptions.map((option) => (
          <button
            key={option.value}
            onClick={() => onChange(option.value)}
            className={cn(
              "px-3 py-1.5 rounded-lg text-xs font-medium transition-all duration-200",
              value === option.value
                ? "bg-gradient-to-r from-seagreen-500 to-mint-500 text-white shadow-sm"
                : "text-muted-foreground hover:text-foreground hover:bg-muted"
            )}
          >
            {option.shortLabel}
          </button>
        ))}
      </div>
    )
  }

  return (
    <div ref={containerRef} className={cn("relative", className)}>
      <button
        onClick={() => setIsOpen(!isOpen)}
        className={cn(
          "flex items-center gap-2 px-4 py-2 rounded-xl",
          "bg-white/80 dark:bg-twilight-800/80 backdrop-blur-sm",
          "border border-seagreen-500/20",
          "text-sm font-medium",
          "hover:border-seagreen-500/40",
          "transition-all duration-200",
          isOpen && "border-seagreen-500/40 ring-2 ring-seagreen-500/20"
        )}
      >
        <Calendar className="h-4 w-4 text-seagreen-500" />
        <span>{selectedOption?.label}</span>
        <ChevronDown className={cn(
          "h-4 w-4 text-muted-foreground transition-transform duration-200",
          isOpen && "rotate-180"
        )} />
      </button>

      {/* Dropdown */}
      {isOpen && (
        <div className={cn(
          "absolute right-0 top-full mt-2 z-50",
          "w-48 p-2 rounded-xl",
          "bg-white/95 dark:bg-twilight-800/95 backdrop-blur-xl",
          "border border-seagreen-500/20",
          "shadow-xl shadow-twilight-800/10",
          "animate-scale-in origin-top-right"
        )}>
          {timeRangeOptions.map((option) => (
            <button
              key={option.value}
              onClick={() => {
                onChange(option.value)
                setIsOpen(false)
              }}
              className={cn(
                "flex items-center justify-between w-full px-3 py-2 rounded-lg",
                "text-sm transition-colors duration-150",
                value === option.value
                  ? "bg-gradient-to-r from-twilight-800/10 to-seagreen-500/10 text-foreground font-medium"
                  : "text-muted-foreground hover:text-foreground hover:bg-muted/50"
              )}
            >
              <span>{option.label}</span>
              {value === option.value && (
                <div className="h-2 w-2 rounded-full bg-gradient-to-r from-seagreen-500 to-mint-500" />
              )}
            </button>
          ))}
        </div>
      )}
    </div>
  )
}

// Helper function to get date range from TimeRange
export function getDateRangeFromTimeRange(timeRange: TimeRange): { from: Date; to: Date } {
  const now = new Date()
  const to = new Date(now)
  let from: Date

  switch (timeRange) {
    case "7d":
      from = new Date(now)
      from.setDate(from.getDate() - 7)
      break
    case "30d":
      from = new Date(now)
      from.setDate(from.getDate() - 30)
      break
    case "90d":
      from = new Date(now)
      from.setDate(from.getDate() - 90)
      break
    case "1y":
      from = new Date(now)
      from.setFullYear(from.getFullYear() - 1)
      break
    case "all":
    default:
      from = new Date(0) // Beginning of time
      break
  }

  return { from, to }
}

// Format date range for display
export function formatDateRange(timeRange: TimeRange): string {
  const { from, to } = getDateRangeFromTimeRange(timeRange)

  if (timeRange === "all") {
    return "All time"
  }

  const formatDate = (date: Date) => {
    return date.toLocaleDateString("en-US", { month: "short", day: "numeric" })
  }

  return `${formatDate(from)} - ${formatDate(to)}`
}
