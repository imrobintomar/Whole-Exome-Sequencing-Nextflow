"use client"

import * as React from "react"
import { cn } from "@/lib/utils"

interface SuccessGaugeProps {
  value: number // 0-100
  label?: string
  size?: "sm" | "md" | "lg"
  showPercentage?: boolean
  className?: string
  animate?: boolean
}

export function SuccessGauge({
  value,
  label = "Success Rate",
  size = "md",
  showPercentage = true,
  className,
  animate = true
}: SuccessGaugeProps) {
  const [displayValue, setDisplayValue] = React.useState(0)
  const clampedValue = Math.max(0, Math.min(100, value))

  // Animate value on mount/change
  React.useEffect(() => {
    if (!animate) {
      setDisplayValue(clampedValue)
      return
    }

    const duration = 1000
    const steps = 60
    const stepValue = clampedValue / steps
    let current = 0

    const interval = setInterval(() => {
      current += stepValue
      if (current >= clampedValue) {
        setDisplayValue(clampedValue)
        clearInterval(interval)
      } else {
        setDisplayValue(Math.round(current))
      }
    }, duration / steps)

    return () => clearInterval(interval)
  }, [clampedValue, animate])

  // Size configurations
  const sizes = {
    sm: { container: "h-24 w-24", stroke: 6, fontSize: "text-lg", labelSize: "text-[10px]" },
    md: { container: "h-36 w-36", stroke: 8, fontSize: "text-2xl", labelSize: "text-xs" },
    lg: { container: "h-48 w-48", stroke: 10, fontSize: "text-4xl", labelSize: "text-sm" },
  }

  const config = sizes[size]
  const radius = 45
  const circumference = 2 * Math.PI * radius
  const strokeDashoffset = circumference - (displayValue / 100) * circumference

  // Color based on value
  const getGradientColors = () => {
    if (displayValue >= 80) return { start: "#00E897", end: "#007F4F" } // Mint to turf
    if (displayValue >= 60) return { start: "#00A0A0", end: "#00E897" } // Seagreen to mint
    if (displayValue >= 40) return { start: "#F2D513", end: "#00A0A0" } // Golden to seagreen
    return { start: "#EF4444", end: "#F2D513" } // Red to golden
  }

  const colors = getGradientColors()
  const gradientId = React.useId()

  return (
    <div className={cn("flex flex-col items-center gap-2", className)}>
      <div className={cn("relative", config.container)}>
        <svg className="w-full h-full -rotate-90" viewBox="0 0 100 100">
          {/* Gradient Definition */}
          <defs>
            <linearGradient id={gradientId} x1="0%" y1="0%" x2="100%" y2="0%">
              <stop offset="0%" stopColor={colors.start} />
              <stop offset="100%" stopColor={colors.end} />
            </linearGradient>
          </defs>

          {/* Background Circle */}
          <circle
            cx="50"
            cy="50"
            r={radius}
            fill="none"
            stroke="currentColor"
            strokeWidth={config.stroke}
            className="text-muted/30"
          />

          {/* Progress Circle */}
          <circle
            cx="50"
            cy="50"
            r={radius}
            fill="none"
            stroke={`url(#${gradientId})`}
            strokeWidth={config.stroke}
            strokeLinecap="round"
            strokeDasharray={circumference}
            strokeDashoffset={strokeDashoffset}
            className={cn(
              "transition-all duration-1000 ease-out",
              "drop-shadow-[0_0_8px_rgba(0,232,151,0.5)]"
            )}
          />
        </svg>

        {/* Center Content */}
        <div className="absolute inset-0 flex flex-col items-center justify-center">
          {showPercentage && (
            <span className={cn(
              "font-bold tabular-nums",
              config.fontSize,
              "bg-gradient-to-r from-seagreen-500 to-mint-500 bg-clip-text text-transparent"
            )}>
              {displayValue}%
            </span>
          )}
        </div>
      </div>

      {/* Label */}
      {label && (
        <span className={cn(
          "font-medium text-muted-foreground",
          config.labelSize
        )}>
          {label}
        </span>
      )}
    </div>
  )
}

// Mini gauge for inline use
export function MiniGauge({
  value,
  className
}: {
  value: number
  className?: string
}) {
  const clampedValue = Math.max(0, Math.min(100, value))

  return (
    <div className={cn("flex items-center gap-2", className)}>
      <div className="relative h-8 w-8">
        <svg className="w-full h-full -rotate-90" viewBox="0 0 32 32">
          <circle
            cx="16"
            cy="16"
            r="12"
            fill="none"
            stroke="currentColor"
            strokeWidth="3"
            className="text-muted/30"
          />
          <circle
            cx="16"
            cy="16"
            r="12"
            fill="none"
            stroke="url(#miniGradient)"
            strokeWidth="3"
            strokeLinecap="round"
            strokeDasharray={`${clampedValue * 0.754} 100`}
            className="transition-all duration-500"
          />
          <defs>
            <linearGradient id="miniGradient" x1="0%" y1="0%" x2="100%" y2="0%">
              <stop offset="0%" stopColor="#00A0A0" />
              <stop offset="100%" stopColor="#00E897" />
            </linearGradient>
          </defs>
        </svg>
      </div>
      <span className="text-sm font-semibold">{clampedValue}%</span>
    </div>
  )
}

// Success Gauge Card with additional stats
export function SuccessGaugeCard({
  successRate,
  totalJobs,
  completedJobs,
  failedJobs,
  className
}: {
  successRate: number
  totalJobs: number
  completedJobs: number
  failedJobs: number
  className?: string
}) {
  return (
    <div className={cn(
      "flex flex-col items-center p-6 rounded-2xl",
      "bg-gradient-to-br from-twilight-800/5 via-seagreen-500/5 to-mint-500/5",
      "border border-seagreen-500/20",
      className
    )}>
      <SuccessGauge value={successRate} size="md" />

      <div className="grid grid-cols-3 gap-4 mt-6 w-full">
        <div className="text-center">
          <p className="text-2xl font-bold text-foreground">{totalJobs}</p>
          <p className="text-xs text-muted-foreground">Total</p>
        </div>
        <div className="text-center">
          <p className="text-2xl font-bold text-mint-500">{completedJobs}</p>
          <p className="text-xs text-muted-foreground">Completed</p>
        </div>
        <div className="text-center">
          <p className="text-2xl font-bold text-red-500">{failedJobs}</p>
          <p className="text-xs text-muted-foreground">Failed</p>
        </div>
      </div>
    </div>
  )
}
