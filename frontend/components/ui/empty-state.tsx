"use client"

import * as React from "react"
import {
  Inbox,
  Search,
  FileX,
  Database,
  Upload,
  FolderOpen,
  AlertCircle,
  Wifi,
  Lock,
  Clock,
  Dna,
  BarChart3,
  Filter,
  LucideIcon,
} from "lucide-react"
import { Button } from "./button"
import { cn } from "@/lib/utils"

// Preset empty state types
type EmptyStatePreset =
  | "no-data"
  | "no-results"
  | "no-files"
  | "no-jobs"
  | "no-analytics"
  | "no-panels"
  | "error"
  | "offline"
  | "unauthorized"
  | "coming-soon"
  | "custom"

interface EmptyStateAction {
  label: string
  onClick: () => void
  variant?: "default" | "outline" | "ghost"
  icon?: LucideIcon
}

interface EmptyStateProps {
  preset?: EmptyStatePreset
  icon?: LucideIcon
  title?: string
  description?: string
  action?: EmptyStateAction
  secondaryAction?: EmptyStateAction
  className?: string
  size?: "sm" | "md" | "lg"
  children?: React.ReactNode
}

// Preset configurations
const presets: Record<
  EmptyStatePreset,
  {
    icon: LucideIcon
    title: string
    description: string
    gradient: string
    iconColor: string
  }
> = {
  "no-data": {
    icon: Database,
    title: "No data yet",
    description: "There's no data to display at the moment. Start by adding some content.",
    gradient: "from-seagreen-500/10 to-mint-500/10",
    iconColor: "text-seagreen-500",
  },
  "no-results": {
    icon: Search,
    title: "No results found",
    description: "We couldn't find any matches for your search. Try adjusting your filters.",
    gradient: "from-twilight-800/10 to-seagreen-500/10",
    iconColor: "text-twilight-800 dark:text-seagreen-400",
  },
  "no-files": {
    icon: FileX,
    title: "No files found",
    description: "This folder is empty. Upload files to get started.",
    gradient: "from-golden-500/10 to-mint-500/10",
    iconColor: "text-golden-600",
  },
  "no-jobs": {
    icon: Dna,
    title: "No analysis jobs yet",
    description: "Start your first WES analysis by uploading FASTQ files.",
    gradient: "from-twilight-800/10 via-seagreen-500/10 to-mint-500/10",
    iconColor: "text-seagreen-500",
  },
  "no-analytics": {
    icon: BarChart3,
    title: "No analytics data",
    description: "Analytics will appear once you have completed jobs to analyze.",
    gradient: "from-mint-500/10 to-seagreen-500/10",
    iconColor: "text-mint-500",
  },
  "no-panels": {
    icon: Filter,
    title: "No gene panels configured",
    description: "Create custom gene panels for targeted variant filtering.",
    gradient: "from-seagreen-500/10 to-turf-500/10",
    iconColor: "text-seagreen-500",
  },
  error: {
    icon: AlertCircle,
    title: "Something went wrong",
    description: "We encountered an error loading this content. Please try again.",
    gradient: "from-red-500/10 to-red-500/5",
    iconColor: "text-red-500",
  },
  offline: {
    icon: Wifi,
    title: "You're offline",
    description: "Check your internet connection and try again.",
    gradient: "from-amber-500/10 to-amber-500/5",
    iconColor: "text-amber-500",
  },
  unauthorized: {
    icon: Lock,
    title: "Access restricted",
    description: "You don't have permission to view this content.",
    gradient: "from-twilight-800/10 to-twilight-800/5",
    iconColor: "text-twilight-800 dark:text-seagreen-400",
  },
  "coming-soon": {
    icon: Clock,
    title: "Coming soon",
    description: "This feature is under development and will be available soon.",
    gradient: "from-seagreen-500/10 to-mint-500/10",
    iconColor: "text-seagreen-500",
  },
  custom: {
    icon: Inbox,
    title: "Nothing here",
    description: "This section is empty.",
    gradient: "from-muted/50 to-muted/30",
    iconColor: "text-muted-foreground",
  },
}

const sizeClasses = {
  sm: {
    container: "py-8",
    iconWrapper: "h-12 w-12",
    icon: "h-6 w-6",
    title: "text-base",
    description: "text-sm",
    button: "text-sm",
  },
  md: {
    container: "py-12",
    iconWrapper: "h-16 w-16",
    icon: "h-8 w-8",
    title: "text-lg",
    description: "text-sm",
    button: "text-sm",
  },
  lg: {
    container: "py-16",
    iconWrapper: "h-20 w-20",
    icon: "h-10 w-10",
    title: "text-xl",
    description: "text-base",
    button: "text-base",
  },
}

export function EmptyState({
  preset = "no-data",
  icon: CustomIcon,
  title,
  description,
  action,
  secondaryAction,
  className,
  size = "md",
  children,
}: EmptyStateProps) {
  const config = presets[preset]
  const Icon = CustomIcon || config.icon
  const displayTitle = title || config.title
  const displayDescription = description || config.description
  const sizes = sizeClasses[size]

  return (
    <div
      className={cn(
        "flex flex-col items-center justify-center text-center",
        sizes.container,
        className
      )}
    >
      {/* Icon with gradient background */}
      <div className="relative mb-6">
        <div
          className={cn(
            "absolute inset-0 rounded-full blur-xl opacity-60",
            `bg-gradient-to-br ${config.gradient}`
          )}
        />
        <div
          className={cn(
            "relative flex items-center justify-center rounded-2xl",
            "bg-gradient-to-br border",
            config.gradient,
            "border-white/20 dark:border-white/10",
            sizes.iconWrapper
          )}
        >
          <Icon className={cn(sizes.icon, config.iconColor)} />
        </div>
      </div>

      {/* Text content */}
      <h3 className={cn("font-semibold mb-2", sizes.title)}>{displayTitle}</h3>
      <p
        className={cn(
          "text-muted-foreground max-w-sm mx-auto mb-6",
          sizes.description
        )}
      >
        {displayDescription}
      </p>

      {/* Actions */}
      {(action || secondaryAction) && (
        <div className="flex flex-col sm:flex-row gap-2">
          {action && (
            <Button
              onClick={action.onClick}
              variant={action.variant || "default"}
              className={cn(
                sizes.button,
                action.variant !== "outline" && action.variant !== "ghost"
                  ? "bg-gradient-to-r from-twilight-800 via-seagreen-500 to-mint-500 text-white hover:opacity-90"
                  : ""
              )}
            >
              {action.icon && <action.icon className="h-4 w-4 mr-2" />}
              {action.label}
            </Button>
          )}
          {secondaryAction && (
            <Button
              onClick={secondaryAction.onClick}
              variant={secondaryAction.variant || "outline"}
              className={sizes.button}
            >
              {secondaryAction.icon && (
                <secondaryAction.icon className="h-4 w-4 mr-2" />
              )}
              {secondaryAction.label}
            </Button>
          )}
        </div>
      )}

      {/* Custom content */}
      {children}
    </div>
  )
}

// Convenience components for common use cases
export function NoJobsEmptyState({
  onUpload,
}: {
  onUpload?: () => void
}) {
  return (
    <EmptyState
      preset="no-jobs"
      action={
        onUpload
          ? { label: "Upload Files", onClick: onUpload, icon: Upload }
          : undefined
      }
    />
  )
}

export function NoResultsEmptyState({
  onClearFilters,
}: {
  onClearFilters?: () => void
}) {
  return (
    <EmptyState
      preset="no-results"
      action={
        onClearFilters
          ? { label: "Clear Filters", onClick: onClearFilters, variant: "outline" }
          : undefined
      }
    />
  )
}

export function NoFilesEmptyState({
  onUpload,
}: {
  onUpload?: () => void
}) {
  return (
    <EmptyState
      preset="no-files"
      action={
        onUpload
          ? { label: "Upload Files", onClick: onUpload, icon: Upload }
          : undefined
      }
    />
  )
}

export function ErrorEmptyState({
  onRetry,
  message,
}: {
  onRetry?: () => void
  message?: string
}) {
  return (
    <EmptyState
      preset="error"
      description={message}
      action={
        onRetry
          ? { label: "Try Again", onClick: onRetry, variant: "outline" }
          : undefined
      }
    />
  )
}

export function ComingSoonEmptyState({
  feature,
}: {
  feature?: string
}) {
  return (
    <EmptyState
      preset="coming-soon"
      title={feature ? `${feature} - Coming Soon` : "Coming Soon"}
    />
  )
}
