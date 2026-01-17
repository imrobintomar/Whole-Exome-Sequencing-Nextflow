"use client"

import * as React from "react"
import { ChevronRight, Home } from "lucide-react"
import { cn } from "@/lib/utils"

interface BreadcrumbItem {
  label: string
  href?: string
  icon?: React.ElementType
  onClick?: () => void
}

interface BreadcrumbProps {
  items: BreadcrumbItem[]
  className?: string
}

export function Breadcrumb({ items, className }: BreadcrumbProps) {
  return (
    <nav
      aria-label="Breadcrumb"
      className={cn("flex items-center text-sm", className)}
    >
      <ol className="flex items-center gap-1">
        {items.map((item, index) => {
          const Icon = item.icon
          const isLast = index === items.length - 1
          const isFirst = index === 0

          return (
            <li key={index} className="flex items-center">
              {index > 0 && (
                <ChevronRight className="mx-1 h-4 w-4 text-muted-foreground/50" />
              )}

              {isLast ? (
                <span className={cn(
                  "flex items-center gap-1.5 font-medium",
                  "text-foreground"
                )}>
                  {Icon && <Icon className="h-4 w-4" />}
                  {item.label}
                </span>
              ) : (
                <button
                  onClick={item.onClick}
                  className={cn(
                    "flex items-center gap-1.5 rounded-md px-2 py-1",
                    "text-muted-foreground",
                    "hover:text-foreground hover:bg-muted/50",
                    "transition-colors duration-150"
                  )}
                >
                  {isFirst && !Icon ? (
                    <Home className="h-4 w-4" />
                  ) : Icon ? (
                    <Icon className="h-4 w-4" />
                  ) : null}
                  <span className={cn(isFirst && !item.label && "sr-only")}>
                    {item.label || "Home"}
                  </span>
                </button>
              )}
            </li>
          )
        })}
      </ol>
    </nav>
  )
}

// Helper to generate breadcrumbs from view
export function getBreadcrumbsForView(
  view: string,
  additionalItems?: BreadcrumbItem[],
  onNavigate?: (view: string) => void
): BreadcrumbItem[] {
  const viewLabels: Record<string, { label: string; icon?: React.ElementType }> = {
    overview: { label: "Dashboard" },
    upload: { label: "Upload" },
    jobs: { label: "Jobs" },
    "job-details": { label: "Job Details" },
    panels: { label: "Gene Panels" },
    analytics: { label: "Analytics" },
    acmg: { label: "ACMG Classification" },
    igv: { label: "IGV Browser" },
    variants: { label: "Variant Visualization" },
    support: { label: "Support" },
    phenotype: { label: "Phenotype Analysis" },
  }

  const items: BreadcrumbItem[] = [
    {
      label: "Home",
      onClick: () => onNavigate?.("overview")
    }
  ]

  // Add parent breadcrumb for nested views
  if (view === "job-details" || view === "acmg" || view === "igv" || view === "variants" || view === "phenotype") {
    items.push({
      label: "Jobs",
      onClick: () => onNavigate?.("jobs")
    })
  }

  // Add current view
  if (view !== "overview") {
    items.push({
      label: viewLabels[view]?.label || view,
    })
  }

  // Add any additional items (like job name)
  if (additionalItems) {
    items.push(...additionalItems)
  }

  return items
}
