"use client"

import { cn } from "@/lib/utils"
import {
  Home,
  Upload,
  ListChecks,
  BarChart3,
  User,
  Filter,
  Plus,
  Sparkles
} from "lucide-react"

interface MobileNavProps {
  currentView: string
  onViewChange: (view: string) => void
}

interface NavItem {
  id: string
  label: string
  icon: React.ElementType
  badge?: number
}

export function MobileNav({ currentView, onViewChange }: MobileNavProps) {
  const navItems: NavItem[] = [
    { id: "overview", label: "Home", icon: Home },
    { id: "jobs", label: "Jobs", icon: ListChecks },
    { id: "upload", label: "Upload", icon: Upload },
    { id: "panels", label: "Panels", icon: Filter },
    { id: "analytics", label: "Analytics", icon: BarChart3 },
  ]

  return (
    <nav className={cn(
      "lg:hidden fixed bottom-0 left-0 right-0 z-50",
      "bg-white/90 dark:bg-slate-900/90",
      "backdrop-blur-xl",
      "border-t border-white/20 dark:border-white/10",
      "shadow-[0_-4px_20px_rgba(0,0,0,0.1)]",
      "safe-area-bottom"
    )}>
      <div className="flex items-center justify-around h-16 px-2">
        {navItems.map((item) => {
          const Icon = item.icon
          const isActive = currentView === item.id
          const isUpload = item.id === "upload"

          // Special styling for upload button (center FAB style)
          if (isUpload) {
            return (
              <button
                key={item.id}
                onClick={() => onViewChange(item.id)}
                className={cn(
                  "relative -mt-6 flex flex-col items-center justify-center",
                  "transition-all duration-300 ease-out"
                )}
              >
                <div className={cn(
                  "flex h-14 w-14 items-center justify-center rounded-2xl",
                  "bg-gradient-to-br from-purple-500 to-cyan-500",
                  "shadow-lg shadow-purple-500/30",
                  "transition-all duration-300",
                  isActive && "scale-110 shadow-xl shadow-purple-500/40"
                )}>
                  <Plus className="h-6 w-6 text-white" />
                </div>
                <span className={cn(
                  "text-[10px] font-medium mt-1",
                  isActive
                    ? "text-purple-600 dark:text-purple-400"
                    : "text-muted-foreground"
                )}>
                  {item.label}
                </span>
              </button>
            )
          }

          return (
            <button
              key={item.id}
              onClick={() => onViewChange(item.id)}
              className={cn(
                "relative flex flex-col items-center justify-center",
                "min-w-[60px] py-2 px-3 rounded-xl",
                "transition-all duration-300 ease-out",
                isActive && "bg-purple-500/10 dark:bg-purple-500/20"
              )}
            >
              {/* Active indicator dot */}
              {isActive && (
                <span className={cn(
                  "absolute top-0 left-1/2 -translate-x-1/2",
                  "h-1 w-6 rounded-full",
                  "bg-gradient-to-r from-purple-500 to-cyan-500",
                  "animate-scale-in"
                )} />
              )}

              {/* Icon container */}
              <div className={cn(
                "relative flex items-center justify-center",
                "h-7 w-7 rounded-lg",
                "transition-all duration-300",
                isActive && "scale-110"
              )}>
                <Icon className={cn(
                  "h-5 w-5 transition-colors duration-300",
                  isActive
                    ? "text-purple-600 dark:text-purple-400"
                    : "text-muted-foreground"
                )} />

                {/* Badge */}
                {item.badge && item.badge > 0 && (
                  <span className={cn(
                    "absolute -top-1 -right-1",
                    "flex h-4 w-4 items-center justify-center",
                    "rounded-full text-[10px] font-bold",
                    "bg-gradient-to-r from-purple-500 to-cyan-500 text-white",
                    "animate-pulse"
                  )}>
                    {item.badge > 9 ? "9+" : item.badge}
                  </span>
                )}
              </div>

              {/* Label */}
              <span className={cn(
                "text-[10px] font-medium mt-0.5",
                "transition-colors duration-300",
                isActive
                  ? "text-purple-600 dark:text-purple-400"
                  : "text-muted-foreground"
              )}>
                {item.label}
              </span>
            </button>
          )
        })}
      </div>
    </nav>
  )
}

// Floating Action Button variant for quick upload
export function MobileUploadFAB({ onClick }: { onClick: () => void }) {
  return (
    <button
      onClick={onClick}
      className={cn(
        "lg:hidden fixed bottom-20 right-4 z-40",
        "flex h-14 w-14 items-center justify-center",
        "rounded-2xl",
        "bg-gradient-to-br from-purple-500 to-cyan-500",
        "shadow-lg shadow-purple-500/30",
        "hover:shadow-xl hover:shadow-purple-500/40",
        "hover:scale-105",
        "active:scale-95",
        "transition-all duration-300 ease-out"
      )}
    >
      <div className="relative">
        <Plus className="h-6 w-6 text-white" />
        <Sparkles className="absolute -top-1 -right-1 h-3 w-3 text-white animate-pulse" />
      </div>
    </button>
  )
}

// Bottom sheet trigger for mobile actions
export function MobileActionSheet({
  isOpen,
  onClose,
  onViewChange
}: {
  isOpen: boolean
  onClose: () => void
  onViewChange: (view: string) => void
}) {
  if (!isOpen) return null

  const actions = [
    { id: "upload", label: "Upload New Sample", icon: Upload, description: "Start a new WES analysis" },
    { id: "panels", label: "Gene Panels", icon: Filter, description: "Filter by gene panels" },
    { id: "phenotype", label: "Phenotype Analysis", icon: Sparkles, description: "HPO-driven prioritization" },
  ]

  return (
    <>
      {/* Backdrop */}
      <div
        className={cn(
          "lg:hidden fixed inset-0 z-50",
          "bg-black/50 backdrop-blur-sm",
          "animate-fade-in"
        )}
        onClick={onClose}
      />

      {/* Sheet */}
      <div className={cn(
        "lg:hidden fixed bottom-0 left-0 right-0 z-50",
        "bg-white dark:bg-slate-900",
        "rounded-t-3xl",
        "shadow-2xl",
        "animate-slide-in-bottom",
        "safe-area-bottom"
      )}>
        {/* Handle */}
        <div className="flex justify-center pt-3 pb-2">
          <div className="h-1.5 w-12 rounded-full bg-muted" />
        </div>

        {/* Header */}
        <div className="px-6 pb-4">
          <h3 className="text-lg font-semibold">Quick Actions</h3>
          <p className="text-sm text-muted-foreground">Choose an action to get started</p>
        </div>

        {/* Actions */}
        <div className="px-4 pb-6 space-y-2">
          {actions.map((action) => {
            const Icon = action.icon
            return (
              <button
                key={action.id}
                onClick={() => {
                  onViewChange(action.id)
                  onClose()
                }}
                className={cn(
                  "flex items-center gap-4 w-full p-4 rounded-2xl",
                  "bg-muted/30 hover:bg-muted/50",
                  "transition-all duration-200",
                  "active:scale-[0.98]"
                )}
              >
                <div className={cn(
                  "flex h-12 w-12 items-center justify-center rounded-xl",
                  "bg-gradient-to-br from-purple-500/10 to-cyan-500/10",
                  "border border-purple-500/20"
                )}>
                  <Icon className="h-6 w-6 text-purple-500" />
                </div>
                <div className="flex-1 text-left">
                  <p className="font-medium">{action.label}</p>
                  <p className="text-sm text-muted-foreground">{action.description}</p>
                </div>
              </button>
            )
          })}
        </div>

        {/* Cancel button */}
        <div className="px-4 pb-6">
          <button
            onClick={onClose}
            className={cn(
              "w-full py-3 rounded-xl",
              "bg-muted/50 hover:bg-muted",
              "font-medium",
              "transition-colors duration-200"
            )}
          >
            Cancel
          </button>
        </div>
      </div>
    </>
  )
}
