"use client"

import { Home, Upload, ListChecks, BarChart3, Dna } from "lucide-react"
import Link from "next/link"
import { usePathname } from "next/navigation"
import { cn } from "@/lib/utils"
import { ThemeToggle } from "./theme-toggle"
import { Separator } from "./ui/separator"

interface SidebarProps {
  currentView: string
  onViewChange: (view: string) => void
}

export function Sidebar({ currentView, onViewChange }: SidebarProps) {
  const menuItems = [
    { id: "overview", label: "Overview", icon: Home },
    { id: "upload", label: "Upload", icon: Upload },
    { id: "jobs", label: "Jobs", icon: ListChecks },
    { id: "analytics", label: "Analytics", icon: BarChart3 },
  ]

  return (
    <div className="flex h-screen w-64 flex-col border-r bg-card">
      {/* Logo/Brand */}
      <div className="flex h-16 items-center gap-2 border-b px-6">
        <Dna className="h-6 w-6 text-primary" />
        <div className="flex flex-col">
          <span className="text-lg font-bold">WES Pipeline</span>
          <span className="text-xs text-muted-foreground">Exome Sequencing</span>
        </div>
      </div>

      {/* Navigation */}
      <nav className="flex-1 space-y-1 p-4">
        {menuItems.map((item) => {
          const Icon = item.icon
          const isActive = currentView === item.id

          return (
            <button
              key={item.id}
              onClick={() => onViewChange(item.id)}
              className={cn(
                "flex w-full items-center gap-3 rounded-lg px-3 py-2 text-sm font-medium transition-colors",
                isActive
                  ? "bg-primary text-primary-foreground"
                  : "text-muted-foreground hover:bg-accent hover:text-accent-foreground"
              )}
            >
              <Icon className="h-5 w-5" />
              {item.label}
            </button>
          )
        })}
      </nav>

      {/* Footer with Theme Toggle */}
      <div className="border-t p-4">
        <div className="flex items-center justify-between">
          <span className="text-xs text-muted-foreground">Theme</span>
          <ThemeToggle />
        </div>
      </div>
    </div>
  )
}
