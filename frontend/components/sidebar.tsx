"use client"

import { useState } from "react"
import {
  Home,
  Upload,
  ListChecks,
  BarChart3,
  Dna,
  Filter,
  MessageCircle,
  ChevronLeft,
  ChevronRight,
  Settings,
  HelpCircle,
  Sparkles
} from "lucide-react"
import { cn } from "@/lib/utils"
import { ThemeToggle } from "./theme-toggle"
import { Badge } from "./ui/badge"

interface SidebarProps {
  currentView: string
  onViewChange: (view: string) => void
}

export function Sidebar({ currentView, onViewChange }: SidebarProps) {
  const [collapsed, setCollapsed] = useState(false)

  const menuItems = [
    { id: "overview", label: "Overview", icon: Home, badge: null },
    { id: "upload", label: "Upload", icon: Upload, badge: null },
    { id: "jobs", label: "Jobs", icon: ListChecks, badge: null },
    { id: "panels", label: "Gene Panels", icon: Filter, badge: null },
    { id: "analytics", label: "Analytics", icon: BarChart3, badge: "New" },
    { id: "support", label: "Support", icon: MessageCircle, badge: null },
  ]

  const bottomItems = [
    { id: "help", label: "Help", icon: HelpCircle },
  ]

  return (
    <div
      className={cn(
        "hidden lg:flex h-screen flex-col",
        "bg-white/80 dark:bg-slate-900/80",
        "backdrop-blur-xl",
        "border-r border-white/20 dark:border-white/10",
        "shadow-[4px_0_24px_rgba(0,0,0,0.05)] dark:shadow-[4px_0_24px_rgba(0,0,0,0.2)]",
        "transition-all duration-300 ease-out",
        collapsed ? "w-20" : "w-64"
      )}
    >
      {/* Logo/Brand */}
      <div className={cn(
        "flex h-16 items-center gap-3 border-b border-seagreen-500/10 px-4",
        collapsed ? "justify-center" : "px-6"
      )}>
        <div className="relative">
          <div className="absolute inset-0 bg-gradient-to-br from-twilight-800 to-mint-500 rounded-xl blur-lg opacity-50" />
          <div className="relative flex h-10 w-10 items-center justify-center rounded-xl bg-gradient-to-br from-twilight-800 to-seagreen-500 shadow-lg">
            <Dna className="h-5 w-5 text-white" />
          </div>
        </div>
        {!collapsed && (
          <div className="flex flex-col animate-fade-in">
            <span className="text-lg font-bold bg-gradient-to-r from-twilight-800 via-seagreen-500 to-mint-500 dark:from-seagreen-400 dark:via-mint-400 dark:to-mint-300 bg-clip-text text-transparent">
              ATGCFlow
            </span>
            <span className="text-[10px] text-muted-foreground font-medium tracking-wide uppercase">
              Whole Exome Sequencing
            </span>
          </div>
        )}
      </div>

      {/* Navigation */}
      <nav className="flex-1 space-y-1 p-3 overflow-y-auto scrollbar-thin">
        <div className="space-y-1">
          {!collapsed && (
            <p className="px-3 py-2 text-xs font-semibold text-muted-foreground uppercase tracking-wider">
              Main Menu
            </p>
          )}
          {menuItems.map((item) => {
            const Icon = item.icon
            const isActive = currentView === item.id

            return (
              <button
                key={item.id}
                onClick={() => onViewChange(item.id)}
                className={cn(
                  "group relative flex w-full items-center gap-3 rounded-xl px-3 py-2.5 text-sm font-medium",
                  "transition-all duration-200 ease-out",
                  isActive
                    ? [
                        "bg-gradient-to-r from-twilight-800/10 via-seagreen-500/10 to-mint-500/10",
                        "text-twilight-800 dark:text-mint-400",
                        "shadow-sm",
                      ]
                    : [
                        "text-muted-foreground",
                        "hover:bg-muted/50 hover:text-foreground",
                      ],
                  collapsed && "justify-center px-0"
                )}
              >
                {/* Active indicator */}
                {isActive && (
                  <div className="absolute left-0 top-1/2 -translate-y-1/2 h-8 w-1 rounded-r-full bg-gradient-to-b from-twilight-800 via-seagreen-500 to-mint-500 animate-fade-in" />
                )}

                <div className={cn(
                  "flex h-9 w-9 items-center justify-center rounded-lg transition-all duration-200",
                  isActive
                    ? "bg-gradient-to-br from-twilight-800 to-seagreen-500 text-white shadow-lg shadow-seagreen-500/25"
                    : "bg-muted/50 group-hover:bg-muted"
                )}>
                  <Icon className="h-5 w-5" />
                </div>

                {!collapsed && (
                  <span className="flex-1 text-left animate-fade-in">{item.label}</span>
                )}

                {!collapsed && item.badge && (
                  <Badge
                    variant="secondary"
                    className="bg-gradient-to-r from-seagreen-500 to-mint-500 text-white text-[10px] px-1.5 py-0 border-0"
                  >
                    {item.badge}
                  </Badge>
                )}
              </button>
            )
          })}
        </div>
      </nav>

      {/* Upgrade Card */}
      {!collapsed && (
        <div className="mx-3 mb-3">
          <div className="relative overflow-hidden rounded-xl bg-gradient-to-br from-twilight-800/10 via-seagreen-500/5 to-mint-500/10 p-4 border border-seagreen-500/20">
            <div className="absolute -right-4 -top-4 h-20 w-20 rounded-full bg-gradient-to-br from-seagreen-500/20 to-mint-500/20 blur-2xl" />
            <div className="relative">
              <div className="flex items-center gap-2 mb-2">
                <Sparkles className="h-4 w-4 text-mint-500" />
                <span className="text-sm font-semibold">Pro Features</span>
              </div>
              <p className="text-xs text-muted-foreground mb-3">
                Unlock advanced analysis tools and priority support.
              </p>
              <button className="w-full py-2 px-3 rounded-lg bg-gradient-to-r from-twilight-800 via-seagreen-500 to-mint-500 text-white text-sm font-medium hover:opacity-90 transition-opacity btn-press">
                Upgrade Now
              </button>
            </div>
          </div>
        </div>
      )}

      {/* Bottom Section */}
      <div className="border-t border-white/10 p-3 space-y-2">
        {/* Help Button */}
        {bottomItems.map((item) => {
          const Icon = item.icon
          return (
            <button
              key={item.id}
              className={cn(
                "flex w-full items-center gap-3 rounded-xl px-3 py-2 text-sm font-medium",
                "text-muted-foreground hover:bg-muted/50 hover:text-foreground",
                "transition-all duration-200",
                collapsed && "justify-center px-0"
              )}
            >
              <div className="flex h-9 w-9 items-center justify-center rounded-lg bg-muted/50">
                <Icon className="h-5 w-5" />
              </div>
              {!collapsed && <span>{item.label}</span>}
            </button>
          )
        })}

        {/* Theme Toggle */}
        <div className={cn(
          "flex items-center rounded-xl px-3 py-2",
          collapsed ? "justify-center" : "justify-between"
        )}>
          {!collapsed && (
            <span className="text-xs text-muted-foreground font-medium">Theme</span>
          )}
          <ThemeToggle />
        </div>

        {/* Collapse Button */}
        <button
          onClick={() => setCollapsed(!collapsed)}
          className={cn(
            "flex w-full items-center gap-3 rounded-xl px-3 py-2 text-sm font-medium",
            "text-muted-foreground hover:bg-muted/50 hover:text-foreground",
            "transition-all duration-200",
            collapsed && "justify-center px-0"
          )}
        >
          <div className="flex h-9 w-9 items-center justify-center rounded-lg bg-muted/50">
            {collapsed ? (
              <ChevronRight className="h-5 w-5" />
            ) : (
              <ChevronLeft className="h-5 w-5" />
            )}
          </div>
          {!collapsed && <span>Collapse</span>}
        </button>
      </div>
    </div>
  )
}
