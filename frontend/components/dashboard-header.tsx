"use client"

import { useState } from "react"
import {
  User,
  Menu,
  Search,
  Bell,
  LogOut,
  Settings,
  ChevronDown,
  Home,
  Upload,
  ListChecks,
  BarChart3,
  Filter,
  X,
  Command
} from "lucide-react"
import { Avatar, AvatarFallback, AvatarImage } from "./ui/avatar"
import { Button } from "./ui/button"
import { User as UserType } from "@/lib/api"
import { cn } from "@/lib/utils"
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuLabel,
  DropdownMenuSeparator,
  DropdownMenuTrigger,
} from "@/components/ui/dropdown-menu""

interface DashboardHeaderProps {
  user: UserType
  onLogout: () => void
  currentView?: string
  onViewChange?: (view: string) => void
}

export function DashboardHeader({ user, onLogout, currentView, onViewChange }: DashboardHeaderProps) {
  const [mobileMenuOpen, setMobileMenuOpen] = useState(false)
  const [searchOpen, setSearchOpen] = useState(false)

  const menuItems = [
    { id: "overview", label: "Overview", icon: Home },
    { id: "upload", label: "Upload", icon: Upload },
    { id: "jobs", label: "Jobs", icon: ListChecks },
    { id: "panels", label: "Gene Panels", icon: Filter },
    { id: "analytics", label: "Analytics", icon: BarChart3 },
  ]

  const getInitials = (name: string | undefined, email: string) => {
    if (name) {
      return name
        .split(" ")
        .map((n) => n[0])
        .join("")
        .toUpperCase()
        .slice(0, 2)
    }
    if (email) {
      return email.charAt(0).toUpperCase()
    }
    return "U"
  }

  const getViewTitle = (view: string | undefined) => {
    const titles: Record<string, string> = {
      overview: "Dashboard",
      upload: "Upload Files",
      jobs: "Analysis Jobs",
      panels: "Gene Panels",
      analytics: "Analytics",
      "job-details": "Job Details",
      acmg: "ACMG Classification",
      igv: "Genome Browser",
      variants: "Variant Visualization",
      "gene-panel": "Gene Panel Analysis",
      support: "Support",
      phenotype: "Phenotype Analysis",
    }
    return titles[view || "overview"] || "Dashboard"
  }

  return (
    <>
      <header className={cn(
        "sticky top-0 z-40",
        "flex h-16 items-center justify-between",
        "px-4 sm:px-6",
        "bg-white/80 dark:bg-slate-900/80",
        "backdrop-blur-xl",
        "border-b border-white/20 dark:border-white/10",
        "shadow-sm"
      )}>
        {/* Left Section - Mobile Menu & Title */}
        <div className="flex items-center gap-4">
          {/* Mobile menu button */}
          <Button
            variant="ghost"
            size="icon"
            className="lg:hidden h-9 w-9 rounded-xl"
            onClick={() => setMobileMenuOpen(!mobileMenuOpen)}
          >
            {mobileMenuOpen ? (
              <X className="h-5 w-5" />
            ) : (
              <Menu className="h-5 w-5" />
            )}
          </Button>

          {/* Page Title */}
          <div className="hidden sm:block">
            <h1 className="text-xl font-semibold">{getViewTitle(currentView)}</h1>
          </div>
        </div>

        {/* Center Section - Search (Desktop) */}
        <div className="hidden md:flex flex-1 max-w-md mx-8">
          <div className="relative w-full">
            <Search className="absolute left-3 top-1/2 -translate-y-1/2 h-4 w-4 text-muted-foreground" />
            <input
              type="text"
              placeholder="Search jobs, samples..."
              className={cn(
                "w-full h-10 pl-10 pr-12 rounded-xl",
                "bg-muted/50 border border-transparent",
                "text-sm placeholder:text-muted-foreground",
                "focus:outline-none focus:ring-2 focus:ring-primary/20 focus:border-primary/50",
                "transition-all duration-200"
              )}
            />
            <div className="absolute right-3 top-1/2 -translate-y-1/2 flex items-center gap-1 text-xs text-muted-foreground">
              <kbd className="px-1.5 py-0.5 rounded bg-background border text-[10px] font-medium">
                <Command className="h-3 w-3 inline" />
              </kbd>
              <kbd className="px-1.5 py-0.5 rounded bg-background border text-[10px] font-medium">K</kbd>
            </div>
          </div>
        </div>

        {/* Right Section - Actions & User */}
        <div className="flex items-center gap-2 sm:gap-3">
          {/* Mobile Search Button */}
          <Button
            variant="ghost"
            size="icon"
            className="md:hidden h-9 w-9 rounded-xl"
            onClick={() => setSearchOpen(!searchOpen)}
          >
            <Search className="h-5 w-5" />
          </Button>

          {/* Notifications */}
          <Button
            variant="ghost"
            size="icon"
            className="relative h-9 w-9 rounded-xl"
          >
            <Bell className="h-5 w-5" />
            {/* Notification badge */}
            <span className="absolute top-1.5 right-1.5 h-2 w-2 rounded-full bg-gradient-to-r from-purple-500 to-cyan-500 animate-pulse" />
          </Button>

          {/* User Menu */}
          <DropdownMenu>
            <DropdownMenuTrigger asChild>
              <button className={cn(
                "flex items-center gap-2 px-2 py-1.5 rounded-xl",
                "hover:bg-muted/50 transition-colors duration-200",
                "focus:outline-none focus:ring-2 focus:ring-primary/20"
              )}>
                <Avatar className="h-8 w-8 ring-2 ring-purple-500/20">
                  <AvatarImage src={user.photoURL || undefined} />
                  <AvatarFallback className="bg-gradient-to-br from-purple-500 to-cyan-500 text-white text-sm font-medium">
                    {getInitials(user.displayName, user.email)}
                  </AvatarFallback>
                </Avatar>
                <div className="hidden sm:block text-left">
                  <p className="text-sm font-medium leading-none">{user.displayName || "User"}</p>
                  <p className="text-xs text-muted-foreground truncate max-w-[120px]">{user.email}</p>
                </div>
                <ChevronDown className="h-4 w-4 text-muted-foreground hidden sm:block" />
              </button>
            </DropdownMenuTrigger>
            <DropdownMenuContent align="end" className="w-56 rounded-xl">
              <DropdownMenuLabel className="font-normal">
                <div className="flex flex-col space-y-1">
                  <p className="text-sm font-medium">{user.displayName || "User"}</p>
                  <p className="text-xs text-muted-foreground">{user.email}</p>
                </div>
              </DropdownMenuLabel>
              <DropdownMenuSeparator />
              <DropdownMenuItem className="cursor-pointer rounded-lg">
                <User className="mr-2 h-4 w-4" />
                Profile
              </DropdownMenuItem>
              <DropdownMenuItem className="cursor-pointer rounded-lg">
                <Settings className="mr-2 h-4 w-4" />
                Settings
              </DropdownMenuItem>
              <DropdownMenuSeparator />
              <DropdownMenuItem
                className="cursor-pointer rounded-lg text-destructive focus:text-destructive"
                onClick={onLogout}
              >
                <LogOut className="mr-2 h-4 w-4" />
                Log out
              </DropdownMenuItem>
            </DropdownMenuContent>
          </DropdownMenu>
        </div>
      </header>

      {/* Mobile Search Bar */}
      {searchOpen && (
        <div className="md:hidden p-4 bg-background/80 backdrop-blur-xl border-b animate-fade-in">
          <div className="relative">
            <Search className="absolute left-3 top-1/2 -translate-y-1/2 h-4 w-4 text-muted-foreground" />
            <input
              type="text"
              placeholder="Search jobs, samples..."
              autoFocus
              className={cn(
                "w-full h-10 pl-10 pr-4 rounded-xl",
                "bg-muted/50 border border-transparent",
                "text-sm placeholder:text-muted-foreground",
                "focus:outline-none focus:ring-2 focus:ring-primary/20 focus:border-primary/50"
              )}
            />
          </div>
        </div>
      )}

      {/* Mobile Navigation Menu */}
      {mobileMenuOpen && (
        <div className={cn(
          "lg:hidden fixed inset-x-0 top-16 bottom-0 z-30",
          "bg-background/95 backdrop-blur-xl",
          "animate-fade-in"
        )}>
          <nav className="p-4 space-y-2">
            {menuItems.map((item) => {
              const Icon = item.icon
              const isActive = currentView === item.id

              return (
                <button
                  key={item.id}
                  onClick={() => {
                    onViewChange?.(item.id)
                    setMobileMenuOpen(false)
                  }}
                  className={cn(
                    "flex w-full items-center gap-3 rounded-xl px-4 py-3 text-sm font-medium",
                    "transition-all duration-200",
                    isActive
                      ? [
                          "bg-gradient-to-r from-purple-500/10 to-cyan-500/10",
                          "text-purple-600 dark:text-purple-400",
                        ]
                      : [
                          "text-muted-foreground",
                          "hover:bg-muted/50 hover:text-foreground",
                        ]
                  )}
                >
                  <div className={cn(
                    "flex h-10 w-10 items-center justify-center rounded-lg",
                    isActive
                      ? "bg-gradient-to-br from-purple-500 to-cyan-500 text-white"
                      : "bg-muted"
                  )}>
                    <Icon className="h-5 w-5" />
                  </div>
                  <span>{item.label}</span>
                </button>
              )
            })}
          </nav>

          {/* User info at bottom */}
          <div className="absolute bottom-0 left-0 right-0 p-4 border-t bg-background/50 backdrop-blur-xl">
            <div className="flex items-center justify-between">
              <div className="flex items-center gap-3">
                <Avatar className="h-10 w-10">
                  <AvatarImage src={user.photoURL || undefined} />
                  <AvatarFallback className="bg-gradient-to-br from-purple-500 to-cyan-500 text-white">
                    {getInitials(user.displayName, user.email)}
                  </AvatarFallback>
                </Avatar>
                <div>
                  <p className="text-sm font-medium">{user.displayName || "User"}</p>
                  <p className="text-xs text-muted-foreground">{user.email}</p>
                </div>
              </div>
              <Button variant="outline" size="sm" onClick={onLogout} className="rounded-lg">
                <LogOut className="h-4 w-4 mr-2" />
                Logout
              </Button>
            </div>
          </div>
        </div>
      )}
    </>
  )
}
