"use client"

import { User, Menu } from "lucide-react"
import { Avatar, AvatarFallback } from "./ui/avatar"
import { Button } from "./ui/button"
import { User as UserType } from "@/lib/api"
import { useState } from "react"
import { Home, Upload, ListChecks, BarChart3, Filter } from "lucide-react"
import { cn } from "@/lib/utils"

interface DashboardHeaderProps {
  user: UserType
  onLogout: () => void
  currentView?: string
  onViewChange?: (view: string) => void
}

export function DashboardHeader({ user, onLogout, currentView, onViewChange }: DashboardHeaderProps) {
  const [mobileMenuOpen, setMobileMenuOpen] = useState(false)

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

  return (
    <>
      <header className="sticky top-0 z-10 flex h-16 items-center justify-between border-b bg-background px-4 sm:px-6">
        {/* Mobile menu button */}
        <Button
          variant="ghost"
          size="sm"
          className="lg:hidden"
          onClick={() => setMobileMenuOpen(!mobileMenuOpen)}
        >
          <Menu className="h-5 w-5" />
        </Button>

        <div className="flex-1 min-w-0 lg:block hidden">
          <h1 className="text-xl sm:text-2xl font-bold truncate">Welcome back, {user.displayName || "User"}</h1>
          <p className="text-xs sm:text-sm text-muted-foreground truncate">{user.email}</p>
        </div>

        {/* Mobile: Just show user avatar */}
        <div className="flex-1 lg:hidden" />

        <div className="flex items-center gap-2 sm:gap-4">
          <Avatar className="h-8 w-8 sm:h-10 sm:w-10">
            <AvatarFallback>{getInitials(user.displayName, user.email)}</AvatarFallback>
          </Avatar>
          <Button variant="outline" onClick={onLogout} size="sm" className="hidden sm:inline-flex">
            Logout
          </Button>
          <Button variant="outline" onClick={onLogout} size="sm" className="sm:hidden">
            <User className="h-4 w-4" />
          </Button>
        </div>
      </header>

      {/* Mobile Navigation Menu */}
      {mobileMenuOpen && (
        <div className="lg:hidden border-b bg-card p-4 space-y-2">
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
        </div>
      )}
    </>
  )
}
