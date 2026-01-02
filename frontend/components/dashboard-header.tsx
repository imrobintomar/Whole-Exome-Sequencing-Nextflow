"use client"

import { User } from "lucide-react"
import { Avatar, AvatarFallback } from "./ui/avatar"
import { Button } from "./ui/button"

interface DashboardHeaderProps {
  user: {
    displayName: string | null
    email: string | null
  }
  onLogout: () => void
}

export function DashboardHeader({ user, onLogout }: DashboardHeaderProps) {
  const getInitials = (name: string | null, email: string | null) => {
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
    <header className="sticky top-0 z-10 flex h-16 items-center justify-between border-b bg-background px-6">
      <div className="flex-1">
        <h1 className="text-2xl font-bold">Welcome back, {user.displayName || "User"}</h1>
        <p className="text-sm text-muted-foreground">{user.email}</p>
      </div>

      <div className="flex items-center gap-4">
        <Avatar>
          <AvatarFallback>{getInitials(user.displayName, user.email)}</AvatarFallback>
        </Avatar>
        <Button variant="outline" onClick={onLogout}>
          Logout
        </Button>
      </div>
    </header>
  )
}
