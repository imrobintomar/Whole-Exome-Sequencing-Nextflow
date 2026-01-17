"use client"

import { cn } from "@/lib/utils"

interface SkipToContentProps {
  contentId?: string
  className?: string
}

export function SkipToContent({
  contentId = "main-content",
  className,
}: SkipToContentProps) {
  const handleClick = (e: React.MouseEvent) => {
    e.preventDefault()
    const content = document.getElementById(contentId)
    if (content) {
      content.tabIndex = -1
      content.focus()
      content.scrollIntoView({ behavior: "smooth" })
    }
  }

  return (
    <a
      href={`#${contentId}`}
      onClick={handleClick}
      className={cn(
        // Hidden by default, shown on focus
        "sr-only focus:not-sr-only",
        // Position and appearance when focused
        "focus:fixed focus:top-4 focus:left-4 focus:z-[9999]",
        "focus:px-4 focus:py-2 focus:rounded-lg",
        "focus:bg-gradient-to-r focus:from-twilight-800 focus:to-seagreen-500",
        "focus:text-white focus:font-medium focus:text-sm",
        "focus:shadow-lg focus:outline-none focus:ring-2 focus:ring-mint-500 focus:ring-offset-2",
        "transition-all duration-200",
        className
      )}
    >
      Skip to main content
    </a>
  )
}

// Skip links for multiple sections
interface SkipLink {
  id: string
  label: string
}

interface SkipLinksProps {
  links: SkipLink[]
  className?: string
}

export function SkipLinks({ links, className }: SkipLinksProps) {
  const handleClick = (e: React.MouseEvent, id: string) => {
    e.preventDefault()
    const element = document.getElementById(id)
    if (element) {
      element.tabIndex = -1
      element.focus()
      element.scrollIntoView({ behavior: "smooth" })
    }
  }

  return (
    <nav
      aria-label="Skip links"
      className={cn(
        "sr-only focus-within:not-sr-only",
        "focus-within:fixed focus-within:top-4 focus-within:left-4 focus-within:z-[9999]",
        "focus-within:flex focus-within:flex-col focus-within:gap-2",
        "focus-within:p-4 focus-within:rounded-xl",
        "focus-within:bg-white/95 focus-within:dark:bg-slate-900/95",
        "focus-within:backdrop-blur-xl focus-within:shadow-xl",
        "focus-within:border focus-within:border-seagreen-500/20",
        className
      )}
    >
      {links.map((link) => (
        <a
          key={link.id}
          href={`#${link.id}`}
          onClick={(e) => handleClick(e, link.id)}
          className={cn(
            "px-3 py-1.5 rounded-lg text-sm font-medium",
            "bg-muted hover:bg-seagreen-500/10",
            "text-foreground hover:text-seagreen-600 dark:hover:text-mint-400",
            "outline-none focus:ring-2 focus:ring-seagreen-500",
            "transition-colors"
          )}
        >
          Skip to {link.label}
        </a>
      ))}
    </nav>
  )
}
