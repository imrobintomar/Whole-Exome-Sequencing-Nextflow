"use client"

import { useState, useEffect, useCallback, useRef } from "react"
import { cn } from "@/lib/utils"
import {
  Home,
  Upload,
  ListChecks,
  BarChart3,
  Filter,
  MessageCircle,
  Search,
  Command,
  ArrowRight,
  FileText,
  Settings,
  HelpCircle,
  Dna,
  Sparkles
} from "lucide-react"

interface CommandPaletteProps {
  isOpen: boolean
  onClose: () => void
  onViewChange: (view: string) => void
}

interface CommandItem {
  id: string
  label: string
  description: string
  icon: React.ElementType
  category: "navigation" | "actions" | "help"
  shortcut?: string
  action?: () => void
}

export function CommandPalette({ isOpen, onClose, onViewChange }: CommandPaletteProps) {
  const [search, setSearch] = useState("")
  const [selectedIndex, setSelectedIndex] = useState(0)
  const inputRef = useRef<HTMLInputElement>(null)
  const listRef = useRef<HTMLDivElement>(null)

  const commands: CommandItem[] = [
    // Navigation
    { id: "overview", label: "Go to Dashboard", description: "View dashboard overview", icon: Home, category: "navigation", shortcut: "G D" },
    { id: "upload", label: "Upload Sample", description: "Start a new WES analysis", icon: Upload, category: "navigation", shortcut: "G U" },
    { id: "jobs", label: "View Jobs", description: "See all analysis jobs", icon: ListChecks, category: "navigation", shortcut: "G J" },
    { id: "panels", label: "Gene Panels", description: "Filter by gene panels", icon: Filter, category: "navigation", shortcut: "G P" },
    { id: "analytics", label: "Analytics", description: "View job statistics", icon: BarChart3, category: "navigation", shortcut: "G A" },
    { id: "support", label: "Support", description: "Get help and support", icon: MessageCircle, category: "navigation", shortcut: "G S" },
    // Actions
    { id: "new-analysis", label: "New Analysis", description: "Start a new WES analysis", icon: Sparkles, category: "actions", shortcut: "N" },
    { id: "search-jobs", label: "Search Jobs", description: "Find a specific job", icon: Search, category: "actions" },
    { id: "export-data", label: "Export Data", description: "Export analysis results", icon: FileText, category: "actions" },
    // Help
    { id: "settings", label: "Settings", description: "Configure your preferences", icon: Settings, category: "help" },
    { id: "help", label: "Help & Documentation", description: "View documentation", icon: HelpCircle, category: "help", shortcut: "?" },
  ]

  const filteredCommands = commands.filter(cmd =>
    cmd.label.toLowerCase().includes(search.toLowerCase()) ||
    cmd.description.toLowerCase().includes(search.toLowerCase())
  )

  const groupedCommands = {
    navigation: filteredCommands.filter(c => c.category === "navigation"),
    actions: filteredCommands.filter(c => c.category === "actions"),
    help: filteredCommands.filter(c => c.category === "help"),
  }

  const flatFilteredCommands = filteredCommands

  const handleSelect = useCallback((command: CommandItem) => {
    if (command.action) {
      command.action()
    } else {
      onViewChange(command.id)
    }
    onClose()
    setSearch("")
    setSelectedIndex(0)
  }, [onViewChange, onClose])

  // Keyboard navigation
  useEffect(() => {
    if (!isOpen) return

    const handleKeyDown = (e: KeyboardEvent) => {
      switch (e.key) {
        case "ArrowDown":
          e.preventDefault()
          setSelectedIndex(i => Math.min(i + 1, flatFilteredCommands.length - 1))
          break
        case "ArrowUp":
          e.preventDefault()
          setSelectedIndex(i => Math.max(i - 1, 0))
          break
        case "Enter":
          e.preventDefault()
          if (flatFilteredCommands[selectedIndex]) {
            handleSelect(flatFilteredCommands[selectedIndex])
          }
          break
        case "Escape":
          e.preventDefault()
          onClose()
          break
      }
    }

    document.addEventListener("keydown", handleKeyDown)
    return () => document.removeEventListener("keydown", handleKeyDown)
  }, [isOpen, selectedIndex, flatFilteredCommands, handleSelect, onClose])

  // Focus input when opened
  useEffect(() => {
    if (isOpen && inputRef.current) {
      inputRef.current.focus()
    }
  }, [isOpen])

  // Reset on close
  useEffect(() => {
    if (!isOpen) {
      setSearch("")
      setSelectedIndex(0)
    }
  }, [isOpen])

  // Scroll selected item into view
  useEffect(() => {
    if (listRef.current) {
      const selectedElement = listRef.current.querySelector(`[data-index="${selectedIndex}"]`)
      if (selectedElement) {
        selectedElement.scrollIntoView({ block: "nearest" })
      }
    }
  }, [selectedIndex])

  if (!isOpen) return null

  return (
    <>
      {/* Backdrop */}
      <div
        className="fixed inset-0 z-50 bg-black/50 backdrop-blur-sm animate-fade-in"
        onClick={onClose}
      />

      {/* Command Palette */}
      <div className={cn(
        "fixed left-1/2 top-[20%] z-50 w-full max-w-lg -translate-x-1/2",
        "animate-scale-in"
      )}>
        <div className={cn(
          "overflow-hidden rounded-2xl",
          "bg-white/95 dark:bg-twilight-900/95",
          "backdrop-blur-xl",
          "border border-seagreen-500/20",
          "shadow-2xl shadow-twilight-800/20"
        )}>
          {/* Search Input */}
          <div className="flex items-center gap-3 border-b border-seagreen-500/10 px-4 py-3">
            <Search className="h-5 w-5 text-seagreen-500" />
            <input
              ref={inputRef}
              type="text"
              placeholder="Search commands..."
              value={search}
              onChange={(e) => {
                setSearch(e.target.value)
                setSelectedIndex(0)
              }}
              className={cn(
                "flex-1 bg-transparent text-base outline-none",
                "placeholder:text-muted-foreground"
              )}
            />
            <kbd className={cn(
              "hidden sm:flex items-center gap-1 px-2 py-1 rounded-lg",
              "bg-muted/50 text-xs text-muted-foreground font-mono"
            )}>
              ESC
            </kbd>
          </div>

          {/* Commands List */}
          <div ref={listRef} className="max-h-[60vh] overflow-y-auto p-2 scrollbar-thin">
            {flatFilteredCommands.length === 0 ? (
              <div className="flex flex-col items-center justify-center py-12 text-center">
                <Dna className="h-12 w-12 text-muted-foreground/50 mb-3" />
                <p className="text-sm text-muted-foreground">No commands found</p>
                <p className="text-xs text-muted-foreground/70 mt-1">Try a different search term</p>
              </div>
            ) : (
              <>
                {groupedCommands.navigation.length > 0 && (
                  <div className="mb-2">
                    <p className="px-3 py-2 text-xs font-semibold text-muted-foreground uppercase tracking-wider">
                      Navigation
                    </p>
                    {groupedCommands.navigation.map((command, index) => {
                      const globalIndex = flatFilteredCommands.indexOf(command)
                      return (
                        <CommandItem
                          key={command.id}
                          command={command}
                          isSelected={selectedIndex === globalIndex}
                          dataIndex={globalIndex}
                          onClick={() => handleSelect(command)}
                          onMouseEnter={() => setSelectedIndex(globalIndex)}
                        />
                      )
                    })}
                  </div>
                )}

                {groupedCommands.actions.length > 0 && (
                  <div className="mb-2">
                    <p className="px-3 py-2 text-xs font-semibold text-muted-foreground uppercase tracking-wider">
                      Actions
                    </p>
                    {groupedCommands.actions.map((command) => {
                      const globalIndex = flatFilteredCommands.indexOf(command)
                      return (
                        <CommandItem
                          key={command.id}
                          command={command}
                          isSelected={selectedIndex === globalIndex}
                          dataIndex={globalIndex}
                          onClick={() => handleSelect(command)}
                          onMouseEnter={() => setSelectedIndex(globalIndex)}
                        />
                      )
                    })}
                  </div>
                )}

                {groupedCommands.help.length > 0 && (
                  <div className="mb-2">
                    <p className="px-3 py-2 text-xs font-semibold text-muted-foreground uppercase tracking-wider">
                      Help
                    </p>
                    {groupedCommands.help.map((command) => {
                      const globalIndex = flatFilteredCommands.indexOf(command)
                      return (
                        <CommandItem
                          key={command.id}
                          command={command}
                          isSelected={selectedIndex === globalIndex}
                          dataIndex={globalIndex}
                          onClick={() => handleSelect(command)}
                          onMouseEnter={() => setSelectedIndex(globalIndex)}
                        />
                      )
                    })}
                  </div>
                )}
              </>
            )}
          </div>

          {/* Footer */}
          <div className="flex items-center justify-between border-t border-seagreen-500/10 px-4 py-2 text-xs text-muted-foreground">
            <div className="flex items-center gap-4">
              <span className="flex items-center gap-1">
                <kbd className="px-1.5 py-0.5 rounded bg-muted/50 font-mono">↑↓</kbd>
                Navigate
              </span>
              <span className="flex items-center gap-1">
                <kbd className="px-1.5 py-0.5 rounded bg-muted/50 font-mono">↵</kbd>
                Select
              </span>
            </div>
            <span className="flex items-center gap-1">
              <Command className="h-3 w-3" />
              <span>K to open</span>
            </span>
          </div>
        </div>
      </div>
    </>
  )
}

// Individual Command Item Component
function CommandItem({
  command,
  isSelected,
  dataIndex,
  onClick,
  onMouseEnter
}: {
  command: CommandItem
  isSelected: boolean
  dataIndex: number
  onClick: () => void
  onMouseEnter: () => void
}) {
  const Icon = command.icon

  return (
    <button
      data-index={dataIndex}
      onClick={onClick}
      onMouseEnter={onMouseEnter}
      className={cn(
        "flex w-full items-center gap-3 rounded-xl px-3 py-2.5",
        "transition-all duration-150",
        isSelected
          ? "bg-gradient-to-r from-twilight-800/10 via-seagreen-500/10 to-mint-500/10 text-foreground"
          : "text-muted-foreground hover:text-foreground"
      )}
    >
      <div className={cn(
        "flex h-9 w-9 items-center justify-center rounded-lg",
        "transition-all duration-150",
        isSelected
          ? "bg-gradient-to-br from-twilight-800 to-seagreen-500 text-white shadow-lg shadow-seagreen-500/25"
          : "bg-muted/50"
      )}>
        <Icon className="h-4 w-4" />
      </div>

      <div className="flex-1 text-left">
        <p className="text-sm font-medium">{command.label}</p>
        <p className="text-xs text-muted-foreground">{command.description}</p>
      </div>

      {command.shortcut && (
        <kbd className={cn(
          "hidden sm:flex px-2 py-1 rounded-lg",
          "bg-muted/50 text-xs font-mono text-muted-foreground"
        )}>
          {command.shortcut}
        </kbd>
      )}

      {isSelected && (
        <ArrowRight className="h-4 w-4 text-seagreen-500" />
      )}
    </button>
  )
}

// Hook for using command palette
export function useCommandPalette() {
  const [isOpen, setIsOpen] = useState(false)

  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      if ((e.metaKey || e.ctrlKey) && e.key === "k") {
        e.preventDefault()
        setIsOpen(prev => !prev)
      }
    }

    document.addEventListener("keydown", handleKeyDown)
    return () => document.removeEventListener("keydown", handleKeyDown)
  }, [])

  return {
    isOpen,
    open: () => setIsOpen(true),
    close: () => setIsOpen(false),
    toggle: () => setIsOpen(prev => !prev)
  }
}
