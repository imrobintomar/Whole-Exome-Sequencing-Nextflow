/**
 * Accessibility utilities for ATGCFlow
 * Provides hooks and helpers for keyboard navigation, focus management, and screen reader support
 */

import { useCallback, useEffect, useRef } from "react"

/**
 * Hook for trapping focus within a container (useful for modals, dialogs)
 */
export function useFocusTrap(isActive: boolean) {
  const containerRef = useRef<HTMLDivElement>(null)

  useEffect(() => {
    if (!isActive || !containerRef.current) return

    const container = containerRef.current
    const focusableElements = container.querySelectorAll<HTMLElement>(
      'button, [href], input, select, textarea, [tabindex]:not([tabindex="-1"])'
    )
    const firstElement = focusableElements[0]
    const lastElement = focusableElements[focusableElements.length - 1]

    // Focus first element when trap activates
    firstElement?.focus()

    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.key !== "Tab") return

      if (e.shiftKey) {
        // Shift + Tab
        if (document.activeElement === firstElement) {
          e.preventDefault()
          lastElement?.focus()
        }
      } else {
        // Tab
        if (document.activeElement === lastElement) {
          e.preventDefault()
          firstElement?.focus()
        }
      }
    }

    container.addEventListener("keydown", handleKeyDown)
    return () => container.removeEventListener("keydown", handleKeyDown)
  }, [isActive])

  return containerRef
}

/**
 * Hook for keyboard list navigation (arrow keys)
 */
export function useListNavigation<T extends HTMLElement>(
  items: T[],
  options: {
    loop?: boolean
    orientation?: "horizontal" | "vertical" | "both"
    onSelect?: (item: T, index: number) => void
  } = {}
) {
  const { loop = true, orientation = "vertical", onSelect } = options
  const currentIndex = useRef(0)

  const getNextIndex = useCallback(
    (direction: "next" | "prev") => {
      const length = items.length
      if (direction === "next") {
        return loop
          ? (currentIndex.current + 1) % length
          : Math.min(currentIndex.current + 1, length - 1)
      }
      return loop
        ? (currentIndex.current - 1 + length) % length
        : Math.max(currentIndex.current - 1, 0)
    },
    [items.length, loop]
  )

  const handleKeyDown = useCallback(
    (e: KeyboardEvent) => {
      const isVertical = orientation === "vertical" || orientation === "both"
      const isHorizontal = orientation === "horizontal" || orientation === "both"

      let newIndex: number | null = null

      switch (e.key) {
        case "ArrowDown":
          if (isVertical) {
            e.preventDefault()
            newIndex = getNextIndex("next")
          }
          break
        case "ArrowUp":
          if (isVertical) {
            e.preventDefault()
            newIndex = getNextIndex("prev")
          }
          break
        case "ArrowRight":
          if (isHorizontal) {
            e.preventDefault()
            newIndex = getNextIndex("next")
          }
          break
        case "ArrowLeft":
          if (isHorizontal) {
            e.preventDefault()
            newIndex = getNextIndex("prev")
          }
          break
        case "Home":
          e.preventDefault()
          newIndex = 0
          break
        case "End":
          e.preventDefault()
          newIndex = items.length - 1
          break
        case "Enter":
        case " ":
          e.preventDefault()
          onSelect?.(items[currentIndex.current], currentIndex.current)
          break
      }

      if (newIndex !== null) {
        currentIndex.current = newIndex
        items[newIndex]?.focus()
      }
    },
    [items, orientation, getNextIndex, onSelect]
  )

  useEffect(() => {
    items.forEach((item) => {
      item.addEventListener("keydown", handleKeyDown)
    })

    return () => {
      items.forEach((item) => {
        item.removeEventListener("keydown", handleKeyDown)
      })
    }
  }, [items, handleKeyDown])

  return {
    focusFirst: () => {
      currentIndex.current = 0
      items[0]?.focus()
    },
    focusLast: () => {
      currentIndex.current = items.length - 1
      items[items.length - 1]?.focus()
    },
    getCurrentIndex: () => currentIndex.current,
  }
}

/**
 * Hook for announcing content to screen readers
 */
export function useAnnounce() {
  const announce = useCallback((message: string, priority: "polite" | "assertive" = "polite") => {
    const announcement = document.createElement("div")
    announcement.setAttribute("role", "status")
    announcement.setAttribute("aria-live", priority)
    announcement.setAttribute("aria-atomic", "true")
    announcement.className = "sr-only"
    announcement.textContent = message

    document.body.appendChild(announcement)

    // Remove after announcement is read
    setTimeout(() => {
      document.body.removeChild(announcement)
    }, 1000)
  }, [])

  return announce
}

/**
 * Hook for managing focus restoration (e.g., after closing a modal)
 */
export function useFocusRestoration() {
  const previousFocus = useRef<HTMLElement | null>(null)

  const saveFocus = useCallback(() => {
    previousFocus.current = document.activeElement as HTMLElement
  }, [])

  const restoreFocus = useCallback(() => {
    previousFocus.current?.focus()
  }, [])

  return { saveFocus, restoreFocus }
}

/**
 * Hook for skip link functionality
 */
export function useSkipToContent(contentId: string = "main-content") {
  const skipToContent = useCallback(() => {
    const content = document.getElementById(contentId)
    if (content) {
      content.tabIndex = -1
      content.focus()
      content.scrollIntoView()
    }
  }, [contentId])

  return skipToContent
}

/**
 * Hook for detecting reduced motion preference
 */
export function useReducedMotion() {
  const prefersReducedMotion =
    typeof window !== "undefined"
      ? window.matchMedia("(prefers-reduced-motion: reduce)").matches
      : false

  return prefersReducedMotion
}

/**
 * Generate unique IDs for ARIA relationships
 */
let idCounter = 0
export function generateId(prefix: string = "atgc"): string {
  return `${prefix}-${++idCounter}`
}

/**
 * Keyboard shortcut helper
 */
export function formatShortcut(keys: string[]): string {
  const isMac = typeof navigator !== "undefined" && navigator.platform.toLowerCase().includes("mac")

  const keyMap: Record<string, string> = {
    mod: isMac ? "⌘" : "Ctrl",
    shift: isMac ? "⇧" : "Shift",
    alt: isMac ? "⌥" : "Alt",
    ctrl: isMac ? "⌃" : "Ctrl",
    enter: "↵",
    escape: "Esc",
    arrowup: "↑",
    arrowdown: "↓",
    arrowleft: "←",
    arrowright: "→",
  }

  return keys
    .map((key) => keyMap[key.toLowerCase()] || key.toUpperCase())
    .join(isMac ? "" : "+")
}

/**
 * Screen reader only text component helper
 */
export const srOnlyStyles = {
  position: "absolute" as const,
  width: "1px",
  height: "1px",
  padding: "0",
  margin: "-1px",
  overflow: "hidden",
  clip: "rect(0, 0, 0, 0)",
  whiteSpace: "nowrap" as const,
  border: "0",
}
