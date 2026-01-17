"use client"

import * as React from "react"
import { createContext, useContext, useState, useCallback } from "react"
import { X, CheckCircle2, XCircle, AlertCircle, Info, Loader2 } from "lucide-react"
import { cn } from "@/lib/utils"

// Toast types
export type ToastType = "success" | "error" | "warning" | "info" | "loading"

export interface Toast {
  id: string
  type: ToastType
  title: string
  description?: string
  duration?: number
  action?: {
    label: string
    onClick: () => void
  }
}

interface ToastContextValue {
  toasts: Toast[]
  addToast: (toast: Omit<Toast, "id">) => string
  removeToast: (id: string) => void
  updateToast: (id: string, toast: Partial<Omit<Toast, "id">>) => void
}

const ToastContext = createContext<ToastContextValue | null>(null)

// Toast Provider
export function ToastProvider({ children }: { children: React.ReactNode }) {
  const [toasts, setToasts] = useState<Toast[]>([])

  const addToast = useCallback((toast: Omit<Toast, "id">) => {
    const id = Math.random().toString(36).substring(2, 9)
    const newToast: Toast = {
      ...toast,
      id,
      duration: toast.duration ?? (toast.type === "loading" ? Infinity : 5000),
    }

    setToasts((prev) => [...prev, newToast])

    // Auto-remove after duration (unless loading or Infinity)
    if (newToast.duration !== Infinity) {
      setTimeout(() => {
        setToasts((prev) => prev.filter((t) => t.id !== id))
      }, newToast.duration)
    }

    return id
  }, [])

  const removeToast = useCallback((id: string) => {
    setToasts((prev) => prev.filter((t) => t.id !== id))
  }, [])

  const updateToast = useCallback((id: string, updates: Partial<Omit<Toast, "id">>) => {
    setToasts((prev) =>
      prev.map((t) => (t.id === id ? { ...t, ...updates } : t))
    )

    // If updating to non-loading type, set auto-remove
    if (updates.type && updates.type !== "loading") {
      const duration = updates.duration ?? 5000
      setTimeout(() => {
        setToasts((prev) => prev.filter((t) => t.id !== id))
      }, duration)
    }
  }, [])

  return (
    <ToastContext.Provider value={{ toasts, addToast, removeToast, updateToast }}>
      {children}
      <ToastContainer />
    </ToastContext.Provider>
  )
}

// Hook to use toast
export function useToast() {
  const context = useContext(ToastContext)
  if (!context) {
    throw new Error("useToast must be used within a ToastProvider")
  }

  const { addToast, removeToast, updateToast } = context

  return {
    toast: addToast,
    dismiss: removeToast,
    update: updateToast,
    // Convenience methods
    success: (title: string, description?: string) =>
      addToast({ type: "success", title, description }),
    error: (title: string, description?: string) =>
      addToast({ type: "error", title, description }),
    warning: (title: string, description?: string) =>
      addToast({ type: "warning", title, description }),
    info: (title: string, description?: string) =>
      addToast({ type: "info", title, description }),
    loading: (title: string, description?: string) =>
      addToast({ type: "loading", title, description }),
    // Promise helper
    promise: async <T,>(
      promise: Promise<T>,
      messages: {
        loading: string
        success: string
        error: string
      }
    ) => {
      const id = addToast({ type: "loading", title: messages.loading })
      try {
        const result = await promise
        updateToast(id, { type: "success", title: messages.success })
        return result
      } catch (err) {
        updateToast(id, { type: "error", title: messages.error })
        throw err
      }
    },
  }
}

// Toast Container
function ToastContainer() {
  const context = useContext(ToastContext)
  if (!context) return null

  const { toasts, removeToast } = context

  return (
    <div
      className={cn(
        "fixed bottom-4 right-4 z-[100]",
        "flex flex-col gap-2",
        "max-w-sm w-full",
        "pointer-events-none"
      )}
      aria-live="polite"
      aria-label="Notifications"
    >
      {toasts.map((toast) => (
        <ToastItem key={toast.id} toast={toast} onDismiss={() => removeToast(toast.id)} />
      ))}
    </div>
  )
}

// Individual Toast Item
function ToastItem({ toast, onDismiss }: { toast: Toast; onDismiss: () => void }) {
  const icons = {
    success: <CheckCircle2 className="h-5 w-5 text-mint-500" />,
    error: <XCircle className="h-5 w-5 text-red-500" />,
    warning: <AlertCircle className="h-5 w-5 text-golden-500" />,
    info: <Info className="h-5 w-5 text-seagreen-500" />,
    loading: <Loader2 className="h-5 w-5 text-seagreen-500 animate-spin" />,
  }

  const bgColors = {
    success: "border-mint-500/30 bg-mint-500/5",
    error: "border-red-500/30 bg-red-500/5",
    warning: "border-golden-500/30 bg-golden-500/5",
    info: "border-seagreen-500/30 bg-seagreen-500/5",
    loading: "border-seagreen-500/30 bg-seagreen-500/5",
  }

  return (
    <div
      className={cn(
        "pointer-events-auto",
        "flex items-start gap-3 p-4",
        "rounded-xl border backdrop-blur-xl",
        "bg-white/90 dark:bg-slate-900/90",
        "shadow-lg shadow-black/10",
        "animate-slide-in-right",
        bgColors[toast.type]
      )}
      role="alert"
    >
      {/* Icon */}
      <div className="flex-shrink-0 mt-0.5">{icons[toast.type]}</div>

      {/* Content */}
      <div className="flex-1 min-w-0">
        <p className="font-medium text-sm">{toast.title}</p>
        {toast.description && (
          <p className="text-sm text-muted-foreground mt-1">{toast.description}</p>
        )}
        {toast.action && (
          <button
            onClick={toast.action.onClick}
            className="mt-2 text-sm font-medium text-seagreen-600 dark:text-seagreen-400 hover:underline"
          >
            {toast.action.label}
          </button>
        )}
      </div>

      {/* Dismiss button */}
      {toast.type !== "loading" && (
        <button
          onClick={onDismiss}
          className={cn(
            "flex-shrink-0 p-1 rounded-lg",
            "text-muted-foreground hover:text-foreground",
            "hover:bg-muted/50 transition-colors"
          )}
          aria-label="Dismiss notification"
        >
          <X className="h-4 w-4" />
        </button>
      )}
    </div>
  )
}

// Add animation to globals.css
// @keyframes slide-in-right {
//   from { transform: translateX(100%); opacity: 0; }
//   to { transform: translateX(0); opacity: 1; }
// }
// .animate-slide-in-right { animation: slide-in-right 0.3s ease-out; }
