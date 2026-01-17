"use client"

import React from "react"
import { AlertTriangle, RefreshCw, Home, Bug, ChevronDown, ChevronUp } from "lucide-react"
import { Button } from "./ui/button"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "./ui/card"
import { cn } from "@/lib/utils"

interface ErrorBoundaryState {
  hasError: boolean
  error: Error | null
  errorInfo: React.ErrorInfo | null
  showDetails: boolean
}

interface ErrorBoundaryProps {
  children: React.ReactNode
  fallback?: React.ReactNode
  onError?: (error: Error, errorInfo: React.ErrorInfo) => void
  onReset?: () => void
}

export class ErrorBoundary extends React.Component<ErrorBoundaryProps, ErrorBoundaryState> {
  constructor(props: ErrorBoundaryProps) {
    super(props)
    this.state = {
      hasError: false,
      error: null,
      errorInfo: null,
      showDetails: false,
    }
  }

  static getDerivedStateFromError(error: Error): Partial<ErrorBoundaryState> {
    return { hasError: true, error }
  }

  componentDidCatch(error: Error, errorInfo: React.ErrorInfo) {
    this.setState({ errorInfo })

    // Call optional error handler
    this.props.onError?.(error, errorInfo)

    // Log to console in development
    if (process.env.NODE_ENV === "development") {
      console.error("Error Boundary caught an error:", error, errorInfo)
    }
  }

  handleReset = () => {
    this.props.onReset?.()
    this.setState({
      hasError: false,
      error: null,
      errorInfo: null,
      showDetails: false,
    })
  }

  toggleDetails = () => {
    this.setState((prev) => ({ showDetails: !prev.showDetails }))
  }

  render() {
    if (this.state.hasError) {
      // Custom fallback provided
      if (this.props.fallback) {
        return this.props.fallback
      }

      // Default error UI
      return (
        <div className="min-h-[400px] flex items-center justify-center p-6 animate-fade-slide-up">
          <Card className="max-w-lg w-full border-red-500/20 bg-gradient-to-br from-red-500/5 via-transparent to-transparent">
            <CardHeader className="text-center">
              <div className="mx-auto mb-4 flex h-16 w-16 items-center justify-center rounded-2xl bg-red-500/10 border border-red-500/20">
                <AlertTriangle className="h-8 w-8 text-red-500" />
              </div>
              <CardTitle className="text-xl">Something went wrong</CardTitle>
              <CardDescription>
                An unexpected error occurred. We apologize for the inconvenience.
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              {/* Error message preview */}
              {this.state.error && (
                <div className="p-3 rounded-lg bg-muted/50 border border-border">
                  <p className="text-sm font-mono text-muted-foreground truncate">
                    {this.state.error.message || "Unknown error"}
                  </p>
                </div>
              )}

              {/* Action buttons */}
              <div className="flex flex-col sm:flex-row gap-2">
                <Button
                  onClick={this.handleReset}
                  className="flex-1 bg-gradient-to-r from-twilight-800 to-seagreen-500 text-white"
                >
                  <RefreshCw className="h-4 w-4 mr-2" />
                  Try Again
                </Button>
                <Button
                  variant="outline"
                  onClick={() => window.location.href = "/"}
                  className="flex-1"
                >
                  <Home className="h-4 w-4 mr-2" />
                  Go Home
                </Button>
              </div>

              {/* Technical details toggle (development mode) */}
              {process.env.NODE_ENV === "development" && this.state.errorInfo && (
                <>
                  <button
                    onClick={this.toggleDetails}
                    className="flex items-center gap-2 text-sm text-muted-foreground hover:text-foreground transition-colors w-full justify-center"
                  >
                    <Bug className="h-4 w-4" />
                    Technical Details
                    {this.state.showDetails ? (
                      <ChevronUp className="h-4 w-4" />
                    ) : (
                      <ChevronDown className="h-4 w-4" />
                    )}
                  </button>

                  {this.state.showDetails && (
                    <div className="p-3 rounded-lg bg-muted/30 border border-border overflow-x-auto">
                      <pre className="text-xs font-mono text-muted-foreground whitespace-pre-wrap break-words">
                        {this.state.error?.stack}
                        {"\n\nComponent Stack:"}
                        {this.state.errorInfo.componentStack}
                      </pre>
                    </div>
                  )}
                </>
              )}
            </CardContent>
          </Card>
        </div>
      )
    }

    return this.props.children
  }
}

// Hook for functional components to trigger error boundary
export function useErrorBoundary() {
  const [error, setError] = React.useState<Error | null>(null)

  if (error) {
    throw error
  }

  return {
    showBoundary: setError,
    resetBoundary: () => setError(null),
  }
}

// Specialized error boundaries for different sections
export function PageErrorBoundary({ children }: { children: React.ReactNode }) {
  return (
    <ErrorBoundary
      fallback={
        <div className="min-h-screen flex items-center justify-center p-6 bg-gradient-to-br from-background to-muted/20">
          <Card className="max-w-md w-full border-red-500/20">
            <CardHeader className="text-center">
              <div className="mx-auto mb-4 flex h-20 w-20 items-center justify-center rounded-2xl bg-red-500/10">
                <AlertTriangle className="h-10 w-10 text-red-500" />
              </div>
              <CardTitle className="text-2xl">Page Error</CardTitle>
              <CardDescription>
                This page encountered an error and couldn&apos;t be displayed.
              </CardDescription>
            </CardHeader>
            <CardContent className="flex flex-col gap-2">
              <Button
                onClick={() => window.location.reload()}
                className="bg-gradient-to-r from-twilight-800 to-seagreen-500 text-white"
              >
                <RefreshCw className="h-4 w-4 mr-2" />
                Reload Page
              </Button>
              <Button
                variant="outline"
                onClick={() => window.location.href = "/"}
              >
                <Home className="h-4 w-4 mr-2" />
                Return to Dashboard
              </Button>
            </CardContent>
          </Card>
        </div>
      }
    >
      {children}
    </ErrorBoundary>
  )
}

export function ComponentErrorBoundary({
  children,
  name = "Component",
}: {
  children: React.ReactNode
  name?: string
}) {
  return (
    <ErrorBoundary
      fallback={
        <div className="p-4 rounded-xl border border-red-500/20 bg-red-500/5">
          <div className="flex items-center gap-3">
            <div className="flex h-10 w-10 items-center justify-center rounded-lg bg-red-500/10">
              <AlertTriangle className="h-5 w-5 text-red-500" />
            </div>
            <div className="flex-1">
              <p className="font-medium text-sm">{name} failed to load</p>
              <p className="text-xs text-muted-foreground">
                Please try refreshing the page
              </p>
            </div>
            <Button
              size="sm"
              variant="outline"
              onClick={() => window.location.reload()}
            >
              <RefreshCw className="h-3 w-3 mr-1" />
              Retry
            </Button>
          </div>
        </div>
      }
    >
      {children}
    </ErrorBoundary>
  )
}
