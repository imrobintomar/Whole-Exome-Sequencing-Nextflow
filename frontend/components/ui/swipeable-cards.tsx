"use client"

import * as React from "react"
import { useState, useRef, useEffect, useCallback } from "react"
import { cn } from "@/lib/utils"
import { ChevronLeft, ChevronRight } from "lucide-react"

interface SwipeableCardsProps {
  children: React.ReactNode[]
  className?: string
  showArrows?: boolean
  showDots?: boolean
  autoPlay?: boolean
  autoPlayInterval?: number
}

export function SwipeableCards({
  children,
  className,
  showArrows = true,
  showDots = true,
  autoPlay = false,
  autoPlayInterval = 5000,
}: SwipeableCardsProps) {
  const [currentIndex, setCurrentIndex] = useState(0)
  const [isDragging, setIsDragging] = useState(false)
  const [startX, setStartX] = useState(0)
  const [translateX, setTranslateX] = useState(0)
  const containerRef = useRef<HTMLDivElement>(null)
  const cardCount = React.Children.count(children)

  // Touch/Mouse handlers
  const handleDragStart = (clientX: number) => {
    setIsDragging(true)
    setStartX(clientX)
  }

  const handleDragMove = (clientX: number) => {
    if (!isDragging) return
    const diff = clientX - startX
    setTranslateX(diff)
  }

  const handleDragEnd = () => {
    if (!isDragging) return
    setIsDragging(false)

    const threshold = 50
    if (translateX > threshold && currentIndex > 0) {
      setCurrentIndex(prev => prev - 1)
    } else if (translateX < -threshold && currentIndex < cardCount - 1) {
      setCurrentIndex(prev => prev + 1)
    }
    setTranslateX(0)
  }

  // Mouse events
  const onMouseDown = (e: React.MouseEvent) => handleDragStart(e.clientX)
  const onMouseMove = (e: React.MouseEvent) => handleDragMove(e.clientX)
  const onMouseUp = () => handleDragEnd()
  const onMouseLeave = () => handleDragEnd()

  // Touch events
  const onTouchStart = (e: React.TouchEvent) => handleDragStart(e.touches[0].clientX)
  const onTouchMove = (e: React.TouchEvent) => handleDragMove(e.touches[0].clientX)
  const onTouchEnd = () => handleDragEnd()

  // Navigation
  const goToSlide = useCallback((index: number) => {
    setCurrentIndex(Math.max(0, Math.min(index, cardCount - 1)))
  }, [cardCount])

  const goToPrevious = useCallback(() => {
    goToSlide(currentIndex - 1)
  }, [currentIndex, goToSlide])

  const goToNext = useCallback(() => {
    goToSlide(currentIndex + 1)
  }, [currentIndex, goToSlide])

  // Auto-play
  useEffect(() => {
    if (!autoPlay) return

    const interval = setInterval(() => {
      setCurrentIndex(prev => (prev + 1) % cardCount)
    }, autoPlayInterval)

    return () => clearInterval(interval)
  }, [autoPlay, autoPlayInterval, cardCount])

  // Keyboard navigation
  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.key === "ArrowLeft") goToPrevious()
      if (e.key === "ArrowRight") goToNext()
    }

    const container = containerRef.current
    if (container) {
      container.addEventListener("keydown", handleKeyDown)
      return () => container.removeEventListener("keydown", handleKeyDown)
    }
  }, [goToPrevious, goToNext])

  return (
    <div
      ref={containerRef}
      className={cn("relative", className)}
      tabIndex={0}
      role="region"
      aria-label="Swipeable cards"
    >
      {/* Cards Container */}
      <div
        className="overflow-hidden rounded-2xl"
        onMouseDown={onMouseDown}
        onMouseMove={onMouseMove}
        onMouseUp={onMouseUp}
        onMouseLeave={onMouseLeave}
        onTouchStart={onTouchStart}
        onTouchMove={onTouchMove}
        onTouchEnd={onTouchEnd}
      >
        <div
          className={cn(
            "flex transition-transform duration-300 ease-out",
            isDragging && "transition-none"
          )}
          style={{
            transform: `translateX(calc(-${currentIndex * 100}% + ${translateX}px))`,
          }}
        >
          {React.Children.map(children, (child, index) => (
            <div
              key={index}
              className="w-full flex-shrink-0 px-1"
              aria-hidden={index !== currentIndex}
            >
              {child}
            </div>
          ))}
        </div>
      </div>

      {/* Navigation Arrows (Desktop/Tablet) */}
      {showArrows && cardCount > 1 && (
        <>
          <button
            onClick={goToPrevious}
            disabled={currentIndex === 0}
            className={cn(
              "hidden sm:flex absolute left-0 top-1/2 -translate-y-1/2 -translate-x-1/2",
              "h-10 w-10 items-center justify-center rounded-full",
              "bg-white/90 dark:bg-twilight-800/90 backdrop-blur-sm",
              "border border-seagreen-500/20",
              "shadow-lg shadow-twilight-800/10",
              "text-foreground",
              "transition-all duration-200",
              "hover:scale-110 hover:shadow-xl",
              "disabled:opacity-50 disabled:cursor-not-allowed disabled:hover:scale-100"
            )}
            aria-label="Previous card"
          >
            <ChevronLeft className="h-5 w-5" />
          </button>

          <button
            onClick={goToNext}
            disabled={currentIndex === cardCount - 1}
            className={cn(
              "hidden sm:flex absolute right-0 top-1/2 -translate-y-1/2 translate-x-1/2",
              "h-10 w-10 items-center justify-center rounded-full",
              "bg-white/90 dark:bg-twilight-800/90 backdrop-blur-sm",
              "border border-seagreen-500/20",
              "shadow-lg shadow-twilight-800/10",
              "text-foreground",
              "transition-all duration-200",
              "hover:scale-110 hover:shadow-xl",
              "disabled:opacity-50 disabled:cursor-not-allowed disabled:hover:scale-100"
            )}
            aria-label="Next card"
          >
            <ChevronRight className="h-5 w-5" />
          </button>
        </>
      )}

      {/* Dots Indicator */}
      {showDots && cardCount > 1 && (
        <div className="flex justify-center gap-2 mt-4">
          {Array.from({ length: cardCount }).map((_, index) => (
            <button
              key={index}
              onClick={() => goToSlide(index)}
              className={cn(
                "h-2 rounded-full transition-all duration-300",
                index === currentIndex
                  ? "w-6 bg-gradient-to-r from-seagreen-500 to-mint-500"
                  : "w-2 bg-muted-foreground/30 hover:bg-muted-foreground/50"
              )}
              aria-label={`Go to card ${index + 1}`}
              aria-current={index === currentIndex ? "true" : "false"}
            />
          ))}
        </div>
      )}

      {/* Swipe Hint (Mobile only, first render) */}
      <div className="sm:hidden flex justify-center mt-2">
        <span className="text-xs text-muted-foreground flex items-center gap-1">
          <ChevronLeft className="h-3 w-3" />
          Swipe
          <ChevronRight className="h-3 w-3" />
        </span>
      </div>
    </div>
  )
}

// Mobile-specific Stats Cards Wrapper
export function MobileStatsCarousel({
  children,
  className
}: {
  children: React.ReactNode[]
  className?: string
}) {
  return (
    <div className={cn("md:hidden", className)}>
      <SwipeableCards showArrows={false}>
        {children}
      </SwipeableCards>
    </div>
  )
}

// Desktop Grid (shown on larger screens)
export function DesktopStatsGrid({
  children,
  className
}: {
  children: React.ReactNode
  className?: string
}) {
  return (
    <div className={cn(
      "hidden md:grid gap-4 md:grid-cols-2 lg:grid-cols-4",
      className
    )}>
      {children}
    </div>
  )
}
