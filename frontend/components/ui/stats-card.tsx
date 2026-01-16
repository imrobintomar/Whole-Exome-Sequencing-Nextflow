'use client';

import * as React from 'react';
import { cn } from '@/lib/utils';
import { LucideIcon } from 'lucide-react';

interface StatsCardProps extends React.HTMLAttributes<HTMLDivElement> {
  title: string;
  value: number | string;
  icon: LucideIcon;
  trend?: {
    value: number;
    direction: 'up' | 'down' | 'neutral';
  };
  gradient: 'purple' | 'green' | 'blue' | 'red' | 'amber' | 'cyan' | 'teal';
  description?: string;
  loading?: boolean;
}

const gradientMap = {
  purple: {
    bg: 'bg-gradient-to-br from-purple-500/10 via-purple-500/5 to-transparent dark:from-purple-500/20 dark:via-purple-500/10',
    iconBg: 'bg-gradient-to-br from-purple-500 to-purple-600',
    iconGlow: 'shadow-[0_0_20px_rgba(139,92,246,0.4)]',
    border: 'border-purple-500/20',
    text: 'text-purple-600 dark:text-purple-400',
  },
  green: {
    bg: 'bg-gradient-to-br from-emerald-500/10 via-emerald-500/5 to-transparent dark:from-emerald-500/20 dark:via-emerald-500/10',
    iconBg: 'bg-gradient-to-br from-emerald-500 to-teal-500',
    iconGlow: 'shadow-[0_0_20px_rgba(16,185,129,0.4)]',
    border: 'border-emerald-500/20',
    text: 'text-emerald-600 dark:text-emerald-400',
  },
  blue: {
    bg: 'bg-gradient-to-br from-blue-500/10 via-blue-500/5 to-transparent dark:from-blue-500/20 dark:via-blue-500/10',
    iconBg: 'bg-gradient-to-br from-blue-500 to-cyan-500',
    iconGlow: 'shadow-[0_0_20px_rgba(59,130,246,0.4)]',
    border: 'border-blue-500/20',
    text: 'text-blue-600 dark:text-blue-400',
  },
  red: {
    bg: 'bg-gradient-to-br from-red-500/10 via-red-500/5 to-transparent dark:from-red-500/20 dark:via-red-500/10',
    iconBg: 'bg-gradient-to-br from-red-500 to-rose-500',
    iconGlow: 'shadow-[0_0_20px_rgba(239,68,68,0.4)]',
    border: 'border-red-500/20',
    text: 'text-red-600 dark:text-red-400',
  },
  amber: {
    bg: 'bg-gradient-to-br from-amber-500/10 via-amber-500/5 to-transparent dark:from-amber-500/20 dark:via-amber-500/10',
    iconBg: 'bg-gradient-to-br from-amber-500 to-orange-500',
    iconGlow: 'shadow-[0_0_20px_rgba(245,158,11,0.4)]',
    border: 'border-amber-500/20',
    text: 'text-amber-600 dark:text-amber-400',
  },
  cyan: {
    bg: 'bg-gradient-to-br from-cyan-500/10 via-cyan-500/5 to-transparent dark:from-cyan-500/20 dark:via-cyan-500/10',
    iconBg: 'bg-gradient-to-br from-cyan-500 to-blue-500',
    iconGlow: 'shadow-[0_0_20px_rgba(6,182,212,0.4)]',
    border: 'border-cyan-500/20',
    text: 'text-cyan-600 dark:text-cyan-400',
  },
  teal: {
    bg: 'bg-gradient-to-br from-teal-500/10 via-teal-500/5 to-transparent dark:from-teal-500/20 dark:via-teal-500/10',
    iconBg: 'bg-gradient-to-br from-teal-500 to-emerald-500',
    iconGlow: 'shadow-[0_0_20px_rgba(20,184,166,0.4)]',
    border: 'border-teal-500/20',
    text: 'text-teal-600 dark:text-teal-400',
  },
};

const TrendIndicator = ({ value, direction }: { value: number; direction: 'up' | 'down' | 'neutral' }) => {
  const colors = {
    up: 'text-emerald-600 dark:text-emerald-400 bg-emerald-100 dark:bg-emerald-900/30',
    down: 'text-red-600 dark:text-red-400 bg-red-100 dark:bg-red-900/30',
    neutral: 'text-gray-600 dark:text-gray-400 bg-gray-100 dark:bg-gray-900/30',
  };

  const icons = {
    up: (
      <svg className="w-3 h-3" fill="none" viewBox="0 0 24 24" stroke="currentColor">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 10l7-7m0 0l7 7m-7-7v18" />
      </svg>
    ),
    down: (
      <svg className="w-3 h-3" fill="none" viewBox="0 0 24 24" stroke="currentColor">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 14l-7 7m0 0l-7-7m7 7V3" />
      </svg>
    ),
    neutral: (
      <svg className="w-3 h-3" fill="none" viewBox="0 0 24 24" stroke="currentColor">
        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 12h14" />
      </svg>
    ),
  };

  return (
    <span className={cn('inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-xs font-medium', colors[direction])}>
      {icons[direction]}
      {Math.abs(value)}%
    </span>
  );
};

const SkeletonLoader = () => (
  <div className="animate-pulse space-y-3">
    <div className="h-4 w-20 bg-muted rounded"></div>
    <div className="h-8 w-24 bg-muted rounded"></div>
    <div className="h-3 w-16 bg-muted rounded"></div>
  </div>
);

const StatsCard = React.forwardRef<HTMLDivElement, StatsCardProps>(
  ({ className, title, value, icon: Icon, trend, gradient, description, loading = false, ...props }, ref) => {
    const colors = gradientMap[gradient];

    return (
      <div
        ref={ref}
        className={cn(
          'relative overflow-hidden rounded-2xl border p-6',
          'backdrop-blur-xl transition-all duration-300 ease-out',
          'hover:-translate-y-1 hover:shadow-xl',
          colors.bg,
          colors.border,
          className
        )}
        {...props}
      >
        {/* Background decoration */}
        <div className="absolute -right-4 -top-4 h-24 w-24 rounded-full bg-gradient-to-br from-white/5 to-transparent" />

        <div className="relative flex items-start justify-between">
          <div className="flex-1">
            {loading ? (
              <SkeletonLoader />
            ) : (
              <>
                <p className="text-sm font-medium text-muted-foreground mb-1">{title}</p>
                <div className="flex items-baseline gap-2">
                  <h3 className="text-3xl font-bold tracking-tight animate-count-up">{value}</h3>
                  {trend && <TrendIndicator value={trend.value} direction={trend.direction} />}
                </div>
                {description && (
                  <p className={cn('text-xs mt-1', colors.text)}>{description}</p>
                )}
              </>
            )}
          </div>

          {/* Icon container */}
          <div
            className={cn(
              'flex h-12 w-12 items-center justify-center rounded-xl',
              colors.iconBg,
              colors.iconGlow,
              'transition-transform duration-300 group-hover:scale-110'
            )}
          >
            <Icon className="h-6 w-6 text-white" />
          </div>
        </div>
      </div>
    );
  }
);

StatsCard.displayName = 'StatsCard';

export { StatsCard };
