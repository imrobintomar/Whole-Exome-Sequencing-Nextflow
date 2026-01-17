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

/* Brand Colors:
 * twilight: #110B52 (deep purple)
 * seagreen: #00A0A0 (light sea green)
 * mint: #00E897 (tropical mint)
 * turf: #007F4F (turf green)
 * golden: #F2D513 (golden glow)
 */
const gradientMap = {
  purple: {
    bg: 'bg-gradient-to-br from-twilight-800/10 via-twilight-800/5 to-transparent dark:from-twilight-800/20 dark:via-twilight-800/10',
    iconBg: 'bg-gradient-to-br from-twilight-800 to-seagreen-500',
    iconGlow: 'shadow-[0_0_20px_rgba(17,11,82,0.4)]',
    border: 'border-twilight-800/20',
    text: 'text-twilight-800 dark:text-seagreen-400',
  },
  green: {
    bg: 'bg-gradient-to-br from-mint-500/10 via-mint-500/5 to-transparent dark:from-mint-500/20 dark:via-mint-500/10',
    iconBg: 'bg-gradient-to-br from-mint-500 to-mint-800',
    iconGlow: 'shadow-[0_0_20px_rgba(0,232,151,0.4)]',
    border: 'border-mint-500/20',
    text: 'text-mint-800 dark:text-mint-400',
  },
  blue: {
    bg: 'bg-gradient-to-br from-seagreen-500/10 via-seagreen-500/5 to-transparent dark:from-seagreen-500/20 dark:via-seagreen-500/10',
    iconBg: 'bg-gradient-to-br from-seagreen-500 to-mint-500',
    iconGlow: 'shadow-[0_0_20px_rgba(0,160,160,0.4)]',
    border: 'border-seagreen-500/20',
    text: 'text-seagreen-600 dark:text-seagreen-400',
  },
  red: {
    bg: 'bg-gradient-to-br from-red-500/10 via-red-500/5 to-transparent dark:from-red-500/20 dark:via-red-500/10',
    iconBg: 'bg-gradient-to-br from-red-500 to-rose-500',
    iconGlow: 'shadow-[0_0_20px_rgba(239,68,68,0.4)]',
    border: 'border-red-500/20',
    text: 'text-red-600 dark:text-red-400',
  },
  amber: {
    bg: 'bg-gradient-to-br from-golden-500/10 via-golden-500/5 to-transparent dark:from-golden-500/20 dark:via-golden-500/10',
    iconBg: 'bg-gradient-to-br from-golden-500 to-mint-800',
    iconGlow: 'shadow-[0_0_20px_rgba(242,213,19,0.4)]',
    border: 'border-golden-500/20',
    text: 'text-golden-600 dark:text-golden-400',
  },
  cyan: {
    bg: 'bg-gradient-to-br from-seagreen-500/10 via-seagreen-500/5 to-transparent dark:from-seagreen-500/20 dark:via-seagreen-500/10',
    iconBg: 'bg-gradient-to-br from-seagreen-500 to-twilight-800',
    iconGlow: 'shadow-[0_0_20px_rgba(0,160,160,0.4)]',
    border: 'border-seagreen-500/20',
    text: 'text-seagreen-600 dark:text-seagreen-400',
  },
  teal: {
    bg: 'bg-gradient-to-br from-mint-500/10 via-mint-500/5 to-transparent dark:from-mint-500/20 dark:via-mint-500/10',
    iconBg: 'bg-gradient-to-br from-mint-500 to-seagreen-500',
    iconGlow: 'shadow-[0_0_20px_rgba(0,232,151,0.4)]',
    border: 'border-mint-500/20',
    text: 'text-mint-600 dark:text-mint-400',
  },
};

const TrendIndicator = ({ value, direction }: { value: number; direction: 'up' | 'down' | 'neutral' }) => {
  const colors = {
    up: 'text-mint-700 dark:text-mint-400 bg-mint-100 dark:bg-mint-900/30',
    down: 'text-red-600 dark:text-red-400 bg-red-100 dark:bg-red-900/30',
    neutral: 'text-seagreen-600 dark:text-seagreen-400 bg-seagreen-100 dark:bg-seagreen-900/30',
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
