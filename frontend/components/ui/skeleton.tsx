'use client';

import * as React from 'react';
import { cn } from '@/lib/utils';

interface SkeletonProps extends React.HTMLAttributes<HTMLDivElement> {
  variant?: 'text' | 'circular' | 'rectangular' | 'card';
  width?: string | number;
  height?: string | number;
  animation?: 'pulse' | 'shimmer' | 'none';
}

const Skeleton = React.forwardRef<HTMLDivElement, SkeletonProps>(
  ({ className, variant = 'rectangular', width, height, animation = 'shimmer', style, ...props }, ref) => {
    const variantClasses = {
      text: 'rounded h-4',
      circular: 'rounded-full',
      rectangular: 'rounded-lg',
      card: 'rounded-xl',
    };

    const animationClasses = {
      pulse: 'animate-pulse',
      shimmer: 'skeleton',
      none: '',
    };

    return (
      <div
        ref={ref}
        className={cn(
          'bg-muted',
          variantClasses[variant],
          animationClasses[animation],
          className
        )}
        style={{
          width: width,
          height: height,
          ...style,
        }}
        {...props}
      />
    );
  }
);

Skeleton.displayName = 'Skeleton';

// Preset skeleton components for common use cases
const SkeletonText = ({ lines = 3, className }: { lines?: number; className?: string }) => (
  <div className={cn('space-y-2', className)}>
    {Array.from({ length: lines }).map((_, i) => (
      <Skeleton
        key={i}
        variant="text"
        width={i === lines - 1 ? '60%' : '100%'}
        height={16}
      />
    ))}
  </div>
);

const SkeletonAvatar = ({ size = 'md', className }: { size?: 'sm' | 'md' | 'lg'; className?: string }) => {
  const sizeMap = {
    sm: 32,
    md: 40,
    lg: 56,
  };

  return (
    <Skeleton
      variant="circular"
      width={sizeMap[size]}
      height={sizeMap[size]}
      className={className}
    />
  );
};

const SkeletonButton = ({ size = 'md', className }: { size?: 'sm' | 'md' | 'lg'; className?: string }) => {
  const sizeMap = {
    sm: { width: 60, height: 32 },
    md: { width: 80, height: 40 },
    lg: { width: 120, height: 48 },
  };

  return (
    <Skeleton
      variant="rectangular"
      width={sizeMap[size].width}
      height={sizeMap[size].height}
      className={cn('rounded-lg', className)}
    />
  );
};

const SkeletonCard = ({ className }: { className?: string }) => (
  <div className={cn('rounded-xl border bg-card p-6 space-y-4', className)}>
    <div className="flex items-center gap-4">
      <SkeletonAvatar size="md" />
      <div className="flex-1 space-y-2">
        <Skeleton variant="text" width="40%" height={16} />
        <Skeleton variant="text" width="60%" height={12} />
      </div>
    </div>
    <SkeletonText lines={3} />
  </div>
);

const SkeletonTable = ({ rows = 5, columns = 4, className }: { rows?: number; columns?: number; className?: string }) => (
  <div className={cn('space-y-3', className)}>
    {/* Header */}
    <div className="flex gap-4">
      {Array.from({ length: columns }).map((_, i) => (
        <Skeleton key={i} variant="text" height={16} className="flex-1" />
      ))}
    </div>
    {/* Rows */}
    {Array.from({ length: rows }).map((_, rowIndex) => (
      <div key={rowIndex} className="flex gap-4">
        {Array.from({ length: columns }).map((_, colIndex) => (
          <Skeleton key={colIndex} variant="text" height={20} className="flex-1" />
        ))}
      </div>
    ))}
  </div>
);

const SkeletonStatsCard = ({ className }: { className?: string }) => (
  <div
    className={cn('rounded-2xl border bg-card p-6', className)}
    role="status"
    aria-label="Loading statistics"
  >
    <div className="flex items-start justify-between">
      <div className="space-y-3">
        <Skeleton variant="text" width={80} height={14} />
        <Skeleton variant="text" width={100} height={32} />
        <Skeleton variant="text" width={60} height={12} />
      </div>
      <Skeleton variant="rectangular" width={48} height={48} className="rounded-xl" />
    </div>
  </div>
);

// Job list item skeleton
const SkeletonJobItem = ({ className }: { className?: string }) => (
  <div
    className={cn('flex items-center gap-4 p-4 border-b', className)}
    role="status"
    aria-label="Loading job"
  >
    <Skeleton variant="rectangular" width={40} height={40} className="rounded-xl" />
    <div className="flex-1 space-y-2">
      <Skeleton variant="text" width="40%" height={16} />
      <Skeleton variant="text" width="25%" height={12} />
    </div>
    <Skeleton variant="rectangular" width={80} height={24} className="rounded-full" />
  </div>
);

// Job list skeleton (multiple items)
const SkeletonJobList = ({ count = 5, className }: { count?: number; className?: string }) => (
  <div className={cn('space-y-0', className)} role="status" aria-label="Loading job list">
    {Array.from({ length: count }).map((_, i) => (
      <SkeletonJobItem key={i} />
    ))}
    <span className="sr-only">Loading jobs...</span>
  </div>
);

// Analytics chart skeleton
const SkeletonChart = ({ type = 'bar', className }: { type?: 'bar' | 'pie' | 'line'; className?: string }) => (
  <div
    className={cn('rounded-xl border bg-card p-6', className)}
    role="status"
    aria-label="Loading chart"
  >
    <div className="space-y-4">
      <div className="flex justify-between">
        <Skeleton variant="text" width={120} height={20} />
        <Skeleton variant="text" width={80} height={20} />
      </div>
      {type === 'pie' ? (
        <div className="flex justify-center py-8">
          <Skeleton variant="circular" width={160} height={160} />
        </div>
      ) : (
        <div className="flex items-end gap-2 h-48 pt-4">
          {Array.from({ length: 7 }).map((_, i) => (
            <div key={i} className="flex-1 flex flex-col justify-end">
              <Skeleton
                variant="rectangular"
                height={`${Math.random() * 60 + 40}%`}
                className="w-full rounded-t-md"
              />
            </div>
          ))}
        </div>
      )}
    </div>
    <span className="sr-only">Loading chart data...</span>
  </div>
);

// Form skeleton
const SkeletonForm = ({ fields = 3, className }: { fields?: number; className?: string }) => (
  <div
    className={cn('space-y-6', className)}
    role="status"
    aria-label="Loading form"
  >
    {Array.from({ length: fields }).map((_, i) => (
      <div key={i} className="space-y-2">
        <Skeleton variant="text" width={100} height={14} />
        <Skeleton variant="rectangular" width="100%" height={40} className="rounded-lg" />
      </div>
    ))}
    <div className="flex gap-2 pt-4">
      <Skeleton variant="rectangular" width={100} height={40} className="rounded-lg" />
      <Skeleton variant="rectangular" width={80} height={40} className="rounded-lg" />
    </div>
    <span className="sr-only">Loading form...</span>
  </div>
);

// Dashboard overview skeleton
const SkeletonDashboard = ({ className }: { className?: string }) => (
  <div
    className={cn('space-y-8', className)}
    role="status"
    aria-label="Loading dashboard"
  >
    {/* Header section */}
    <div className="rounded-2xl border bg-card p-8">
      <Skeleton variant="text" width={200} height={12} />
      <Skeleton variant="text" width={300} height={32} className="mt-2" />
      <Skeleton variant="text" width={400} height={16} className="mt-2" />
      <div className="flex gap-3 mt-6">
        <Skeleton variant="rectangular" width={140} height={40} className="rounded-xl" />
        <Skeleton variant="rectangular" width={140} height={40} className="rounded-xl" />
      </div>
    </div>

    {/* Stats grid */}
    <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-4">
      {Array.from({ length: 4 }).map((_, i) => (
        <SkeletonStatsCard key={i} />
      ))}
    </div>

    {/* Recent jobs section */}
    <div className="rounded-xl border bg-card">
      <div className="p-6 border-b">
        <Skeleton variant="text" width={120} height={20} />
        <Skeleton variant="text" width={180} height={14} className="mt-1" />
      </div>
      <SkeletonJobList count={5} />
    </div>

    <span className="sr-only">Loading dashboard...</span>
  </div>
);

export {
  Skeleton,
  SkeletonText,
  SkeletonAvatar,
  SkeletonButton,
  SkeletonCard,
  SkeletonTable,
  SkeletonStatsCard,
  SkeletonJobItem,
  SkeletonJobList,
  SkeletonChart,
  SkeletonForm,
  SkeletonDashboard,
};
