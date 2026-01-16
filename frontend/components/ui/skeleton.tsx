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
  <div className={cn('rounded-2xl border bg-card p-6', className)}>
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

export {
  Skeleton,
  SkeletonText,
  SkeletonAvatar,
  SkeletonButton,
  SkeletonCard,
  SkeletonTable,
  SkeletonStatsCard,
};
