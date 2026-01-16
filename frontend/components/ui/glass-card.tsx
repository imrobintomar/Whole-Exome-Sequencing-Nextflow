'use client';

import * as React from 'react';
import { cn } from '@/lib/utils';

interface GlassCardProps extends React.HTMLAttributes<HTMLDivElement> {
  variant?: 'default' | 'primary' | 'success' | 'warning' | 'danger' | 'info';
  glow?: boolean;
  hover?: boolean;
  blur?: 'sm' | 'md' | 'lg' | 'xl';
}

const GlassCard = React.forwardRef<HTMLDivElement, GlassCardProps>(
  ({ className, variant = 'default', glow = false, hover = true, blur = 'lg', children, ...props }, ref) => {
    const blurClasses = {
      sm: 'backdrop-blur-sm',
      md: 'backdrop-blur-md',
      lg: 'backdrop-blur-lg',
      xl: 'backdrop-blur-xl',
    };

    const variantClasses = {
      default: 'bg-white/80 dark:bg-slate-900/80 border-white/20 dark:border-white/10',
      primary: 'bg-purple-500/10 dark:bg-purple-500/20 border-purple-500/20',
      success: 'bg-emerald-500/10 dark:bg-emerald-500/20 border-emerald-500/20',
      warning: 'bg-amber-500/10 dark:bg-amber-500/20 border-amber-500/20',
      danger: 'bg-red-500/10 dark:bg-red-500/20 border-red-500/20',
      info: 'bg-cyan-500/10 dark:bg-cyan-500/20 border-cyan-500/20',
    };

    const glowClasses = {
      default: 'shadow-glow-sm',
      primary: 'shadow-glow',
      success: 'shadow-glow-teal',
      warning: 'shadow-[0_0_30px_rgba(245,158,11,0.2)]',
      danger: 'shadow-[0_0_30px_rgba(239,68,68,0.2)]',
      info: 'shadow-glow-cyan',
    };

    return (
      <div
        ref={ref}
        className={cn(
          'rounded-xl border shadow-lg',
          blurClasses[blur],
          variantClasses[variant],
          hover && 'transition-all duration-300 ease-out hover:-translate-y-1 hover:shadow-xl',
          glow && glowClasses[variant],
          className
        )}
        {...props}
      >
        {children}
      </div>
    );
  }
);

GlassCard.displayName = 'GlassCard';

const GlassCardHeader = React.forwardRef<
  HTMLDivElement,
  React.HTMLAttributes<HTMLDivElement>
>(({ className, ...props }, ref) => (
  <div
    ref={ref}
    className={cn('flex flex-col space-y-1.5 p-6', className)}
    {...props}
  />
));
GlassCardHeader.displayName = 'GlassCardHeader';

const GlassCardTitle = React.forwardRef<
  HTMLParagraphElement,
  React.HTMLAttributes<HTMLHeadingElement>
>(({ className, ...props }, ref) => (
  <h3
    ref={ref}
    className={cn('text-xl font-semibold leading-none tracking-tight', className)}
    {...props}
  />
));
GlassCardTitle.displayName = 'GlassCardTitle';

const GlassCardDescription = React.forwardRef<
  HTMLParagraphElement,
  React.HTMLAttributes<HTMLParagraphElement>
>(({ className, ...props }, ref) => (
  <p
    ref={ref}
    className={cn('text-sm text-muted-foreground', className)}
    {...props}
  />
));
GlassCardDescription.displayName = 'GlassCardDescription';

const GlassCardContent = React.forwardRef<
  HTMLDivElement,
  React.HTMLAttributes<HTMLDivElement>
>(({ className, ...props }, ref) => (
  <div ref={ref} className={cn('p-6 pt-0', className)} {...props} />
));
GlassCardContent.displayName = 'GlassCardContent';

const GlassCardFooter = React.forwardRef<
  HTMLDivElement,
  React.HTMLAttributes<HTMLDivElement>
>(({ className, ...props }, ref) => (
  <div
    ref={ref}
    className={cn('flex items-center p-6 pt-0', className)}
    {...props}
  />
));
GlassCardFooter.displayName = 'GlassCardFooter';

export {
  GlassCard,
  GlassCardHeader,
  GlassCardTitle,
  GlassCardDescription,
  GlassCardContent,
  GlassCardFooter,
};
