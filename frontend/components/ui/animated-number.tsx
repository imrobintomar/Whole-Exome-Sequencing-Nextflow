'use client';

import * as React from 'react';
import { cn } from '@/lib/utils';

interface AnimatedNumberProps extends React.HTMLAttributes<HTMLSpanElement> {
  value: number;
  duration?: number;
  format?: (n: number) => string;
  prefix?: string;
  suffix?: string;
  decimals?: number;
}

const AnimatedNumber = React.forwardRef<HTMLSpanElement, AnimatedNumberProps>(
  ({ value, duration = 1000, format, prefix = '', suffix = '', decimals = 0, className, ...props }, ref) => {
    const [displayValue, setDisplayValue] = React.useState(0);
    const previousValueRef = React.useRef(0);
    const animationRef = React.useRef<number>();

    React.useEffect(() => {
      const startValue = previousValueRef.current;
      const endValue = value;
      const startTime = performance.now();

      const animate = (currentTime: number) => {
        const elapsed = currentTime - startTime;
        const progress = Math.min(elapsed / duration, 1);

        // Easing function (ease-out)
        const easeOut = 1 - Math.pow(1 - progress, 3);

        const currentValue = startValue + (endValue - startValue) * easeOut;
        setDisplayValue(currentValue);

        if (progress < 1) {
          animationRef.current = requestAnimationFrame(animate);
        } else {
          setDisplayValue(endValue);
          previousValueRef.current = endValue;
        }
      };

      animationRef.current = requestAnimationFrame(animate);

      return () => {
        if (animationRef.current) {
          cancelAnimationFrame(animationRef.current);
        }
      };
    }, [value, duration]);

    const formattedValue = React.useMemo(() => {
      if (format) {
        return format(displayValue);
      }
      return displayValue.toFixed(decimals);
    }, [displayValue, format, decimals]);

    return (
      <span ref={ref} className={cn('tabular-nums', className)} {...props}>
        {prefix}
        {formattedValue}
        {suffix}
      </span>
    );
  }
);

AnimatedNumber.displayName = 'AnimatedNumber';

// Preset formatters
const formatCompact = (n: number): string => {
  if (n >= 1000000) {
    return (n / 1000000).toFixed(1) + 'M';
  }
  if (n >= 1000) {
    return (n / 1000).toFixed(1) + 'K';
  }
  return Math.round(n).toString();
};

const formatWithCommas = (n: number): string => {
  return Math.round(n).toLocaleString();
};

const formatPercent = (n: number): string => {
  return n.toFixed(1) + '%';
};

const formatDuration = (seconds: number): string => {
  if (seconds < 60) {
    return `${Math.round(seconds)}s`;
  }
  if (seconds < 3600) {
    const mins = Math.floor(seconds / 60);
    const secs = Math.round(seconds % 60);
    return `${mins}m ${secs}s`;
  }
  const hours = Math.floor(seconds / 3600);
  const mins = Math.round((seconds % 3600) / 60);
  return `${hours}h ${mins}m`;
};

export { AnimatedNumber, formatCompact, formatWithCommas, formatPercent, formatDuration };
