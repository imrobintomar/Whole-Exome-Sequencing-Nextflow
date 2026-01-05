'use client';

import { Badge } from './ui/badge';
import { cn } from '@/lib/utils';

interface ACMGBadgeProps {
  classification: 'Pathogenic' | 'Likely Pathogenic' | 'Uncertain Significance' | 'Likely Benign' | 'Benign';
  className?: string;
  size?: 'sm' | 'md' | 'lg';
}

export default function ACMGBadge({ classification, className, size = 'md' }: ACMGBadgeProps) {
  const getVariant = () => {
    switch (classification) {
      case 'Pathogenic':
        return 'destructive';
      case 'Likely Pathogenic':
        return 'default';
      case 'Uncertain Significance':
        return 'secondary';
      case 'Likely Benign':
        return 'outline';
      case 'Benign':
        return 'outline';
      default:
        return 'secondary';
    }
  };

  const getColors = () => {
    switch (classification) {
      case 'Pathogenic':
        return 'bg-red-600 text-white hover:bg-red-700';
      case 'Likely Pathogenic':
        return 'bg-orange-500 text-white hover:bg-orange-600';
      case 'Uncertain Significance':
        return 'bg-yellow-500 text-black hover:bg-yellow-600';
      case 'Likely Benign':
        return 'bg-green-100 text-green-800 border-green-300 hover:bg-green-200';
      case 'Benign':
        return 'bg-green-600 text-white hover:bg-green-700';
      default:
        return '';
    }
  };

  const getSizeClass = () => {
    switch (size) {
      case 'sm':
        return 'text-xs px-2 py-0.5';
      case 'lg':
        return 'text-base px-4 py-2';
      default:
        return 'text-sm px-3 py-1';
    }
  };

  const getAbbreviation = () => {
    switch (classification) {
      case 'Pathogenic':
        return 'P';
      case 'Likely Pathogenic':
        return 'LP';
      case 'Uncertain Significance':
        return 'VUS';
      case 'Likely Benign':
        return 'LB';
      case 'Benign':
        return 'B';
      default:
        return '?';
    }
  };

  return (
    <Badge
      className={cn(getColors(), getSizeClass(), 'font-semibold', className)}
      title={classification}
    >
      {size === 'sm' ? getAbbreviation() : classification}
    </Badge>
  );
}
