import { clsx, type ClassValue } from "clsx"
import { twMerge } from "tailwind-merge"

export function cn(...inputs: ClassValue[]) {
  return twMerge(clsx(inputs))
}

export function formatJobId(jobId: string): string {
  if (!jobId) return '';
  // Show first 8 characters of UUID for cleaner display
  return jobId.length > 8 ? jobId.substring(0, 8) : jobId;
}
