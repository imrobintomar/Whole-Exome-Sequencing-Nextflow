import { type ClassValue, clsx } from "clsx"
import { twMerge } from "tailwind-merge"

export function cn(...inputs: ClassValue[]) {
  return twMerge(clsx(inputs))
}

/**
 * Format a job UUID into a human-readable job ID
 * Converts: b705070f-6dea-45b0-a570-40312d659a1e
 * To: ATGC2026-001234
 */
export function formatJobId(uuid: string, jobNumber?: number): string {
  // If job number is provided, use it
  if (jobNumber !== undefined) {
    return `ATGC2026-${String(jobNumber).padStart(6, '0')}`;
  }

  // Otherwise, generate a number from the UUID
  // Take first 8 characters of UUID and convert to a number
  const hexPart = uuid.replace(/-/g, '').substring(0, 8);
  const num = parseInt(hexPart, 16) % 999999; // Keep it within 6 digits
  return `ATGC2026-${String(num).padStart(6, '0')}`;
}
