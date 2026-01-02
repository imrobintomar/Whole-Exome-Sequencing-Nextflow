# Upload Progress Feature

## Overview

Added real-time upload progress indicator to the frontend when users submit FASTQ files for analysis.

## Changes Made

### 1. Backend API (`lib/api.ts`)

Updated the `submitJob` function to support upload progress callbacks:

```typescript
submitJob: async (
  sampleName: string,
  fastqR1: File,
  fastqR2: File,
  onProgress?: (progress: number) => void  // NEW: Optional progress callback
) => {
  // ... FormData setup ...
  
  const response = await api.post<Job>('/jobs/submit', formData, {
    headers: {
      'Content-Type': 'multipart/form-data',
    },
    onUploadProgress: (progressEvent) => {  // NEW: Track upload progress
      if (progressEvent.total) {
        const percentCompleted = Math.round((progressEvent.loaded * 100) / progressEvent.total);
        onProgress?.(percentCompleted);
      }
    },
  });
  return response.data;
}
```

### 2. Upload Form Component (`components/UploadForm.tsx`)

#### New State Variables

- `uploadProgress: number` - Tracks current upload percentage (0-100)

#### New Features

1. **Progress Bar Display**
   - Shows animated progress bar when uploading
   - Displays percentage in real-time
   - Blue-themed progress indicator with spinner

2. **File Size Display**
   - Shows file size for each uploaded FASTQ file
   - Formats bytes to KB/MB/GB automatically

3. **UI Improvements**
   - Disabled state for all inputs during upload
   - Visual feedback with loading spinner
   - Progress percentage in button text

#### Visual Components

**Progress Bar Section:**
```tsx
{uploading && (
  <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
    <div className="flex items-center justify-between mb-2">
      <div className="flex items-center space-x-2">
        <Loader className="h-5 w-5 text-blue-600 animate-spin" />
        <span>Uploading files...</span>
      </div>
      <span>{uploadProgress}%</span>
    </div>
    <div className="w-full bg-blue-200 rounded-full h-2.5">
      <div
        className="bg-blue-600 h-2.5 rounded-full transition-all"
        style={{ width: `${uploadProgress}%` }}
      ></div>
    </div>
    <p className="text-xs text-blue-600 mt-2">
      Please wait while your files are being uploaded...
    </p>
  </div>
)}
```

**Submit Button with Progress:**
```tsx
<button disabled={uploading}>
  {uploading ? (
    <>
      <Loader className="h-5 w-5 animate-spin" />
      <span>Uploading {uploadProgress}%...</span>
    </>
  ) : (
    <span>Submit Analysis</span>
  )}
</button>
```

## User Experience Flow

1. User fills in sample name and selects R1/R2 FASTQ files
2. User sees file names and sizes displayed
3. User clicks "Submit Analysis"
4. **NEW**: Progress bar appears showing upload percentage
5. **NEW**: Button shows "Uploading X%..." with spinner
6. **NEW**: All form inputs are disabled during upload
7. On completion, success message appears
8. User is redirected to "My Jobs" tab

## Technical Details

### Progress Tracking

- Uses Axios `onUploadProgress` event
- Calculates percentage: `(loaded / total) * 100`
- Updates React state in real-time
- Smooth transitions via CSS

### File Size Formatting

New utility function:
```typescript
const formatFileSize = (bytes: number) => {
  if (bytes === 0) return '0 Bytes';
  const k = 1024;
  const sizes = ['Bytes', 'KB', 'MB', 'GB'];
  const i = Math.floor(Math.log(bytes) / Math.log(k));
  return Math.round(bytes / Math.pow(k, i) * 100) / 100 + ' ' + sizes[i];
};
```

## Benefits

1. **Better UX**: Users know exactly how much of the upload is complete
2. **Reduced Anxiety**: Large files (GB size) take time - progress indicator reassures users
3. **Professional Look**: Modern web app standard feature
4. **Error Prevention**: Disabled form prevents accidental double-submissions

## Testing

To test the feature:

1. Start the frontend: `cd frontend && npm run dev`
2. Login to the application
3. Go to "Upload Files" tab
4. Select two FASTQ files (larger files show progress better)
5. Click "Submit Analysis"
6. Observe:
   - Progress bar animating from 0% to 100%
   - Percentage updating in real-time
   - Button showing "Uploading X%..."
   - Form inputs disabled

## Browser Compatibility

Works on all modern browsers that support:
- Axios XMLHttpRequest progress events
- CSS transitions
- Flexbox

## Future Enhancements

Possible improvements:
- Show upload speed (MB/s)
- Show estimated time remaining
- Show individual file progress (R1 vs R2)
- Pause/resume upload capability
- Upload queue for multiple samples

