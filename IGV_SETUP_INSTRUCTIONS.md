# IGV.js Browser Setup Instructions

## Issue: "Not authenticated" Error

If you see "Not authenticated" when trying to view the IGV browser, follow these steps:

## 1. Restart the Backend

The backend has been updated to support serving BAM/VCF index files. You need to restart it:

```bash
# Stop the current backend (Ctrl+C in the terminal running it)
# Or kill the process:
pkill -f "python main.py"

# Start the backend again
cd /media/drprabudh/m3/Nextflow-Script/WholeExome/backend
python main.py
```

## 2. Ensure Index Files Exist

IGV.js requires index files for BAM and VCF files:

### For BAM files:
```bash
# Check if .bai index exists
ls -lh /path/to/output/Mapsam/*_recall.bam.bai

# If missing, create it with samtools:
samtools index /path/to/output/Mapsam/sample_recall.bam
```

### For VCF files:
```bash
# Check if .tbi index exists (for bgzipped VCF)
ls -lh /path/to/output/Germline_VCF/*.vcf.gz.tbi

# If VCF is bgzipped (.vcf.gz), create index with tabix:
tabix -p vcf /path/to/output/Germline_VCF/sample.vcf.gz

# If VCF is NOT gzipped, need to bgzip it first:
bgzip /path/to/output/Germline_VCF/sample.vcf
tabix -p vcf /path/to/output/Germline_VCF/sample.vcf.gz
```

## 3. Check Annotated VCF Format

The annotated VCF needs to be bgzipped for IGV.js:

```bash
cd /media/drprabudh/m3/Nextflow-Script/WholeExome

# Find annotated VCF files
find . -name "*annovar*.vcf" -type f

# For each annotated VCF that's not gzipped:
# 1. Bgzip it
bgzip -c sample.annovar.hg38_multianno.vcf > sample.annovar.hg38_multianno.vcf.gz

# 2. Index it
tabix -p vcf sample.annovar.hg38_multianno.vcf.gz
```

## 4. Update Pipeline to Create Indexes

To automatically create indexes for future runs, update your Nextflow pipeline:

### Add to your Nextflow workflow:

```nextflow
// After BAM generation
process INDEX_BAM {
    input:
    path bam

    output:
    path "${bam}.bai"

    script:
    """
    samtools index ${bam}
    """
}

// After VCF annotation
process BGZIP_AND_INDEX_VCF {
    input:
    path vcf

    output:
    tuple path("${vcf}.gz"), path("${vcf}.gz.tbi")

    script:
    """
    bgzip -c ${vcf} > ${vcf}.gz
    tabix -p vcf ${vcf}.gz
    """
}
```

## 5. Backend Changes Made

The download endpoint now supports an `index` query parameter:

- **BAM file**: `GET /jobs/{job_id}/download/bam`
- **BAM index**: `GET /jobs/{job_id}/download/bam?index=true`
- **VCF file**: `GET /jobs/{job_id}/download/annotated_vcf`
- **VCF index**: `GET /jobs/{job_id}/download/annotated_vcf?index=true`

## 6. Frontend Changes Made

Updated `IGVBrowser.tsx` to use `oauthToken` for authentication:

```typescript
{
  name: `${sampleName} - Aligned Reads`,
  type: 'alignment',
  format: 'bam',
  url: `${API_URL}/jobs/${jobId}/download/bam`,
  indexURL: `${API_URL}/jobs/${jobId}/download/bam?index=true`,
  oauthToken: token  // ← Uses Firebase token
}
```

## 7. Testing

After restarting the backend and ensuring index files exist:

1. Login to the app
2. Navigate to Jobs page
3. Click "IGV" button on a completed job
4. The genome browser should load with:
   - BAM track showing aligned reads
   - VCF track showing variants
   - No "Not authenticated" error

## 8. Troubleshooting

### Still seeing "Not authenticated"?

**Check browser console:**
```
F12 → Console tab → Look for error messages
```

**Check backend logs:**
```
# In the terminal running the backend, look for:
# - 401 Unauthorized errors
# - 404 Not Found errors
# - 500 Internal Server errors
```

### Common Issues:

1. **Token expired**: Logout and login again
2. **Files not found**: Check file paths in database match actual files
3. **Index files missing**: Create them with samtools/tabix (see above)
4. **CORS issues**: Ensure backend allows requests from frontend origin

### Verify Authentication:

```bash
# Get token from browser localStorage
# Open browser console and run:
localStorage.getItem('token')

# Test download endpoint with curl:
TOKEN="your-token-here"
curl -H "Authorization: Bearer $TOKEN" \
  http://localhost:8000/jobs/{job-id}/download/bam \
  -I

# Should return 200 OK, not 401 Unauthorized
```

## 9. Alternative: Local File Mode

If authentication continues to be problematic, you can use IGV.js in local file mode:

1. Download BAM and VCF files manually
2. Use IGV desktop app or serve files from local web server
3. This bypasses authentication entirely

## 10. Next Steps

Once IGV is working:

- Navigate to known variant positions
- Inspect read coverage and quality
- Verify variant calls visually
- Use quick gene buttons (BRCA1, BRCA2, TP53, EGFR)

---

## Summary Checklist

- [ ] Restart backend to enable index file support
- [ ] Create .bai indexes for BAM files
- [ ] Create .tbi indexes for VCF files (must be bgzipped)
- [ ] Update pipeline to auto-create indexes
- [ ] Test IGV browser loads without errors
- [ ] Verify BAM and VCF tracks display correctly
