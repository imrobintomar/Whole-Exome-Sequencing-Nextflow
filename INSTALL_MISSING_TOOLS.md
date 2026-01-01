# Installing Missing Tools

## Issue: bgzip and tabix not found

Your system has these tools installed in conda, but they're not in your current PATH.

## Solution Options

### Option 1: Install System-Wide (Recommended for Pipeline)

```bash
# Install htslib package which contains bgzip and tabix
sudo apt-get update
sudo apt-get install -y tabix

# Verify installation
which bgzip
which tabix
bgzip --version
```

### Option 2: Use Conda Environment

You already have bgzip and tabix in your conda environment at:
- `/home/drprabudh/miniconda3/envs/lncrna/bin/bgzip`
- `/home/drprabudh/miniconda3/envs/lncrna/bin/tabix`

**To use them:**

```bash
# Activate the conda environment
conda activate lncrna

# Verify tools are available
which bgzip
which tabix

# Run the validation script again
./validate_setup.sh

# Run the pipeline in this environment
nextflow run main.nf
```

### Option 3: Add Conda Tools to PATH (Quick Fix)

```bash
# Add to current session
export PATH="/home/drprabudh/miniconda3/envs/lncrna/bin:$PATH"

# Verify
which bgzip
which tabix

# Run validation
./validate_setup.sh
```

To make permanent, add to `~/.bashrc`:
```bash
echo 'export PATH="/home/drprabudh/miniconda3/envs/lncrna/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### Option 4: Create Symbolic Links

```bash
# Create symlinks in a directory that's in your PATH
sudo ln -s /home/drprabudh/miniconda3/envs/lncrna/bin/bgzip /usr/local/bin/bgzip
sudo ln -s /home/drprabudh/miniconda3/envs/lncrna/bin/tabix /usr/local/bin/tabix

# Verify
which bgzip
which tabix
```

## Recommended: Option 1 (System Install)

For production pipelines, it's best to have all tools installed system-wide:

```bash
sudo apt-get update
sudo apt-get install -y tabix
```

This package includes both `bgzip` and `tabix`.

## After Installation

Run the validation script again:
```bash
./validate_setup.sh
```

All checks should pass now.

## Alternative: Modify Pipeline to Not Use bgzip/tabix

If you prefer not to install these tools, we can modify the pipeline to work without VCF compression. However, this is **not recommended** because:
- VCF files are large and compression saves significant disk space
- Many tools expect compressed VCF files
- Indexing with tabix enables fast random access

But if you need to proceed without them temporarily, let me know and I can adjust the pipeline.
