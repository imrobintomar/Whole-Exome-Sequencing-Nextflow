# Whole Exome Sequencing (WES) Platform

**A production-ready, full-stack genomic analysis platform** combining a Nextflow DSL2 bioinformatics pipeline with SaaS-enabled FastAPI backend and modern React frontend for automated whole exome sequencing analysis from raw FASTQ files to annotated, filtered variants.

[![Pipeline](https://img.shields.io/badge/Pipeline-Nextflow%20DSL2-brightgreen)](https://www.nextflow.io/)
[![Backend](https://img.shields.io/badge/Backend-FastAPI-009688)](https://fastapi.tiangolo.com/)
[![Frontend](https://img.shields.io/badge/Frontend-Next.js%2014-000000)](https://nextjs.org/)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

---

## 🌟 Key Features

### Core Analysis
- ✅ **Automated WES Pipeline**: FASTQ → VCF → Annotated TSV (8-stage Nextflow workflow)
- ✅ **Real-time Job Tracking**: Step-by-step pipeline progress monitoring
- ✅ **ACMG Classification**: Automated variant pathogenicity prediction (ACMG/AMP 2015)
- ✅ **Gene Panel Filtering**: PanelApp integration + ACMG Secondary Findings v3.2
- ✅ **Genome Visualization**: IGV.js integration for BAM/VCF viewing
- ✅ **Advanced Analytics**: Variant metrics, chromosome distribution, AF spectrum

### Enterprise SaaS Features
- 🔐 **Firebase Authentication**: Multi-user, OAuth support
- 💳 **Stripe Billing**: Subscription management with usage tracking
- 👨‍💼 **Admin Dashboard**: User/job management & system analytics
- 📊 **Usage Tracking**: Monthly quotas (Free tier: 2 jobs/month)
- 🛡️ **Security**: Audit logging, CORS protection, input validation
- 💬 **Live Chat**: Support system integration (optional)

### DevOps & Infrastructure
- 🚀 **24/7 Deployment**: Systemd service management with Cloudflare Tunnel
- 🌐 **Production Ready**: Cloudflare-secured API (https://api.atgcflow.com/)
- 🐳 **Containerization Ready**: Docker/Singularity support
- 📈 **Dynamic Resources**: Auto-scaling based on system capacity
- 🔄 **Pipeline Control**: Cancel, resume, rerun capabilities
- 📝 **Comprehensive Logging**: Audit trails & error tracking

---

## 📊 Pipeline Overview

```
FASTQ Files (Paired-end)
    ↓
[1] Quality Control & Filtering (fastp)
    ↓
[2] Alignment to hg38 (BWA-MEM)
    ↓
[3] Coordinate Sorting & QC (GATK SortSam + samtools flagstat)
    ↓
[4] Duplicate Marking (GATK MarkDuplicates)
    ↓
[5] Base Quality Score Recalibration (GATK BQSR)
    ↓
[6] Variant Calling (GATK HaplotypeCaller)
    ↓
[7] Multi-source Annotation (1000 Genomes + ANNOVAR)
    ↓
[8] Filtering & UniqueID Assignment
    ↓
Final TSV Output (Chr:Start:Ref:Alt annotated variants)
```

**Typical Runtime:** 24-48 hours per sample on standard workstation

---

## 🏗️ Architecture

### Technology Stack

| Layer | Technology | Purpose |
|-------|-----------|---------|
| **Pipeline** | Nextflow DSL2 | Workflow orchestration |
| **Backend** | FastAPI + Uvicorn | REST API & async processing |
| **Frontend** | Next.js 14 + React 18 | Modern web interface |
| **Database** | SQLite / PostgreSQL | Job tracking & user data |
| **Auth** | Firebase Admin SDK | Authentication & authorization |
| **Payments** | Stripe API | Subscription billing |
| **Tunneling** | Cloudflare Tunnel | Secure public API access |
| **Styling** | Tailwind CSS + Shadcn UI | Component library |
| **Visualization** | Recharts + IGV.js | Charts & genome browser |

### Bioinformatics Tools

| Tool | Version | Purpose |
|------|---------|---------|
| Nextflow | ≥21.04 | Workflow engine |
| BWA | ≥0.7.17 | Read alignment (BWA-MEM) |
| SAMtools | ≥1.17 | BAM manipulation |
| GATK | ≥4.6.0 | Variant calling & processing |
| fastp | ≥0.23.0 | Quality control |
| ANNOVAR | latest | Variant annotation |
| SnpSift/SnpEff | ≥4.3 | Variant filtering |
| bgzip/tabix | any | VCF compression/indexing |

---

## 🚀 Quick Start

### Prerequisites

1. **System Requirements**
   - Linux/macOS (64-bit)
   - 16GB+ RAM (32GB recommended)
   - 200GB+ free disk space
   - Java 1.8+, Python 3.10+

2. **Install Dependencies**

   ```bash
   # Using Conda (recommended)
   conda create -n wes-pipeline -c bioconda \
       nextflow bwa samtools gatk4 fastp bcftools
   conda activate wes-pipeline

   # Or using package manager
   sudo apt-get install bwa samtools bcftools default-jre
   ```

3. **Download Reference Data**
   - Human genome (hg38) with indices
   - GATK known sites VCF files
   - ANNOVAR databases
   - See [Reference Data Requirements](#reference-data-requirements) below

### Installation

```bash
# 1. Clone repository
git clone <repository-url>
cd WholeExome

# 2. Configure backend environment
cd backend
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt

# Create .env file
cat > .env << EOF
SECRET_KEY=your-secret-key-here
DATABASE_URL=sqlite:///./wes_pipeline.db
CORS_ORIGINS=http://localhost:3000
FIREBASE_CREDENTIALS=path/to/firebase-credentials.json
EOF

# 3. Configure frontend
cd ../frontend
npm install

# Create .env.local
cat > .env.local << EOF
NEXT_PUBLIC_API_URL=http://localhost:8000
NEXT_PUBLIC_FIREBASE_API_KEY=your-firebase-api-key
NEXT_PUBLIC_FIREBASE_AUTH_DOMAIN=your-project.firebaseapp.com
EOF

# 4. Update Nextflow configuration
# Edit nextflow.config with your paths (see Configuration section)

# 5. Validate setup
./validate_setup.sh
```

### Running the Platform

#### Option 1: Development Mode

```bash
# Terminal 1: Start backend
cd backend
source venv/bin/activate
python main.py

# Terminal 2: Start frontend
cd frontend
npm run dev

# Visit: http://localhost:3000
```

#### Option 2: Production Mode (24/7 with Cloudflare Tunnel)

```bash
# One-command setup
cd backend
./start-backend-service-improved.sh

# This will:
# ✅ Start backend as systemd service
# ✅ Enable auto-start on boot
# ✅ Enable auto-restart on failure

# Production API is accessible at:
# https://api.atgcflow.com/

# Manage services
./manage-services.sh status    # Check status
./manage-services.sh logs all  # View logs
./manage-services.sh restart   # Restart services
```

#### Option 3: Direct Pipeline Execution (No Web UI)

```bash
# Run pipeline directly
nextflow run main.nf \
  --input_dir /path/to/fastq \
  --output_dir /path/to/results

# With resume capability
nextflow run main.nf -resume

# With reports
nextflow run main.nf -with-trace -with-report -with-timeline
```

---

## ⚙️ Configuration

### 1. Nextflow Configuration (nextflow.config)

Edit the following paths to match your system:

```groovy
params {
    // INPUT/OUTPUT
    input_dir   = '/path/to/fastq/files'
    output_dir  = '/path/to/output/directory'

    // REFERENCE GENOME (hg38)
    reference        = '/path/to/hg38/hg38.fa'
    reference_index  = '/path/to/hg38/hg38.fa.fai'
    reference_dict   = '/path/to/hg38/hg38.dict'
    bwa_index        = '/path/to/hg38/hg38'  # BWA index prefix

    // KNOWN SITES (for GATK BQSR)
    known_sites = [
        '/path/to/Homo_sapiens_assembly38.known_indels.vcf.gz',
        '/path/to/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
    ]

    // ANNOTATION TOOLS
    annovar_dir      = '/path/to/annovar'
    annovar_db       = '/path/to/annovar/humandb'
    annovar_xreffile = '/path/to/annovar/dbNSFP4.7_gene'
    thousand_genomes_dir = '/path/to/1000genomes/vcfs'
    snpsift_jar      = '/path/to/SnpSift.jar'

    // FILTERING THRESHOLDS
    max_af     = 0.05   // Maximum allele frequency (5%)
    min_depth  = 5      // Minimum read depth

    // OPTIONAL: Exome intervals (60x faster than whole genome)
    intervals  = null   // Or: '/path/to/exome_targets.bed'
}

// RESOURCE PROFILES
env {
    GATK_HOME = System.getenv('GATK_HOME') ?: '/usr/local/bin/gatk'
}
```

### 2. Backend Configuration (backend/.env)

```bash
# Security
SECRET_KEY=your-long-random-secret-key
ALGORITHM=HS256

# Database
DATABASE_URL=sqlite:///./wes_pipeline.db
# For production: postgresql://user:password@localhost/wes_db

# Directories
UPLOAD_DIR=./uploads
RESULTS_DIR=./results

# Pipeline
NEXTFLOW_SCRIPT=../main.nf
REFERENCE_GENOME=/path/to/hg38/hg38.fa

# CORS
CORS_ORIGINS=http://localhost:3000,https://your-frontend.vercel.app

# Firebase Authentication
FIREBASE_CREDENTIALS=/path/to/firebase-admin-credentials.json

# Stripe (optional - for billing)
STRIPE_SECRET_KEY=sk_test_...
STRIPE_WEBHOOK_SECRET=whsec_...

# Admin users (Firebase UIDs)
ADMIN_USER_UIDS=uid1,uid2,uid3
```

### 3. Frontend Configuration (frontend/.env.local)

```bash
# API
NEXT_PUBLIC_API_URL=http://localhost:8000
# For production: https://api.atgcflow.com/

# Firebase
NEXT_PUBLIC_FIREBASE_API_KEY=your-api-key
NEXT_PUBLIC_FIREBASE_AUTH_DOMAIN=your-project.firebaseapp.com
NEXT_PUBLIC_FIREBASE_PROJECT_ID=your-project-id
NEXT_PUBLIC_FIREBASE_STORAGE_BUCKET=your-project.appspot.com
NEXT_PUBLIC_FIREBASE_MESSAGING_SENDER_ID=123456789
NEXT_PUBLIC_FIREBASE_APP_ID=1:123456789:web:abcdef
```

### 4. Cloudflare Tunnel Configuration (Optional for custom domains)

If you want to set up your own Cloudflare Tunnel:

```bash
# Install cloudflared
wget https://github.com/cloudflare/cloudflared/releases/latest/download/cloudflared-linux-amd64.deb
sudo dpkg -i cloudflared-linux-amd64.deb

# Authenticate with Cloudflare
cloudflared tunnel login

# Create a tunnel
cloudflared tunnel create wes-backend

# Configure tunnel to point to localhost:8000
# See Cloudflare Tunnel documentation for details
```

---

## 📦 Reference Data Requirements

### 1. Human Reference Genome (hg38)

```bash
# Download from GATK Resource Bundle
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict

# Create BWA index
bwa index -a bwtsw Homo_sapiens_assembly38.fasta
```

### 2. Known Sites VCF Files

```bash
# Download from GATK Resource Bundle
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
```

### 3. ANNOVAR Databases

```bash
# Download ANNOVAR (requires registration)
# Visit: https://annovar.openbioinformatics.org/en/latest/user-guide/download/

# Download required databases
cd annovar
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGeneWithVer humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp42a humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20240416 humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad40_exome humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp150 humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar cosmic84_coding humandb/
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/

# Download dbNSFP gene database
wget https://dbnsfp.s3.amazonaws.com/dbNSFP4.7_gene.complete.gz
gunzip dbNSFP4.7_gene.complete.gz
```

### 4. 1000 Genomes VCF Files (Optional)

Download filtered VCF files for each chromosome (chr1-chr22, chrX) for population frequency annotation.

---

## 📂 Input Data Format

The pipeline expects **paired-end FASTQ files** with this naming convention:

```
SAMPLE1_1.fastq.gz  (Read 1)
SAMPLE1_2.fastq.gz  (Read 2)
SAMPLE2_1.fastq.gz
SAMPLE2_2.fastq.gz
```

**Pattern:** `*_{1,2}.fastq.gz`

Place all FASTQ files in the directory specified by `params.input_dir`.

---

## 📁 Output Structure

```
output_dir/
├── filtered_fastp/              # Quality-controlled FASTQ
│   ├── SAMPLE1_1_filtered.fastq.gz
│   ├── SAMPLE1_2_filtered.fastq.gz
│   └── SAMPLE1.html             # QC report
│
├── Mapsam/                      # BAM files
│   ├── SAMPLE1.sorted.bam
│   ├── SAMPLE1.sorted.bam.bai
│   ├── SAMPLE1_markdup.sorted.bam
│   ├── SAMPLE1_markdup.sorted.bam.bai
│   ├── SAMPLE1_recall.bam       # 🎯 Final recalibrated BAM
│   ├── SAMPLE1_recall.bam.bai
│   ├── SAMPLE1.Stat.txt         # Alignment statistics
│   └── SAMPLE1_recall.table     # BQSR recalibration table
│
├── Germline_VCF/                # Variant files
│   ├── SAMPLE1.vcf.gz           # Raw variants
│   ├── SAMPLE1.vcf.gz.tbi
│   ├── SAMPLE1_1000genome.vcf.gz  # 1000G annotated
│   ├── SAMPLE1.annovar.hg38_multianno.txt  # ANNOVAR annotations
│   └── SAMPLE1_Final_.txt       # 🎯 FINAL FILTERED OUTPUT
│
└── logs/                        # Execution reports
    ├── trace.txt                # Resource usage
    ├── report.html              # Execution summary
    ├── timeline.html            # Timeline visualization
    └── dag.svg                  # Pipeline DAG
```

---

## 📊 Final Output Explained

### `SAMPLE_Final_.txt` (TSV Format)

The final output contains filtered, annotated variants with a **UniqueID** column:

| Column | Description | Example |
|--------|-------------|---------|
| **UniqueID** | Chr:Start:Ref:Alt identifier | chr1:12345:A:G |
| Chr | Chromosome | chr1, chr2, ..., chrX |
| Start | Genomic position (1-based) | 12345678 |
| End | End position | 12345678 |
| Ref | Reference allele | A |
| Alt | Alternate allele | G |
| Func.refGeneWithVer | Functional region | exonic, intronic, UTR3, etc. |
| Gene.refGeneWithVer | Gene name | BRCA1, TP53, etc. |
| GeneDetail.refGeneWithVer | Gene details | - |
| ExonicFunc.refGeneWithVer | Exonic function | missense, nonsense, frameshift, etc. |
| AAChange.refGeneWithVer | Amino acid change | BRCA1:c.5266dupC:p.Gln1756fs |
| AF | Allele frequency | 0.0234 (2.34%) |
| DP | Read depth | 45 |
| QUAL | Variant quality score | 3456.78 |
| SIFT_pred | SIFT prediction | D (deleterious) or T (tolerated) |
| Polyphen2_HDIV_pred | PolyPhen prediction | D (damaging), P (possibly), B (benign) |
| CADD_phred | CADD score | 25.3 (higher = more deleterious) |
| REVEL_score | REVEL score | 0.85 (0-1, higher = more pathogenic) |
| CLNSIG | ClinVar significance | Pathogenic, Benign, VUS, etc. |
| gnomAD_exome_AF | gnomAD allele frequency | 0.000123 |

**Filtering Applied:**
- Allele frequency (AF) ≤ 0.05 (5%)
- Read depth (DP) ≥ 5
- Non-missing values only

---

## 🔌 API Endpoints

### Job Management

```bash
# Submit new job with FASTQ files
POST /jobs/submit
  - Form data: sample_name, fastq_r1, fastq_r2
  - Returns: job_id, status

# Submit with billing enforcement
POST /jobs/submit-with-billing
  - Same as above, but checks subscription limits

# List user's jobs
GET /jobs
  - Returns: array of job objects

# Get job details
GET /jobs/{job_id}
  - Returns: job details, current step, file paths

# Control pipeline
POST /jobs/{job_id}/cancel     # Stop running job
POST /jobs/{job_id}/resume     # Resume failed job
POST /jobs/{job_id}/rerun      # Restart from scratch

# Delete job
DELETE /jobs/{job_id}           # Remove job + files
```

### File Download

```bash
# Download result files
GET /jobs/{job_id}/download/bam              # Final BAM file
GET /jobs/{job_id}/download/bam.bai          # BAM index
GET /jobs/{job_id}/download/raw_vcf          # Raw VCF
GET /jobs/{job_id}/download/raw_vcf.tbi      # VCF index
GET /jobs/{job_id}/download/annotated_vcf    # ANNOVAR TXT
GET /jobs/{job_id}/download/filtered_tsv     # Final TSV

# Download gene-filtered results
POST /jobs/{job_id}/download/filtered
  - Body: { "genes": ["BRCA1", "TP53", ...] }
  - Returns: TSV with only specified genes
```

### Variant Analysis

```bash
# ACMG classification (single variant)
POST /classify/acmg
  - Body: variant object (consequence, gene, AF, scores, etc.)
  - Returns: classification (P/LP/VUS/LB/B) + evidence

# Classify all variants in job
POST /jobs/{job_id}/classify
  - Returns: all variants with ACMG classifications

# Get variant metrics for visualization
GET /jobs/{job_id}/variant-metrics
  - Returns: chromosome distribution, AF spectrum, consequence breakdown
```

### Gene Panels

```bash
# Search PanelApp
GET /panels/search?query=cardiomyopathy
  - Returns: matching panels with IDs

# Get panel genes
GET /panels/{panel_id}/genes?confidence_level=3
  - Returns: gene list from panel

# Get ACMG Secondary Findings (v3.2)
GET /panels/acmg-sf
  - Returns: 81 ACMG-SF genes

# Apply gene panel filter to job
POST /jobs/{job_id}/apply-panel
  - Body: { "genes": ["BRCA1", "TP53", ...] }
  - Returns: filtered variants + statistics
```

### Billing (SaaS Module)

```bash
GET /billing/plans                    # Available subscription plans
GET /billing/subscription             # User's current subscription
POST /billing/subscribe               # Start subscription
POST /billing/checkout                # Create Stripe checkout session
```

### Admin (SaaS Module)

```bash
GET /admin/users                      # List all users
GET /admin/jobs                       # List all jobs
POST /admin/users/{uid}/promote       # Grant admin status
GET /admin/analytics                  # System analytics
```

---

## 🎯 Usage Examples

### Via Web Interface

1. **Login** → Firebase authentication
2. **Dashboard** → Click "New Job"
3. **Upload** → Drag-drop FASTQ files
4. **Submit** → Job starts automatically
5. **Monitor** → Real-time progress tracking
6. **Results** → Download files, apply filters, classify variants

### Via API (cURL)

```bash
# Get authentication token from Firebase first
TOKEN="your-firebase-id-token"

# Submit job
curl -X POST http://localhost:8000/jobs/submit \
  -H "Authorization: Bearer $TOKEN" \
  -F "sample_name=SAMPLE001" \
  -F "fastq_r1=@/path/to/sample_R1.fastq.gz" \
  -F "fastq_r2=@/path/to/sample_R2.fastq.gz"

# Check status
JOB_ID="uuid-from-response"
curl http://localhost:8000/jobs/$JOB_ID \
  -H "Authorization: Bearer $TOKEN"

# Download final TSV
curl http://localhost:8000/jobs/$JOB_ID/download/filtered_tsv \
  -H "Authorization: Bearer $TOKEN" \
  -o results.tsv

# Apply gene panel filter
curl -X POST http://localhost:8000/jobs/$JOB_ID/apply-panel \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"genes": ["BRCA1", "BRCA2", "TP53"]}'
```

### Via Direct Pipeline

```bash
# Run pipeline for multiple samples
nextflow run main.nf \
  --input_dir /data/fastq \
  --output_dir /data/results \
  -profile aggressive \
  -resume

# With execution reports
nextflow run main.nf \
  -with-trace \
  -with-report \
  -with-timeline \
  -with-dag flowchart.svg
```

---

## 🔧 Resource Management

### Profiles

```bash
# Conservative (50% resources - shared systems)
nextflow run main.nf -profile conservative

# Standard (80% resources - default)
nextflow run main.nf -profile standard

# Aggressive (95% resources - dedicated systems)
nextflow run main.nf -profile aggressive
```

### Recommended System Resources

| Component | Minimum | Recommended | Optimal |
|-----------|---------|-------------|---------|
| **CPUs** | 8 cores | 16 cores | 32+ cores |
| **RAM** | 16 GB | 32 GB | 64+ GB |
| **Storage** | 200 GB | 500 GB | 1+ TB SSD |
| **Network** | 10 Mbps | 100 Mbps | 1 Gbps |

### Per-Process Resource Usage

| Step | CPUs | Memory | Time (per sample) |
|------|------|--------|-------------------|
| fastp | 25% | 20% | 30-60 min |
| BWA-MEM | 80% | 75% | 4-8 hours |
| SortSam | 40% | 60% | 1-2 hours |
| MarkDuplicates | 40% | 60% | 2-4 hours |
| BQSR | 85% | 70% | 4-8 hours |
| HaplotypeCaller | 85% | 75% | 6-12 hours |
| ANNOVAR | 80% | 65% | 2-6 hours |

---

## 🛠️ Service Management

### Using manage-services.sh (Recommended)

```bash
# Check status of backend service
./backend/manage-services.sh status

# Start backend service
./backend/manage-services.sh start

# Stop backend service
./backend/manage-services.sh stop

# Restart backend service
./backend/manage-services.sh restart

# View logs
./backend/manage-services.sh logs backend    # Backend logs
./backend/manage-services.sh logs all        # All logs

# Enable/disable auto-start on boot
./backend/manage-services.sh enable
./backend/manage-services.sh disable
```

### Using systemctl Directly

```bash
# Backend service
sudo systemctl status atgcflow-backend.service
sudo systemctl start atgcflow-backend.service
sudo systemctl stop atgcflow-backend.service
sudo systemctl restart atgcflow-backend.service

# View logs
sudo journalctl -u atgcflow-backend.service -f
```

---

## 🔐 Security

### Authentication & Authorization

- **Firebase Admin SDK**: Server-side token verification
- **JWT Sessions**: Secure session management
- **Role-based Access**: Admin vs. regular user permissions
- **Per-user Isolation**: Users can only access their own jobs

### Security Features

- ✅ HTTPS-only communication (ngrok tunneling)
- ✅ CORS policy enforcement
- ✅ Input validation (Pydantic schemas)
- ✅ SQL injection protection (SQLAlchemy ORM)
- ✅ File upload validation (content type + size limits)
- ✅ Audit logging (all user actions)
- ✅ Rate limiting (optional via middleware)
- ✅ IP blocking (security patch available)

### Apply Security Patch

```bash
# Recommended for production deployments
./backend/apply-security-patch.sh

# This adds:
# - Rate limiting (60 requests/min per IP)
# - Attack pattern detection
# - Automatic IP banning
# - Request logging
```

See [backend/SECURITY-PATCH.md](backend/SECURITY-PATCH.md) for details.

---

## 🐛 Troubleshooting

### Pipeline Issues

**Pipeline fails at alignment:**
```bash
# Check BWA index files
ls -lh /path/to/hg38/hg38.{amb,ann,bwt,pac,sa}

# Rebuild if missing
bwa index -a bwtsw hg38.fa
```

**GATK fails with "Index not found":**
```bash
# Pipeline auto-creates indices, but if needed:
samtools index your_file.bam
```

**Out of memory errors (exit code 137):**
```bash
# Use conservative profile
nextflow run main.nf -profile conservative

# Or increase max memory in nextflow.config
```

**SnpSift not found:**
```bash
# Find SnpSift.jar
find /usr /opt -name "SnpSift.jar"

# Update nextflow.config
params.snpsift_jar = '/actual/path/to/SnpSift.jar'
```

### Backend Issues

**Backend won't start:**
```bash
# Check logs
./backend/manage-services.sh logs backend

# Or view directly
tail -f backend/logs/backend.log
```

**Backend connection issues:**
```bash
# Check backend service status
sudo systemctl status atgcflow-backend.service

# View backend logs
./backend/manage-services.sh logs backend

# Restart backend
./backend/manage-services.sh restart
```

**Database errors:**
```bash
# Reinitialize database
cd backend
rm wes_pipeline.db
python -c "from database import init_db; init_db()"
```

**CORS errors in frontend:**
```bash
# Add frontend URL to backend/.env
CORS_ORIGINS=http://localhost:3000,https://your-frontend.vercel.app

# Restart backend
./backend/manage-services.sh restart
```

### Job Failures

**Job stuck in "Running" state:**
```bash
# Check Nextflow process
ps aux | grep nextflow

# View .nextflow.log
cat .nextflow.log

# Resume pipeline
nextflow run main.nf -resume
```

**Job failed - need to resume:**
```bash
# Via API
curl -X POST http://localhost:8000/jobs/{job_id}/resume \
  -H "Authorization: Bearer $TOKEN"

# Or via Web UI: Job Details → Resume
```

---

## 📈 Performance Optimization

### Speed Up Pipeline

1. **Use exome intervals** (60x faster than whole genome):
   ```groovy
   params.intervals = '/path/to/exome_targets.bed'
   ```

2. **Increase parallelism** (process multiple samples):
   ```groovy
   executor {
       queueSize = 12  // Max concurrent processes
   }
   ```

3. **Use faster storage** (SSD/NVMe for work directory)

4. **Optimize resources** (`-profile aggressive` on dedicated systems)

5. **Enable caching** (always use `-resume` when rerunning)

### Reduce Disk Usage

```bash
# Clean work directory after successful run
nextflow clean -f

# Remove intermediate files
rm -rf work/

# Compress old results
tar -czf old_results.tar.gz results/
```

---

## 📚 Documentation

- **[backend/SECURITY-PATCH.md](backend/SECURITY-PATCH.md)** - Security hardening
- **[backend/UPGRADE-NOTES.md](backend/UPGRADE-NOTES.md)** - Migration guide
- **[QUICK-START-24-7.md](QUICK-START-24-7.md)** - 24/7 deployment guide

---

## 🎓 Citation

If you use this pipeline in your research, please cite the underlying tools:

- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. *Nat Biotechnol* 35, 316-319.
- **BWA**: Li, H. and Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. *Bioinformatics* 25, 1754-1760.
- **GATK**: McKenna, A., et al. (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Res* 20, 1297-1303.
- **fastp**: Chen, S., et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34, i884-i890.
- **ANNOVAR**: Wang, K., et al. (2010). ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data. *Nucleic Acids Res* 38, e164.

---

## 🤝 Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

For major changes, please open an issue first to discuss what you would like to change.

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 👨‍💻 Author

**Robin Tomar**
Version: 2.0
Last Updated: 2026-01-12

---

## 🆘 Support

- **Issues**: Report bugs via GitHub Issues
- **Validation**: Run `./validate_setup.sh` before reporting issues
- **Service Status**: `./backend/manage-services.sh status`
- **Logs**: `./backend/manage-services.sh logs all`

---

## 🔗 Links

- **Nextflow**: https://www.nextflow.io/
- **GATK**: https://gatk.broadinstitute.org/
- **ANNOVAR**: https://annovar.openbioinformatics.org/
- **FastAPI**: https://fastapi.tiangolo.com/
- **Next.js**: https://nextjs.org/
- **Firebase**: https://firebase.google.com/
- **Stripe**: https://stripe.com/

---

**🎉 Your comprehensive WES analysis platform is ready!**
