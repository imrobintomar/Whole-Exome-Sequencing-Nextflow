from fastapi import FastAPI, Depends, HTTPException, status, UploadFile, File, Form, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from sqlalchemy.orm import Session
from pathlib import Path
from contextlib import asynccontextmanager
import uuid
import shutil
import os

from database import get_db, init_db, User, Job, JobStatus
from firebase_auth import get_current_user
from schemas import JobResponse, JobDetailResponse
from config import settings
from pipeline import run_job_async, PipelineRunner
from gene_panels import GenePanelManager
from acmg_classifier import ACMGClassifier
from constraint_data import get_constraint_db
from variant_analyzer import VariantAnalyzer
import threading

# Lifespan event handler
@asynccontextmanager
async def lifespan(app: FastAPI):
    # Startup
    init_db()
    Path(settings.UPLOAD_DIR).mkdir(parents=True, exist_ok=True)
    Path(settings.RESULTS_DIR).mkdir(parents=True, exist_ok=True)
    yield
    # Shutdown (cleanup if needed)

app = FastAPI(title="WES Pipeline API", version="1.0.0", lifespan=lifespan)

# CORS configuration
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.get_cors_origins(),
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Health check endpoint
@app.get("/")
def read_root():
    return {"status": "ok", "message": "WES Pipeline API is running"}

# Submit new job with FASTQ files
@app.post("/jobs/submit", response_model=JobResponse)
async def submit_job(
    background_tasks: BackgroundTasks,
    sample_name: str = Form(...),
    fastq_r1: UploadFile = File(...),
    fastq_r2: UploadFile = File(...),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    # Validate file extensions
    if not fastq_r1.filename.endswith(('.fastq.gz', '.fq.gz')):
        raise HTTPException(status_code=400, detail="R1 file must be FASTQ (.fastq.gz or .fq.gz)")
    if not fastq_r2.filename.endswith(('.fastq.gz', '.fq.gz')):
        raise HTTPException(status_code=400, detail="R2 file must be FASTQ (.fastq.gz or .fq.gz)")

    # Generate unique job ID
    job_id = str(uuid.uuid4())

    # Create upload directory for this job
    upload_dir = Path(settings.UPLOAD_DIR) / job_id
    upload_dir.mkdir(parents=True, exist_ok=True)

    # Save uploaded files
    r1_path = upload_dir / f"{sample_name}_R1.fastq.gz"
    r2_path = upload_dir / f"{sample_name}_R2.fastq.gz"

    with open(r1_path, "wb") as f:
        shutil.copyfileobj(fastq_r1.file, f)
    with open(r2_path, "wb") as f:
        shutil.copyfileobj(fastq_r2.file, f)

    # Create job record in database
    job = Job(
        job_id=job_id,
        user_id=current_user.id,
        sample_name=sample_name,
        status=JobStatus.PENDING,
        fastq_r1_path=str(r1_path),
        fastq_r2_path=str(r2_path)
    )
    db.add(job)
    db.commit()
    db.refresh(job)

    # Start pipeline in background
    background_tasks.add_task(run_job_async, job_id, db)

    return job

# Get all jobs for current user
@app.get("/jobs", response_model=list[JobResponse])
def get_user_jobs(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    jobs = db.query(Job).filter(Job.user_id == current_user.id).order_by(Job.created_at.desc()).all()
    return jobs

# Get job details
@app.get("/jobs/{job_id}", response_model=JobDetailResponse)
def get_job(
    job_id: str,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    return job

# Pipeline Control Endpoints
@app.post("/jobs/{job_id}/cancel")
def cancel_job(
    job_id: str,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Cancel a running pipeline job"""
    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    try:
        runner = PipelineRunner(job, db)
        success = runner.cancel_pipeline()
        return {"message": "Pipeline cancelled successfully" if success else "Pipeline already completed", "job_id": job_id}
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/jobs/{job_id}/rerun")
def rerun_job(
    job_id: str,
    background_tasks: BackgroundTasks,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Rerun a pipeline from scratch"""
    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # Run pipeline in background thread
    def run_rerun():
        runner = PipelineRunner(job, db)
        runner.rerun_pipeline()

    thread = threading.Thread(target=run_rerun)
    thread.daemon = True
    thread.start()

    return {"message": "Pipeline rerun started", "job_id": job_id}

@app.post("/jobs/{job_id}/resume")
def resume_job(
    job_id: str,
    background_tasks: BackgroundTasks,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Resume a failed pipeline"""
    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status != JobStatus.FAILED:
        raise HTTPException(status_code=400, detail=f"Cannot resume job with status: {job.status}")

    # Run pipeline in background thread
    def run_resume():
        runner = PipelineRunner(job, db)
        runner.resume_pipeline()

    thread = threading.Thread(target=run_resume)
    thread.daemon = True
    thread.start()

    return {"message": "Pipeline resume started", "job_id": job_id}

# Download result file
@app.get("/jobs/{job_id}/download/{file_type}")
def download_file(
    job_id: str,
    file_type: str,
    index: bool = False,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status != JobStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Job not completed yet")

    # Map file type to path and custom filename
    file_map = {
        "bam": (job.bam_path, f"{job.sample_name}_recall.bam"),
        "raw_vcf": (job.raw_vcf_path, f"{job.sample_name}.vcf.gz"),
        "annotated_vcf": (job.annotated_vcf_path, f"{job.sample_name}_annotated.vcf"),
        "filtered_tsv": (job.filtered_tsv_path, f"{job.sample_name}_final_annotated.tsv")
    }

    file_info = file_map.get(file_type)
    if not file_info:
        raise HTTPException(status_code=404, detail=f"Invalid file type: {file_type}")

    file_path, download_filename = file_info

    # If requesting index file, append appropriate extension
    if index:
        if file_type == "bam":
            file_path = f"{file_path}.bai"
            download_filename = f"{download_filename}.bai"
        elif file_type in ["raw_vcf", "annotated_vcf"]:
            file_path = f"{file_path}.tbi"
            download_filename = f"{download_filename}.tbi"
        else:
            raise HTTPException(status_code=400, detail=f"Index not available for {file_type}")

    if not file_path or not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail=f"{file_type} {'index' if index else 'file'} not found")

    # Set appropriate media type
    media_type = "application/octet-stream"
    if file_type == "bam" or (index and file_type == "bam"):
        media_type = "application/octet-stream"
    elif file_type in ["raw_vcf", "annotated_vcf"]:
        media_type = "application/gzip" if not index else "application/octet-stream"

    return FileResponse(
        path=file_path,
        filename=download_filename,
        media_type=media_type
    )

# Download filtered variants by gene list
@app.post("/jobs/{job_id}/download/filtered")
async def download_filtered_variants(
    job_id: str,
    genes: list[str],
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Download TSV file filtered to specific genes"""
    import pandas as pd
    import tempfile

    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status != JobStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Job not completed yet")

    if not job.filtered_tsv_path or not os.path.exists(job.filtered_tsv_path):
        raise HTTPException(status_code=404, detail="TSV file not found")

    # Read and filter TSV
    try:
        df = pd.read_csv(job.filtered_tsv_path, sep='\t')

        # Filter by genes (handle multiple genes per row)
        def gene_in_list(gene_str):
            if pd.isna(gene_str):
                return False
            variant_genes = [g.strip() for g in str(gene_str).split(',')]
            return any(g in genes for g in variant_genes)

        gene_column = 'Gene.refGeneWithVer'
        if gene_column in df.columns:
            filtered_df = df[df[gene_column].apply(gene_in_list)]
        else:
            raise HTTPException(status_code=500, detail="Gene column not found in TSV")

        # Write to temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as tmp:
            filtered_df.to_csv(tmp.name, sep='\t', index=False)
            temp_path = tmp.name

        return FileResponse(
            path=temp_path,
            filename=f"{job.sample_name}_filtered_{len(genes)}genes.tsv",
            media_type="text/tab-separated-values",
            background=lambda: os.unlink(temp_path)  # Cleanup after sending
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to filter variants: {str(e)}")

# Gene Panel Endpoints
@app.get("/panels/search")
def search_gene_panels(query: str, current_user: User = Depends(get_current_user)):
    """Search for gene panels from PanelApp"""
    panels = GenePanelManager.search_panels(query)
    return {"results": panels}

@app.get("/panels/{panel_id}/genes")
def get_panel_genes(
    panel_id: int,
    confidence_level: str = "3",
    current_user: User = Depends(get_current_user)
):
    """Get genes from a specific panel"""
    genes = GenePanelManager.get_panel_genes(panel_id, confidence_level)
    return {"panel_id": panel_id, "genes": genes, "count": len(genes)}

@app.get("/panels/acmg-sf")
def get_acmg_secondary_findings(current_user: User = Depends(get_current_user)):
    """Get ACMG Secondary Findings v3.2 gene list"""
    genes = GenePanelManager.get_acmg_secondary_findings_genes()
    return {"genes": genes, "count": len(genes), "version": "3.2"}

# ACMG Classification Endpoint
@app.post("/classify/acmg")
async def classify_variant_acmg(
    variant: dict,
    current_user: User = Depends(get_current_user)
):
    """
    Classify a variant using ACMG/AMP 2015 guidelines

    Request body should contain variant information:
    {
        "consequence": "missense_variant",
        "gene": "BRCA1",
        "af_gnomad": 0.00001,
        "cadd_phred": 25,
        "revel_score": 0.7,
        "sift_pred": "deleterious",
        "polyphen_pred": "probably_damaging",
        "clinvar_sig": "Pathogenic",
        "spliceai_max": 0.1
    }
    """
    try:
        # Add gene constraint data
        gene = variant.get("gene", "")
        if gene:
            constraint_db = get_constraint_db()
            constraint = constraint_db.get_gene_constraint(gene)
            variant["pli"] = constraint.get("pli", 0)
            variant["loeuf"] = constraint.get("loeuf", 1.0)

        # Classify variant
        classifier = ACMGClassifier()
        result = classifier.classify_variant(variant)

        return result

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Classification failed: {str(e)}")

@app.post("/jobs/{job_id}/classify")
async def classify_job_variants(
    job_id: str,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Apply ACMG classification to all variants in a job's TSV file
    """
    import pandas as pd

    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status != JobStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Job not completed yet")

    if not job.filtered_tsv_path or not os.path.exists(job.filtered_tsv_path):
        raise HTTPException(status_code=404, detail="TSV file not found")

    try:
        # Read TSV
        df = pd.read_csv(job.filtered_tsv_path, sep='\t')

        # Get constraint DB
        constraint_db = get_constraint_db()

        # Classify each variant
        classifications = []
        classifier = ACMGClassifier()

        for idx, row in df.iterrows():
            # Helper function to safely get values
            def safe_get(column_name, default=""):
                if column_name in df.columns:
                    val = row[column_name]
                    if pd.isna(val):
                        return default
                    return val
                return default

            def safe_float(column_name, default=0.0):
                if column_name in df.columns:
                    val = row[column_name]
                    if pd.isna(val) or val == "" or val == ".":
                        return default
                    try:
                        return float(val)
                    except (ValueError, TypeError):
                        return default
                return default

            # Map TSV columns to classifier input
            variant = {
                "consequence": safe_get("ExonicFunc.refGeneWithVer", ""),
                "gene": safe_get("Gene.refGeneWithVer", ""),
                "af_gnomad": safe_float("AF", 0.0),
                "cadd_phred": safe_float("CADD_phred", 0.0),
                "revel_score": safe_float("REVEL_score", 0.0),
                "sift_pred": safe_get("SIFT_pred", ""),
                "polyphen_pred": safe_get("Polyphen2_HDIV_pred", ""),
                "clinvar_sig": safe_get("CLNSIG", ""),  # May not exist in TSV
            }

            # Add constraint data
            gene = variant["gene"]
            if gene and gene != ".":
                constraint = constraint_db.get_gene_constraint(gene)
                variant["pli"] = constraint.get("pli", 0)
                variant["loeuf"] = constraint.get("loeuf", 1.0)

            # Classify
            result = classifier.classify_variant(variant)

            classifications.append({
                "position": f"{safe_get('CHROM', '')}:{safe_get('POS', '')}",
                "gene": gene,
                "consequence": variant["consequence"],
                "classification": result["classification"],
                "evidence_summary": result["evidence_summary"],
                "met_criteria": result["met_criteria"]
            })

        return {
            "job_id": job_id,
            "sample_name": job.sample_name,
            "total_variants": len(classifications),
            "classifications": classifications,
            "summary": {
                "pathogenic": sum(1 for c in classifications if c["classification"] == "Pathogenic"),
                "likely_pathogenic": sum(1 for c in classifications if c["classification"] == "Likely Pathogenic"),
                "vus": sum(1 for c in classifications if c["classification"] == "Uncertain Significance"),
                "likely_benign": sum(1 for c in classifications if c["classification"] == "Likely Benign"),
                "benign": sum(1 for c in classifications if c["classification"] == "Benign"),
            }
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Classification failed: {str(e)}")

# Variant Visualization Endpoint
@app.get("/jobs/{job_id}/variant-metrics")
async def get_variant_metrics(
    job_id: str,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get variant analysis metrics for visualization
    """
    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status != JobStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Job not completed yet")

    if not job.filtered_tsv_path or not os.path.exists(job.filtered_tsv_path):
        raise HTTPException(status_code=404, detail="TSV file not found")

    try:
        analyzer = VariantAnalyzer(job.filtered_tsv_path)
        metrics = analyzer.get_all_metrics()

        return {
            "job_id": job_id,
            "sample_name": job.sample_name,
            "metrics": metrics
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Analysis failed: {str(e)}")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
