import subprocess
import os
import shutil
from datetime import datetime
from pathlib import Path
from typing import Optional
from sqlalchemy.orm import Session
from database import Job, JobStatus
from config import settings

class PipelineRunner:
    def __init__(self, job: Job, db: Session):
        self.job = job
        self.db = db
        self.job_dir = Path(settings.RESULTS_DIR) / job.job_id

    def prepare_job_directory(self):
        """Create job-specific directory structure"""
        self.job_dir.mkdir(parents=True, exist_ok=True)
        (self.job_dir / "input").mkdir(exist_ok=True)
        (self.job_dir / "output").mkdir(exist_ok=True)

        # Copy uploaded FASTQ files to job input directory
        input_dir = self.job_dir / "input"
        shutil.copy(self.job.fastq_r1_path, input_dir / f"{self.job.sample_name}_1.fastq.gz")
        shutil.copy(self.job.fastq_r2_path, input_dir / f"{self.job.sample_name}_2.fastq.gz")

        return input_dir, self.job_dir / "output"

    def run_pipeline(self):
        """Execute Nextflow pipeline"""
        try:
            # Update job status to running
            self.job.status = JobStatus.RUNNING
            self.job.started_at = datetime.utcnow()
            self.db.commit()

            # Prepare directories
            input_dir, output_dir = self.prepare_job_directory()

            # Build Nextflow command
            cmd = [
                "nextflow", "run", settings.NEXTFLOW_SCRIPT,
                "--input_dir", str(input_dir),
                "--output_dir", str(output_dir),
                "--reference", settings.REFERENCE_GENOME,
                "-resume"
            ]

            # Run pipeline
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            stdout, stderr = process.communicate()

            # Check if pipeline completed successfully
            if process.returncode == 0:
                self.update_job_success(output_dir)
            else:
                self.update_job_failed(stderr)

        except Exception as e:
            self.update_job_failed(str(e))

    def update_job_success(self, output_dir: Path):
        """Update job with output file paths"""
        sample = self.job.sample_name

        # Find output files
        bam_file = self.find_file(output_dir / "Mapsam", f"{sample}*_recal.bam")
        raw_vcf = self.find_file(output_dir / "Germline_VCF", f"{sample}*.vcf.gz")
        annotated_vcf = self.find_file(output_dir / "Germline_VCF", f"{sample}*_annovar_annotated.vcf")
        filtered_tsv = self.find_file(output_dir / "Germline_VCF", f"{sample}*_filtered.tsv")

        self.job.status = JobStatus.COMPLETED
        self.job.completed_at = datetime.utcnow()
        self.job.bam_path = str(bam_file) if bam_file else None
        self.job.raw_vcf_path = str(raw_vcf) if raw_vcf else None
        self.job.annotated_vcf_path = str(annotated_vcf) if annotated_vcf else None
        self.job.filtered_tsv_path = str(filtered_tsv) if filtered_tsv else None

        self.db.commit()

    def update_job_failed(self, error_message: str):
        """Update job status to failed"""
        self.job.status = JobStatus.FAILED
        self.job.completed_at = datetime.utcnow()
        self.job.error_message = error_message[:500]
        self.db.commit()

    @staticmethod
    def find_file(directory: Path, pattern: str) -> Optional[Path]:
        """Find file matching pattern in directory"""
        if not directory.exists():
            return None
        files = list(directory.glob(pattern))
        return files[0] if files else None

def run_job_async(job_id: str, db: Session):
    """Run pipeline job asynchronously"""
    job = db.query(Job).filter(Job.job_id == job_id).first()
    if not job:
        return

    runner = PipelineRunner(job, db)
    runner.run_pipeline()
