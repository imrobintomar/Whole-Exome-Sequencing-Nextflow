import subprocess
import os
import shutil
import signal
from datetime import datetime, timezone
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

    def update_pipeline_step(self, step: str):
        """Update the current pipeline step in the database with error handling"""
        try:
            self.job.current_step = step
            self.db.commit()
        except Exception as e:
            print(f"Warning: Failed to update pipeline step to '{step}': {e}")
            # Try to rollback and continue
            try:
                self.db.rollback()
            except:
                pass

    def parse_nextflow_output(self, line: str) -> Optional[str]:
        """Parse Nextflow output to identify current step"""
        step_keywords = {
            'fastpQC': 'Quality Control',
            'bwaMem': 'Alignment',
            'sortSam': 'Sorting',
            'flagstat': 'Quality Statistics',
            'markDuplicates': 'Deduplication',
            'sortSamPostDedup': 'Post-dedup Sorting',
            'baseRecalibrator': 'Base Recalibration',
            'applyBQSR': 'Applying BQSR',
            'haplotypeCaller': 'Variant Calling',
            'snpsiftAnnotate1000G': '1000 Genomes Annotation',
            'annovarAnnotate': 'ANNOVAR Annotation',
            'extractFilterFields': 'Filtering & Extraction'
        }

        for keyword, step_name in step_keywords.items():
            if keyword in line:
                return step_name
        return None

    def run_pipeline(self):
        """Execute Nextflow pipeline"""
        try:
            # Update job status to running
            self.job.status = JobStatus.RUNNING
            self.job.started_at = datetime.now(timezone.utc)
            self.job.current_step = "Initializing"
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

            # Run pipeline with real-time output monitoring
            # Use process group for proper cancellation
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                preexec_fn=os.setsid  # Create new process group
            )

            # Store process ID for cancellation
            self.job.process_id = process.pid
            self.db.commit()

            output_lines = []
            # Monitor output in real-time with error handling
            try:
                for line in iter(process.stdout.readline, ''):
                    if line:
                        output_lines.append(line)
                        # Try to detect current step from output
                        try:
                            step = self.parse_nextflow_output(line)
                            if step:
                                self.update_pipeline_step(step)
                        except Exception as parse_error:
                            print(f"Warning: Error parsing pipeline output: {parse_error}")
                            # Continue monitoring even if parsing fails
            except Exception as monitor_error:
                print(f"Error monitoring pipeline output: {monitor_error}")
                # Continue to wait for process to complete

            process.wait()
            full_output = ''.join(output_lines)

            # Check if pipeline completed successfully
            if process.returncode == 0:
                self.update_job_success(output_dir)
            else:
                self.update_job_failed(full_output)

        except Exception as e:
            self.update_job_failed(str(e))

    def update_job_success(self, output_dir: Path):
        """Update job with output file paths"""
        sample = self.job.sample_name

        # Find output files
        bam_file = self.find_file(output_dir / "Mapsam", f"{sample}_recall.bam")
        raw_vcf = self.find_file(output_dir / "Germline_VCF", f"{sample}.vcf.gz")
        annotated_vcf = self.find_file(output_dir / "Germline_VCF", f"{sample}.annovar.hg38_multianno.vcf")
        filtered_tsv = self.find_file(output_dir / "Germline_VCF", f"{sample}_final_annotated.tsv")

        self.job.status = JobStatus.COMPLETED
        self.job.current_step = "Completed"
        self.job.completed_at = datetime.now(timezone.utc)
        self.job.process_id = None  # Clear process ID
        self.job.bam_path = str(bam_file) if bam_file else None
        self.job.raw_vcf_path = str(raw_vcf) if raw_vcf else None
        self.job.annotated_vcf_path = str(annotated_vcf) if annotated_vcf else None
        self.job.filtered_tsv_path = str(filtered_tsv) if filtered_tsv else None

        self.db.commit()

    def update_job_failed(self, error_message: str):
        """Update job status to failed"""
        self.job.status = JobStatus.FAILED
        self.job.completed_at = datetime.now(timezone.utc)
        self.job.process_id = None  # Clear process ID
        self.job.error_message = error_message[:500]
        # Keep current_step if set, otherwise set to "Initializing"
        if not self.job.current_step:
            self.job.current_step = "Initializing"
        self.db.commit()

    @staticmethod
    def find_file(directory: Path, pattern: str) -> Optional[Path]:
        """Find file matching pattern in directory"""
        if not directory.exists():
            return None
        files = list(directory.glob(pattern))
        return files[0] if files else None

    def cancel_pipeline(self):
        """Cancel a running pipeline"""
        if not self.job.process_id:
            raise ValueError("No process ID found for this job")

        if self.job.status not in [JobStatus.RUNNING, JobStatus.PENDING]:
            raise ValueError(f"Cannot cancel job with status: {self.job.status}")

        try:
            # Kill the process and its children
            os.killpg(os.getpgid(self.job.process_id), signal.SIGTERM)

            # Update job status
            self.job.status = JobStatus.FAILED
            self.job.error_message = "Pipeline cancelled by user"
            self.job.completed_at = datetime.now(timezone.utc)
            self.job.process_id = None
            self.db.commit()

            return True
        except ProcessLookupError:
            # Process already finished
            self.job.process_id = None
            self.db.commit()
            return False
        except Exception as e:
            raise Exception(f"Failed to cancel pipeline: {str(e)}")

    def rerun_pipeline(self):
        """Rerun a pipeline from scratch"""
        # Reset job status
        self.job.status = JobStatus.PENDING
        self.job.current_step = None
        self.job.started_at = None
        self.job.completed_at = None
        self.job.error_message = None
        self.job.process_id = None
        self.job.bam_path = None
        self.job.raw_vcf_path = None
        self.job.annotated_vcf_path = None
        self.job.filtered_tsv_path = None
        self.db.commit()

        # Clean up old output directory
        if self.job_dir.exists():
            shutil.rmtree(self.job_dir)

        # Run pipeline
        self.run_pipeline()

    def resume_pipeline(self):
        """Resume a failed pipeline using Nextflow -resume"""
        if self.job.status not in [JobStatus.FAILED]:
            raise ValueError(f"Can only resume failed jobs, current status: {self.job.status}")

        # Update status to running
        self.job.status = JobStatus.RUNNING
        self.job.started_at = datetime.now(timezone.utc)
        self.job.error_message = None
        self.db.commit()

        # Run pipeline (already has -resume flag in run_pipeline)
        self.run_pipeline()

def run_job_async(job_id: str):
    """Run pipeline job asynchronously"""
    # Create a new database session for this background task
    from database import SessionLocal
    db = SessionLocal()

    try:
        job = db.query(Job).filter(Job.job_id == job_id).first()
        if not job:
            return

        runner = PipelineRunner(job, db)
        runner.run_pipeline()
    finally:
        db.close()
