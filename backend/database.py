from sqlalchemy import create_engine, Column, Integer, String, DateTime, Enum as SQLEnum, Boolean, Text
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from datetime import datetime, timezone
import enum
from config import settings

engine = create_engine(
    settings.DATABASE_URL,
    connect_args={"check_same_thread": False} if "sqlite" in settings.DATABASE_URL else {}
)

SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
Base = declarative_base()

class JobStatus(enum.Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"

class User(Base):
    __tablename__ = "users"

    id = Column(Integer, primary_key=True, index=True)
    firebase_uid = Column(String, unique=True, index=True, nullable=False)
    email = Column(String, index=True)
    username = Column(String)
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc))

    # Ban/suspend fields
    is_active = Column(Boolean, default=True)
    is_banned = Column(Boolean, default=False)
    ban_reason = Column(Text, nullable=True)
    banned_at = Column(DateTime, nullable=True)
    banned_by = Column(String, nullable=True)  # Admin UID who banned

class Job(Base):
    __tablename__ = "jobs"

    id = Column(Integer, primary_key=True, index=True)
    job_id = Column(String, unique=True, index=True)
    user_id = Column(Integer, index=True)
    sample_name = Column(String)
    status = Column(SQLEnum(JobStatus), default=JobStatus.PENDING)
    current_step = Column(String, nullable=True)

    fastq_r1_path = Column(String)
    fastq_r2_path = Column(String)

    bam_path = Column(String, nullable=True)
    raw_vcf_path = Column(String, nullable=True)
    annotated_vcf_path = Column(String, nullable=True)
    filtered_tsv_path = Column(String, nullable=True)

    error_message = Column(String, nullable=True)
    process_id = Column(Integer, nullable=True)  # Nextflow process ID for cancellation
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc))
    started_at = Column(DateTime, nullable=True)
    completed_at = Column(DateTime, nullable=True)

def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()

def init_db():
    Base.metadata.create_all(bind=engine)
