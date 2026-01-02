from pydantic import BaseModel, EmailStr, ConfigDict
from typing import Optional
from datetime import datetime
from database import JobStatus

class UserCreate(BaseModel):
    email: EmailStr
    username: str
    password: str

class UserResponse(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: int
    email: str
    username: str
    created_at: datetime

class Token(BaseModel):
    access_token: str
    token_type: str

class JobCreate(BaseModel):
    sample_name: str

class JobResponse(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: int
    job_id: str
    sample_name: str
    status: JobStatus
    created_at: datetime
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    error_message: Optional[str] = None

class JobDetailResponse(JobResponse):
    bam_path: Optional[str] = None
    raw_vcf_path: Optional[str] = None
    annotated_vcf_path: Optional[str] = None
    filtered_tsv_path: Optional[str] = None
