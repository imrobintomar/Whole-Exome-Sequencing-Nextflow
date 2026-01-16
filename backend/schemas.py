from pydantic import BaseModel, EmailStr, ConfigDict, Field, field_validator
from typing import Optional, List
from datetime import datetime
from database import JobStatus
import re

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
    current_step: Optional[str] = None
    created_at: datetime
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    error_message: Optional[str] = None

class JobDetailResponse(JobResponse):
    bam_path: Optional[str] = None
    raw_vcf_path: Optional[str] = None
    annotated_vcf_path: Optional[str] = None
    filtered_tsv_path: Optional[str] = None
    phenotype_augmented_path: Optional[str] = None
    hpo_terms: Optional[str] = None

class GeneListFilter(BaseModel):
    """Request model for filtering variants by gene list"""
    genes: List[str] = Field(
        ...,
        min_length=1,
        max_length=1000,
        description="List of gene symbols (1-1000 genes)"
    )

    @field_validator('genes')
    @classmethod
    def validate_gene_names(cls, genes: List[str]) -> List[str]:
        """Validate gene names are alphanumeric with allowed characters"""
        validated_genes = []
        gene_pattern = re.compile(r'^[A-Z0-9\-_.]+$', re.IGNORECASE)

        for gene in genes:
            gene = gene.strip()
            if not gene:
                continue
            if not gene_pattern.match(gene):
                raise ValueError(f"Invalid gene name: {gene}. Only alphanumeric characters, hyphens, underscores, and periods allowed.")
            if len(gene) > 50:
                raise ValueError(f"Gene name too long: {gene}. Maximum 50 characters.")
            validated_genes.append(gene.upper())

        if not validated_genes:
            raise ValueError("Gene list cannot be empty after validation")

        # Remove duplicates while preserving order
        seen = set()
        unique_genes = []
        for gene in validated_genes:
            if gene not in seen:
                seen.add(gene)
                unique_genes.append(gene)

        return unique_genes

class HPOPhenotypeInput(BaseModel):
    """Request model for HPO phenotype-driven analysis"""
    hpo_terms: List[str] = Field(
        ...,
        min_length=1,
        max_length=50,
        description="List of HPO term IDs (e.g., HP:0000001)"
    )

    @field_validator('hpo_terms')
    @classmethod
    def validate_hpo_terms(cls, hpo_terms: List[str]) -> List[str]:
        """Validate HPO term IDs"""
        validated_terms = []
        hpo_pattern = re.compile(r'^HP:\d{7}$', re.IGNORECASE)

        for term in hpo_terms:
            term = term.strip().upper()
            if not term:
                continue
            if not hpo_pattern.match(term):
                raise ValueError(f"Invalid HPO term: {term}. Format: HP:0000000")
            validated_terms.append(term)

        if not validated_terms:
            raise ValueError("HPO term list cannot be empty")

        return list(set(validated_terms))
