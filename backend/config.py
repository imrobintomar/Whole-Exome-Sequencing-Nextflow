from pydantic_settings import BaseSettings
from typing import List, Optional
import os

class Settings(BaseSettings):
    # Security - REQUIRED from environment
    SECRET_KEY: str  # No default - must be set in .env
    ALGORITHM: str = "HS256"
    ACCESS_TOKEN_EXPIRE_MINUTES: int = 30

    # Directories - Use relative paths or environment variables
    UPLOAD_DIR: str = "./uploads"
    RESULTS_DIR: str = "./results"

    # Nextflow configuration - Can use defaults but should be overridden in .env
    NEXTFLOW_SCRIPT: str = "../main.nf"  # Relative to backend directory
    REFERENCE_GENOME: str = "/media/drprabudh/m1/hg38/hg38.fa"  # Keep your path as default
    KNOWN_SITES_1: str = "/media/drprabudh/m1/vcf_file/Homo_sapiens_assembly38.known_indels.vcf.gz"
    KNOWN_SITES_2: str = "/media/drprabudh/m1/vcf_file/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

    # Database
    DATABASE_URL: str = "sqlite:///./wes_pipeline.db"

    # CORS - Support multiple origins for ngrok and Vercel
    CORS_ORIGINS: str = "http://localhost:3000"

    # SaaS Extension - Optional (graceful degradation if not set)
    FRONTEND_URL: Optional[str] = "http://localhost:3000"
    STRIPE_SECRET_KEY: Optional[str] = None
    STRIPE_PUBLISHABLE_KEY: Optional[str] = None
    STRIPE_WEBHOOK_SECRET: Optional[str] = None
    ADMIN_USER_UIDS: Optional[str] = ""

    class Config:
        env_file = ".env"
        case_sensitive = True

    def get_cors_origins(self) -> List[str]:
        """Parse CORS origins from comma-separated string"""
        return [origin.strip() for origin in self.CORS_ORIGINS.split(",")]

# Initialize settings
try:
    settings = Settings()
    print(f"✅ Configuration loaded successfully")
    print(f"   CORS Origins: {settings.CORS_ORIGINS}")
    print(f"   Upload Dir: {settings.UPLOAD_DIR}")
    print(f"   Results Dir: {settings.RESULTS_DIR}")
except Exception as e:
    print(f"❌ Configuration Error: {e}")
    print("⚠️  Make sure SECRET_KEY is set in .env file")
    print("   Generate one with: python -c 'import secrets; print(secrets.token_urlsafe(32))'")
    raise
