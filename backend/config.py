from pydantic_settings import BaseSettings
from typing import List

class Settings(BaseSettings):
    SECRET_KEY: str = "your-secret-key-change-this"
    ALGORITHM: str = "HS256"
    ACCESS_TOKEN_EXPIRE_MINUTES: int = 30

    UPLOAD_DIR: str = "/media/drprabudh/m3/Nextflow-Script/WholeExome/uploads"
    RESULTS_DIR: str = "/media/drprabudh/m3/Nextflow-Script/WholeExome/results"

    NEXTFLOW_SCRIPT: str = "/media/drprabudh/m3/Nextflow-Script/WholeExome/main.nf"
    REFERENCE_GENOME: str = "/media/drprabudh/m1/hg38/hg38.fa"
    KNOWN_SITES_1: str = "/media/drprabudh/m1/vcf_file/Homo_sapiens_assembly38.known_indels.vcf.gz"
    KNOWN_SITES_2: str = "/media/drprabudh/m1/vcf_file/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

    DATABASE_URL: str = "sqlite:///./wes_pipeline.db"
    CORS_ORIGINS: str = "http://localhost:3000"

    class Config:
        env_file = ".env"

    def get_cors_origins(self) -> List[str]:
        return [origin.strip() for origin in self.CORS_ORIGINS.split(",")]

settings = Settings()
