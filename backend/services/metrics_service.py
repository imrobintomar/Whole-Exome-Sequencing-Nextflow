"""
System Metrics Service
Collects CPU, memory, disk usage
"""

import psutil
import os
from pathlib import Path

class MetricsService:
    @staticmethod
    def get_system_metrics():
        """
        Get real-time system metrics
        """
        try:
            cpu_percent = psutil.cpu_percent(interval=1)
            memory = psutil.virtual_memory()
            disk = psutil.disk_usage('/')

            return {
                "cpu": {
                    "usage_percent": cpu_percent,
                    "count": psutil.cpu_count()
                },
                "memory": {
                    "total_gb": round(memory.total / (1024**3), 2),
                    "used_gb": round(memory.used / (1024**3), 2),
                    "available_gb": round(memory.available / (1024**3), 2),
                    "percent": memory.percent
                },
                "disk": {
                    "total_gb": round(disk.total / (1024**3), 2),
                    "used_gb": round(disk.used / (1024**3), 2),
                    "free_gb": round(disk.free / (1024**3), 2),
                    "percent": disk.percent
                }
            }
        except Exception as e:
            print(f"Error collecting metrics: {e}")
            # Return default values instead of None
            return {
                "cpu": {
                    "usage_percent": 0,
                    "count": 0
                },
                "memory": {
                    "total_gb": 0,
                    "used_gb": 0,
                    "available_gb": 0,
                    "percent": 0
                },
                "disk": {
                    "total_gb": 0,
                    "used_gb": 0,
                    "free_gb": 0,
                    "percent": 0
                }
            }

    @staticmethod
    def get_nextflow_processes():
        """
        Count active Nextflow processes
        """
        try:
            count = 0
            for proc in psutil.process_iter(['name', 'cmdline']):
                try:
                    if proc.info['name'] == 'nextflow' or \
                       (proc.info['cmdline'] and 'nextflow' in ' '.join(proc.info['cmdline'])):
                        count += 1
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
            return count
        except Exception as e:
            print(f"Error counting Nextflow processes: {e}")
            return 0

    @staticmethod
    def get_storage_usage(base_path: str):
        """
        Get storage usage for uploads/results directories
        """
        try:
            upload_dir = Path(base_path) / "uploads"
            results_dir = Path(base_path) / "results"

            def get_dir_size(path):
                if not path.exists():
                    return 0
                total = 0
                for entry in path.rglob('*'):
                    if entry.is_file():
                        total += entry.stat().st_size
                return total

            upload_size = get_dir_size(upload_dir)
            results_size = get_dir_size(results_dir)

            return {
                "uploads_gb": round(upload_size / (1024**3), 2),
                "results_gb": round(results_size / (1024**3), 2),
                "total_gb": round((upload_size + results_size) / (1024**3), 2)
            }
        except Exception as e:
            print(f"Error calculating storage: {e}")
            return None
