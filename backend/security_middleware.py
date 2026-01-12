"""
Security middleware for protecting against attacks
"""
from fastapi import Request, HTTPException
from fastapi.responses import JSONResponse
from starlette.middleware.base import BaseHTTPMiddleware
from collections import defaultdict
from datetime import datetime, timedelta
import time
import re

class SecurityMiddleware(BaseHTTPMiddleware):
    """
    Comprehensive security middleware that provides:
    - Rate limiting per IP
    - Attack pattern detection
    - Suspicious path blocking
    - Request logging for security analysis
    """

    def __init__(self, app, max_requests_per_minute=60, ban_duration_minutes=15):
        super().__init__(app)
        self.request_counts = defaultdict(list)
        self.banned_ips = {}
        self.max_requests = max_requests_per_minute
        self.ban_duration = timedelta(minutes=ban_duration_minutes)

        # Known attack patterns (IoT/camera exploits, shells, etc.)
        self.attack_patterns = [
            r'/cgi-bin/',
            r'/cgi/',
            r'/shell',
            r'\.asp$',
            r'\.cgi$',
            r'/PSIA/',
            r'/current_config/',
            r'/webs/',
            r'passwd',
            r'\.\./',  # Path traversal
            r'/admin',
            r'/config',
            r'/login\.js',
            r'/device',
            r'cmd=',
            r'get_.*_conf',
            r'/api/v1/system',  # Generic IoT API
        ]
        self.attack_regex = re.compile('|'.join(self.attack_patterns), re.IGNORECASE)

    def is_ip_banned(self, ip: str) -> bool:
        """Check if IP is currently banned"""
        if ip in self.banned_ips:
            ban_time = self.banned_ips[ip]
            if datetime.now() - ban_time < self.ban_duration:
                return True
            else:
                # Ban expired, remove from list
                del self.banned_ips[ip]
        return False

    def ban_ip(self, ip: str):
        """Ban an IP address"""
        self.banned_ips[ip] = datetime.now()
        print(f"🚨 BANNED IP: {ip} for {self.ban_duration.seconds // 60} minutes")

    def is_attack_path(self, path: str) -> bool:
        """Detect if request path matches known attack patterns"""
        return bool(self.attack_regex.search(path))

    def check_rate_limit(self, ip: str) -> bool:
        """Check if IP has exceeded rate limit"""
        now = time.time()
        minute_ago = now - 60

        # Clean old requests
        self.request_counts[ip] = [
            req_time for req_time in self.request_counts[ip]
            if req_time > minute_ago
        ]

        # Add current request
        self.request_counts[ip].append(now)

        # Check limit
        return len(self.request_counts[ip]) <= self.max_requests

    async def dispatch(self, request: Request, call_next):
        """Process each request through security checks"""

        # Get client IP
        client_ip = request.client.host

        # 1. Check if IP is banned
        if self.is_ip_banned(client_ip):
            print(f"🚫 Blocked banned IP: {client_ip} → {request.url.path}")
            return JSONResponse(
                status_code=403,
                content={"detail": "Access denied"}
            )

        # 2. Check for attack patterns
        if self.is_attack_path(request.url.path):
            print(f"🚨 ATTACK DETECTED from {client_ip}:")
            print(f"   Path: {request.url.path}")
            print(f"   Method: {request.method}")
            print(f"   User-Agent: {request.headers.get('user-agent', 'Unknown')}")

            # Ban IP immediately for attack attempts
            self.ban_ip(client_ip)

            return JSONResponse(
                status_code=404,
                content={"detail": "Not found"}
            )

        # 3. Rate limiting
        if not self.check_rate_limit(client_ip):
            print(f"⚠️  Rate limit exceeded: {client_ip}")
            # Ban IP for excessive requests
            self.ban_ip(client_ip)

            return JSONResponse(
                status_code=429,
                content={"detail": "Too many requests"}
            )

        # 4. Log legitimate requests (optional, for debugging)
        if request.url.path not in ["/", "/docs", "/openapi.json"]:
            print(f"✅ {client_ip} → {request.method} {request.url.path}")

        # Process request normally
        response = await call_next(request)
        return response


class IPWhitelistMiddleware(BaseHTTPMiddleware):
    """
    Optional: Only allow specific IPs (for maximum security)
    Use this if you want to restrict access to known IPs only
    """

    def __init__(self, app, allowed_ips=None):
        super().__init__(app)
        self.allowed_ips = set(allowed_ips or [])
        # Add localhost by default
        self.allowed_ips.update(['127.0.0.1', '::1', 'localhost'])

    async def dispatch(self, request: Request, call_next):
        client_ip = request.client.host

        if self.allowed_ips and client_ip not in self.allowed_ips:
            print(f"🚫 Blocked non-whitelisted IP: {client_ip}")
            return JSONResponse(
                status_code=403,
                content={"detail": "Access denied"}
            )

        response = await call_next(request)
        return response


# Utility function to get banned IPs
def get_banned_ips(middleware: SecurityMiddleware) -> dict:
    """Get list of currently banned IPs"""
    now = datetime.now()
    banned = {}
    for ip, ban_time in middleware.banned_ips.items():
        remaining = middleware.ban_duration - (now - ban_time)
        if remaining.total_seconds() > 0:
            banned[ip] = {
                "banned_at": ban_time.isoformat(),
                "remaining_seconds": int(remaining.total_seconds())
            }
    return banned
