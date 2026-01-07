#!/usr/bin/env python3
"""
Setup script for WES Pipeline Backend
Generates secure SECRET_KEY and creates .env file
"""

import secrets
import os
from pathlib import Path

def generate_secret_key():
    """Generate a secure random secret key"""
    return secrets.token_urlsafe(32)

def setup_env_file():
    """Create or update .env file with required configuration"""
    env_file = Path(__file__).parent / ".env"
    env_example = Path(__file__).parent / ".env.example"

    # Check if .env already exists
    if env_file.exists():
        print(f"⚠️  .env file already exists at {env_file}")
        response = input("Do you want to regenerate SECRET_KEY? (yes/no): ")
        if response.lower() not in ['yes', 'y']:
            print("Keeping existing .env file")
            return

    # Read .env.example as template
    if not env_example.exists():
        print(f"❌ .env.example not found at {env_example}")
        return

    with open(env_example, 'r') as f:
        content = f.read()

    # Generate new SECRET_KEY
    secret_key = generate_secret_key()

    # Replace placeholder SECRET_KEY
    content = content.replace(
        'SECRET_KEY=your-secret-key-change-this-to-something-secure',
        f'SECRET_KEY={secret_key}'
    )

    # Write to .env file
    with open(env_file, 'w') as f:
        f.write(content)

    print(f"✅ Created .env file with secure SECRET_KEY")
    print(f"📄 File location: {env_file}")
    print(f"\n⚠️  IMPORTANT NEXT STEPS:")
    print(f"1. Edit {env_file} and update:")
    print(f"   - CORS_ORIGINS with your ngrok and Vercel URLs")
    print(f"   - REFERENCE_GENOME path if different")
    print(f"   - KNOWN_SITES paths if different")
    print(f"2. If using Firebase, download service account JSON and set path in .env")
    print(f"\n🔐 Your SECRET_KEY: {secret_key}")
    print(f"   (Save this securely! It's also in your .env file)")

def create_directories():
    """Create required directories"""
    backend_dir = Path(__file__).parent

    directories = [
        backend_dir / "uploads",
        backend_dir / "results",
        backend_dir / "data",
    ]

    for directory in directories:
        directory.mkdir(parents=True, exist_ok=True)
        print(f"✅ Created directory: {directory}")

def main():
    print("=" * 60)
    print("  WES Pipeline Backend Setup")
    print("=" * 60)
    print()

    # Create .env file
    setup_env_file()
    print()

    # Create required directories
    print("Creating required directories...")
    create_directories()
    print()

    print("=" * 60)
    print("  Setup Complete! ")
    print("=" * 60)
    print()
    print("Next steps:")
    print("1. Review and update backend/.env file")
    print("2. Install dependencies: pip install -r requirements.txt")
    print("3. Download Firebase service account JSON (if using Firebase)")
    print("4. Start the server: python main.py")
    print()

if __name__ == "__main__":
    main()
