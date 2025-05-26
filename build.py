#!/usr/bin/env python3
"""
Build script for creating SPZ Python distributions.
"""

import subprocess
import sys
import os

def run_command(cmd):
    """Run a command and check for errors."""
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
        sys.exit(1)
    return result.stdout

def main():
    """Build the package."""
    # Clean previous builds
    if os.path.exists("dist"):
        import shutil
        shutil.rmtree("dist")
    
    if os.path.exists("build"):
        import shutil
        shutil.rmtree("build")
    
    # Build source distribution
    run_command([sys.executable, "-m", "build", "--sdist"])
    
    # Build wheel
    run_command([sys.executable, "-m", "build", "--wheel"])
    
    print("Build completed successfully!")
    print("Distribution files created in dist/")

if __name__ == "__main__":
    main() 