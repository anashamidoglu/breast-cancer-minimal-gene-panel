"""Check our tools without downloading data or changing the environment."""

import importlib.metadata
import sys

# Show the actual Python being used, rather than relying on a Windows shortcut.
print(f"Python version: {sys.version.split()[0]}")
print(f"Python location: {sys.executable}")

# A package is an extra set of tools that Python can use.
# Report missing packages so we can plan setup before starting the analysis.
packages = ["numpy", "pandas", "scikit-learn", "matplotlib", "seaborn", "jupyter", "xgboost", "streamlit"]
for package in packages:
    try:
        version = importlib.metadata.version(package)
        print(f"{package}: {version}")
    except importlib.metadata.PackageNotFoundError:
        print(f"{package}: not installed")
