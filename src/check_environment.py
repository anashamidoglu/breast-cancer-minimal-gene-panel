"""Report dependency versions without downloading data or changing the environment."""

import importlib.metadata
import sys

print(f"Python version: {sys.version.split()[0]}")
print(f"Python location: {sys.executable}")

packages = ["numpy", "pandas", "scikit-learn", "matplotlib", "seaborn", "jupyterlab", "xgboost"]
for package in packages:
    try:
        version = importlib.metadata.version(package)
        print(f"{package}: {version}")
    except importlib.metadata.PackageNotFoundError:
        print(f"{package}: not installed")
