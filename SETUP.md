# Run locally

Python 3.12 is recommended.

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt
```

See [data/README.md](data/README.md) for preparation and analysis order. The results notebook is in `notebooks/results.ipynb`.

Raw data, processed patient-level files, and model artifacts stay local. Analysis stages protect frozen or evaluated outputs. Use the isolated workspace command in [data/README.md](data/README.md) for a fresh reproduction without changing the published records. The dependency list pins the directly used packages; their transitive dependencies are resolved by pip.
