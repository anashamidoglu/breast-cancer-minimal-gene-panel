# Run locally

Python 3.12 is recommended.

## Results demo

The demo reads public aggregate results already in the repository. No dataset download is needed.

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
pip install -r requirements-demo.txt
streamlit run app.py
```

## Full analysis

```bash
pip install -r requirements.txt
```

See [data/README.md](data/README.md) for preparation and analysis order. Raw data, processed patient-level files, and model artifacts stay local. Analysis stages refuse to overwrite frozen or evaluated outputs; use a separate checkout with a new output directory for a fresh reproduction rather than editing the published results.
