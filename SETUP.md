# Run the project on Windows

Use Python 3.12 for this project. The local environment is in `.venv/` and is excluded from Git: it can be recreated rather than uploaded.

From the project folder in PowerShell, a new user can run:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe src/check_environment.py
```

To open the notebooks:

```powershell
.\.venv\Scripts\python.exe -m jupyterlab
```

Open `notebooks/00_getting_started.ipynb`. Run each code cell with Shift + Enter.

The first notebook uses invented practice data, not patient measurements.

On the original computer, the environment was created using Codex's bundled Python 3.12.14 because the Windows Python shortcuts were unavailable. Other users should install Python 3.12 normally before running the commands above.
