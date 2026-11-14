# Analysis runtime

`scripts/analysis.py` uses only the Python standard library. It requires a
POSIX/Linux environment with `Rscript` available on `PATH` (or `--rscript` /
`RSCRIPT` configured). The preserved R modules check and load their own R
packages; the module manifest lists the explicit packages declared by each
runner. No package installation, data download, or model fit occurs during
`list` or `preflight`.

Use the same R package library that was used to validate the preserved scripts.
Set `R_LIBS_USER` in the environment when that library is outside R's default
library paths. The launcher caps common BLAS/OpenMP thread variables at one;
the original runner's `--workers` option remains an explicit module argument.
No R/package lockfile is available in this package, so package versions are
not pinned by the launcher.

From the repository root:

```sh
python3 scripts/analysis.py list
python3 scripts/analysis.py prepare
python3 scripts/analysis.py preflight core_primary_legacy_mx
python3 scripts/analysis.py run core_primary_legacy_mx --workers 1
```

The default isolated workspace is `build/analysis_workspace`; use
`--workspace PATH` or `ANALYSIS_WORKSPACE` to select another location. `prepare`
delegates exact-byte input assembly to `scripts/prepare_inputs.py`, using
`data/analysis_inputs/manifest.json`. `run` copies source templates into the
workspace once and refuses to replace a differing existing script. Use a new
workspace to refresh source copies. Outputs are written only under the selected
workspace.

The two primary Mexico branches are distinct: `core_primary_legacy_mx` uses
`data/mx_calendar.csv.gz`; `final_mexico_complete_source_models` uses the
supplied `data/mx_updated_calendar.csv.gz` and original runner's
`--stage=models`. The wrapper does not silently build or substitute one
calendar for the other. Modules marked `manual` retain a source runner's
authorization or external-QA prerequisite that cannot be proven as an
automated package step.

`preflight` checks the selected source's R syntax, declared R packages, and
required workspace files. It never sources the scripts or runs statistical
models. A successful preflight does not guarantee that inputs satisfy the
original runner's scientific/data validations.
