# Short-term heat and stone-related admission records in Brazil and Mexico

Reproducibility materials for the manuscript by Shuo Feng, Chao Shi, Tianming Liu, Haoran Ding, Fenghe Hang, Zhenghao Huang and Xianzhi Meng. Correspondence: mengxianzhi@hrbmu.edu.cn.

Repository: https://github.com/FENG1567/latam-heat-stone-admissions. Release: v1.0.0. See [Code availability](CODE_AVAILABILITY.md).

This repository contains scientific analysis code, aggregate analysis inputs, the checksum-verified accepted numerical freeze, figure/table source data and a manuscript snapshot. It contains no patient-level source files, credentials, or large fitted-model caches.

## Reproduce figures and tables

Use Python 3.12. From the repository root:

```sh
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r environment/requirements-publication.txt
python scripts/verify_repo.py
python code/publication/run_publication.py --frozen data/frozen/frozen_numbers.md --outdir build/publication
```

The runner regenerates Figures 1–4 and S1, Tables 1–4 and supplementary Tables S1–S16. It selects frozen estimates and does **not** refit models. Rendering fonts can change pixels; numerical table values are the comparison target. The authoritative final prose and edited captions are in `paper/`; a generated statistical table is not a newly written manuscript.

## Analysis inputs and model refits

```sh
python scripts/prepare_inputs.py --verify-only
python scripts/analysis.py prepare
python scripts/analysis.py list
python scripts/analysis.py preflight
```

See `docs/REPRODUCIBILITY.md`, the analysis entrypoint help and `config/module_manifest.json` before starting model refits. Supplied inputs restore byte-identically from their manifests; large gzip inputs use binary parts to meet GitHub file limits. Brazil and the historic Mexican comparison are fitted by the core branch; the final Mexican primary is the **complete-source 2015–2022 branch**. Its original runner reconstructs and audits the calendar from the supplied legacy calendar, complete-source/2019 state-day aggregates, clinical state-day totals, 2023 weather and amendment metadata before fitting, even in models or smoke mode. The updated calendar in the package is a reference snapshot, not a way to bypass that reconstruction. Do not substitute the historic Mexican cohort for the final one.

The restored workspace contains working copies: source runners can reserialize calendars there. To repeat byte-identical restoration after a model run, use a new `--workspace` rather than overwriting changed files. Repository input parts remain immutable.

Model refits create new caches; the original multi-gigabyte caches are not needed to inspect frozen outputs and are not distributed here. Full national fresh refitting of every analysis was not repeated during repository preparation. Some source-reconstruction or downstream modules have documented manual/external prerequisites; `preflight` and the module registry identify them. This is not a claim of a tested, one-command raw-download-to-paper pipeline.

## Files and rights

- `code/publication/`: portable frozen-result renderers.
- `analysis/templates/`: preserved scientific source implementations; runtime path adaptation happens in a separate workspace.
- `data/analysis_inputs/`: aggregate inputs, binary parts and input manifest.
- `data/frozen/`: accepted result freeze and adjacent SHA256.
- `data/source_data_index.csv`: maps the paper's 87 neutrally named CSVs to scientific output identities.
- `data/input_headers.json`, `data/source_data_dictionary.json`: input/schema inventory.
- `paper/`: exact final manuscript Markdown, scientific supplements, source CSVs and reference figures.
- `environment/`: recorded Python/R versions and dependency specifications.
- `docs/`: data sources, reproducibility limits and upload instructions.
- `provenance/`: scientific source hashes and tested-package identity.

Self-authored code: MIT. Author-owned summary results and figures: CC BY 4.0. Third-party source records, provider-derived analysis inputs and external geography retain their providers' terms; see `LICENSE-DATA.md` and `THIRD_PARTY_NOTICES.md`. Use `CITATION.cff` to cite the code release. Data and figure licences are described separately in `LICENSE-DATA.md`.

[中文说明](README_中文.md)
