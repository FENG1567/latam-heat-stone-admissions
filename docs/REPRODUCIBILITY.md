# Reproducibility levels and limits

## Frozen-result publication reproduction

The accepted freeze is SHA256 284886236fa12d4ab8f6968bcb7f85244cc2f0337dbf07688e4b9e1181355b47. The publication entrypoint verifies the adjacent checksum before reading its embedded result registry. Numerical cells in all four generated main tables were compared with the final manuscript's tables; the plotting step does not estimate, pool, derive intervals or divide separately pooled RRs. Main Fig2 uses aligned shared-median marginal effects, and Fig4 ends at30°C. The final paper remains authoritative for editorial prose/captions and the strict statement that independent-period refits are not prediction validation.

## Aggregate-input reconstruction

The bundled input manifest identifies compressed aggregate files and exact SHA256 values. `prepare_inputs.py` verifies every part, concatenates original bytes, checks the full file hash and refuses to overwrite a different existing file. It does not deserialize/reserialize CSV values. Large files restore only inside the ignored build workspace. The corrected server Brazil calendar is canonical; the earlier local transfer copy was not used.

## Statistical refits

`analysis/templates` preserves implementations. `analysis.py` creates a separate workspace and adapts historical path prefixes in runtime copies only. See `config/module_manifest.json` for exact commands and dependencies. Model refits are computationally substantial and produce caches excluded from the repository. A successful syntax check or preflight is not proof that every statistical result has been refitted. This packaging run validates source syntax, input checksums/restoration and frozen-result publication output; it does not claim a complete fresh national refit.

Primary branch order: core Brazil plus historic-Mexico comparator; final complete-source Mexican models; covariance/clinical/downstream checks; revision temporal, adult bridge, strict-code, marginal/curve and follow-up branches. The two Mexican calendars and their output identities must not be interchanged.

## Original-source reconstruction

Extraction/retrieval code is retained, but raw SIH files, Mexican discharge-source archives, weather grids and source credentials are external. The municipality crosswalk is included when present in the original archive. Historical Mexican weather request strings are unavailable, so exact original API retrieval cannot be asserted. Original patient/episode identity cannot be recovered from these aggregates. These limits do not prevent rerendering frozen results or restoring the supplied analysis inputs.

## Environment

`environment/requirements-publication.txt` and the matching lock file pin the complete dependency set installed and tested in a clean Python3.12 environment. `environment/R_packages.tsv` and `R_sessionInfo.txt` record the available analysis runtime; these are evidence, not a fabricated full historical renv lock or prebuilt container. Optional acquisition dependencies and credentials are separate. The primary accepted server analyses reported R4.6.1. Use a clean environment, record local versions, keep one native BLAS/OpenMP thread per process and set the number of R workers explicitly.

The Mexican complete-source runner calls its cohort builder in models and smoke stages. All required source aggregates are therefore enumerated in the module registry. It can rewrite compressed working-calendar bytes; that behavior does not alter preserved repository inputs, but it means a used workspace is not a pristine input-restoration target. Use a fresh workspace when repeating restoration.

All released files use `.gitattributes` byte-preservation (`* -text`) so Git line-ending conversion does not invalidate checksum-bound CSVs, source templates or the frozen Markdown.

## Release identity

The public v1.0.0 release is available at https://github.com/FENG1567/latam-heat-stone-admissions. Its current freeze is the checksum stated above and verified by scripts/verify_repo.py. provenance/build_inventory.json and provenance/reproduction_verification.json preserve the earlier packaging-test record; the current release identity is recorded in provenance/release_1.0.0.json. The manuscript snapshots include the subsequent editorial and figure-label updates. No numerical analyses were added by those updates.
