# Publication artifact runner

This directory contains a portable, publication-only runner. It reads the accepted `frozen_numbers.md` payload and its adjacent `frozen_numbers.sha256`, then writes figures, tables, source rows, legends and manifests below `--outdir/outputs/`. It does not rebuild the freeze or run statistical estimation.

From a repository root with the frozen file in `00_admin/`:

```sh
python3 code/publication/run_publication.py \
  --frozen 00_admin/frozen_numbers.md \
  --outdir .
```

Python 3, NumPy and Matplotlib are required. The four main figures and Figure S1 are saved as PNG, PDF and SVG. Main tables 1–4 are written as Markdown plus CSV; supplementary tables S1–S16 and their source CSV rows are included. Table 1 retains the full-health-calendar admission-count note.

Fresh renders can differ from archived figures in font pixels because available system fonts vary. Tabular estimates are selected directly from the verified freeze and are not recomputed.


## Current submission layout

The default runner applies the compact figure layout after selecting the same
verified numerical rows. Four main figures and Figure S1 are exported as PNG,
editable PDF/SVG and 600-dpi LZW TIFF. To render only these current figures:

```sh
python3 code/publication/render_readable_figures.py \
  --frozen data/frozen/frozen_numbers.md \
  --outdir build/publication/outputs
```

The numerical table exports remain comprehensive. The manuscript snapshots in
`paper/` contain the compact presentation: Table 3 groups the same sixteen
estimates by grid and rolling window; Table 4 retains all twenty-one effects,
with fit and sample-status detail in Table S8. The supplementary reader and
CSV index in `paper/` provide S1–S17 and links to all 87 complete source files.
The table reader formats values for publication; full-precision values remain
in the CSV files. No model fitting is performed by the display scripts.


## Manuscript text revised on 2026-11-15

The English and Chinese snapshots in `paper/` use the revised manuscript title,
41 references and supplementary table numbers ordered by first citation.
`supplement_table_numbering_20261115.json` maps the earlier labels to the current
ones. Numerical CSV filenames and results are unchanged. The numerical table
renderer retains the original source-table order; the manuscript reader and
CSV index supply the current narrative order. No statistical models were refitted
for this text revision.


## Figure typography and layout (2026-11-15)

The figure renderer now uses the companion `layout_result_panels.py` module.
Multipart figure labels are uppercase A/B/C at 18 pt in the native exported
PDF; text and numeric labels are bold. Figure 3 places the sensitivity grid
on the left and country estimates and period refits on the right. Figure 4
uses a separate shared legend above the two curve panels and a support panel
below. Source selections and numerical values are unchanged.

For exact Arial Bold typography, provide a locally licensed font file using
`LATAM_ARIAL_BOLD=/absolute/path/to/arialbd.ttf` when running the command above.
The font file is not distributed with this repository. Without that setting,
the renderer uses installed Arial when available and a sans-serif fallback
otherwise. Font appearance may therefore differ between systems. PNG exports
are 300 dpi; TIFF exports are 600 dpi with LZW compression. PDF and SVG retain
editable text. The standalone result-panel module relies on the parent
renderer for shared styling and exports; execute the parent renderer.
