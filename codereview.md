# Code Review Report — `spatial_tls_manuscript` (publication readiness)

**Date:** 2026-06-26

## Context

Goal: get the repo into a publication-ready state for a public GitHub release tied to
Bangs et al., 2025. This is a **report only** — no edits and no git actions are performed.
You will implement the changes yourself (and may later ask me to make specific edits).

**Method:** read-only audit (Explore agents + targeted file reads). The built-in
`code-review`/`security-review` skills were intentionally *not* used as the primary tool —
they review a diff, which is the wrong fit for a whole-repo readiness audit. Run
`security-review` at the end on your resulting diff if you want a final gate.

**Scope.** Reviewed: root docs (`README.md`, `CLAUDE.md`, `LICENSE.txt`, `.gitignore`),
`config/paths.py` + `.R`, `conda_envs/`, `inputs/`, `utils/*.py`, `initial_processing/00–08`,
and `downstream_analysis/` root notebooks (`00_sq_interactions`, `01_DE`,
`02_distance_tls_bronchi`, `03_distance_tls_bronchi_vis`, `04_distance_tls_radial_vis`,
`save_final_adata`).

**Explicitly excluded** (per your instruction): `downstream_analysis/prev_organization`,
`downstream_analysis/prev_organization_2`, `downstream_analysis/cellchat` (entire folder),
`initial_processing/prev_organization`.

Findings are tiered **P0 (must-fix before public)** → **P3 (polish)**. Each has *what / where / fix*.

---

## P0 — Must fix before going public

### P0.1 — Private GCP bucket ID leaked in an active notebook
- **Where:** `downstream_analysis/save_final_adata.ipynb` contains
  `gsutil -m cp ... gs://wftissueimmaifi-219b3400-2a33-483a-a457-d5272e41d7db`.
  (Also present in excluded `prev_organization/` copies and in `.ipynb_checkpoints/` — see P1.1.)
- **Why:** exposes a unique internal cloud bucket identifier.
- **Fix:** remove the `gsutil`/`gs://` line(s), or replace the bucket with a placeholder
  (e.g. `gs://<YOUR_BUCKET>`). Verify with `git grep -n 'wftissueimmaifi'` → no hits in
  files you intend to publish.
- **Note (your call):** the bucket ID also lives in **git history**. You said no git actions,
  so flagging only: a clean public release would require history rewrite (BFG/`filter-repo`)
  or publishing from a fresh commit. Implement however you prefer — not done here.

### P0.2 — 79 tracked `.ipynb_checkpoints/` files
- **Where:** 79 files across the repo (root, `config/`, `conda_envs/`, `utils/`,
  `initial_processing/`, `downstream_analysis/`). `.gitignore` already lists
  `.ipynb_checkpoints/`, but these were committed earlier so they remain tracked.
- **Why:** dev artifacts; they also re-introduce the leaks above and stale path strings.
- **Fix:** untrack them (`git rm -r --cached **/.ipynb_checkpoints` — your action, not mine).
  Verify with `git ls-files | grep -c '\.ipynb_checkpoints/'` → expect 0.

---

## P1 — Reproducibility / portability

### P1.1 — Hardcoded `/home/workspace/...` paths
- **Where (canonical source, the only spot that *needs* editing per design):**
  `config/paths.py` and `config/paths.R` — `DATA_DIR`, `BASE_OUTDIR`, `REF_DIR`,
  `INPUTS_DIR`, `FUNCTIONS_DIR` all hardcode `/home/workspace/...`.
- **Also in:** `README.md` run instructions; many notebook cells/outputs (148 tracked files
  contain `/home/workspace`, mostly checkpoints/outputs).
- **Fix options (recommend a + b):**
  a. Refactor `config/paths.*` to derive from the repo root + env var override, e.g.
     `BASE_DIR = Path(os.environ.get("TLS_DATA_DIR", Path(__file__).resolve().parents[1]))`,
     and build the rest relative to it. Keep `.py` and `.R` in sync (CLAUDE.md convention).
  b. Make `INPUTS_DIR`/`FUNCTIONS_DIR` relative to the repo automatically (they live *inside*
     the repo, so they should never be absolute): `Path(__file__).resolve().parents[1] / "inputs"`.
  c. Embedded paths inside notebook **cells/outputs** are cosmetic once `config/paths.*` is
     fixed — addressed by the output-stripping decision in P3.3.

### P1.2 — README env-setup commands are broken / machine-specific
- **Where:** `README.md` lines 44–48.
- **Issues:** (1) all three `conda env create` commands hardcode
  `--prefix /home/workspace/environment/space2`; (2) the **`r_seurat`** and **`cellcharter2`**
  commands both point `--prefix` at `.../space2` (copy-paste bug — they'd collide).
- **Fix:** drop the absolute `--prefix`, or document it as a user-edited placeholder; give each
  env its own prefix/name. Example: `conda env create -f conda_envs/space2_20250604.yml`.

### P1.3 — README notebook→env table references notebooks that no longer exist
- **Where:** `README.md` lines 59–60 list `downstream_analysis/00_downstream_analysis.ipynb`
  and `01_downstream_analysis_hypothesis-testing.ipynb` — those names are gone; current files
  are `00_sq_interactions`, `01_DE`, `02_distance_tls_bronchi`, `03_distance_tls_bronchi_vis`,
  `04_distance_tls_radial_vis`, `save_final_adata`.
- **Fix:** rebuild the table to match the actual `initial_processing/00–08` and
  `downstream_analysis/` notebooks, with the correct env per notebook.

### P1.4 — Data-directory name mismatch between README and config
- **Where:** `README.md` (lines 27–37) shows a `spatial_tls_manuscript_data/` tree; but
  `config/paths.py`/`.R` point `DATA_DIR` at `/home/workspace/data/temp/mouse_lung` and
  `REF_DIR` at `.../hurskainen_ref`.
- **Fix:** make the README tree and the config defaults describe the *same* layout so a new
  user can place data once and run.

### P1.5 — `conda_envs/` not fully documented
- **Where:** 5 YAMLs exist — `space2_20250604`, `cellcharter2_20250604`, `r_seurat_20250604`,
  `cellchat`, `space2_light`. README documents only 3; `cellchat` and `space2_light` are
  unmentioned.
- **Fix:** either document all env files (what each is for) or remove the ones not needed by
  the in-scope, published pipeline. (Note: `cellchat.yml` supports the excluded `cellchat/`
  folder — decide whether it ships.)

---

## P2 — `utils/` code quality (publication-grade correctness)

### P2.1 — `de_utils.py` silences **all** warnings globally at import time
- **Where:** `utils/de_utils.py:47–48` — `import warnings; warnings.simplefilter(action='ignore', category=Warning)`
  runs on import, so *any* notebook that imports `de_utils` loses all warnings repo-wide.
  `warnings` is also imported twice (lines 5 and 47).
- **Fix:** remove the global `simplefilter`, or scope it with a `warnings.catch_warnings()`
  context inside the specific function that needs it. De-dupe the import.

### P2.2 — `plotting_utils.py` uses `colorcet` but never imports it
- **Where:** `utils/plotting_utils.py:117` — `sns.color_palette(colorcet.glasbey_dark, ...)`
  in `zone_composition`'s default-palette branch, but `colorcet` is not imported (top of file
  imports math/anndata/squidpy/pandas/scanpy/scvi/numpy/matplotlib/seaborn only).
- **Why:** calling `zone_composition` without `palette_map` raises `NameError`. Latent bug.
- **Fix:** add `import colorcet` (and to the env spec), or replace with an already-imported
  palette source.

### P2.3 — Unbound-variable branches (latent `NameError`s)
- **Where:**
  - `plotting_utils.py:81–82` — final `else` of the legend logic references `handles, labels`,
    which are only assigned in the `legend_outside` branch.
  - `plotting_utils.py:204–207` (`zone_composition_heatmap_2`) — `summary_frac_scaled` is only
    set when `scale_by` is `'rows'` or `'cols'`; any other value → `NameError`.
- **Fix:** initialize/guard these variables, or `raise ValueError` for invalid `scale_by`.

### P2.4 — Large blocks of unused heavy imports
- **Where:** `utils/de_utils.py:1–48` imports many modules the file never uses —
  `os, time, re, math, scipy.spatial.distance_matrix/cdist, scipy.stats, scipy.ndimage,
  gaussian_filter, find_peaks, shapely (Point, MultiPolygon), alphashape, optimizealpha,
  geopandas, sklearn (MinMaxScaler, PowerTransformer), pickle, joblib, random, mcolors`.
  `plotting_utils.py` similarly imports `scvi`, `squidpy`, `anndata` but doesn't use them.
- **Why:** slow imports (scvi/geopandas/alphashape are heavy) and inflated apparent
  dependencies; reads as in-progress code.
- **Fix:** prune to what each module actually uses. `de_utils` effectively needs
  `numpy, pandas, matplotlib, scanpy, adjustText`.

### P2.5 — Smaller correctness / cleanliness items
- **`de_utils.py:357–358`** — `dot_save_path` save is at function scope, not inside the
  `if plot and gene_order_filt:` block, so it can save an empty/stale figure when no dotplot
  was drawn. Move the save under the plot block.
- **`de_utils.py:50`** — bare `sys.version_info` statement does nothing (leftover).
- **`de_utils.py:105`** — `from adjustText import adjust_text` mid-file; move to top with
  other imports.
- **`distance_utils.py:250`** — bare `except:` in `scatter_with_gaussian_kde`; narrow to the
  expected exception.
- **Docstrings/naming:** `zone_composition_heatmap_2` (plotting_utils) lacks a docstring and
  the `_2` suffix reads as a dev iteration — consider renaming/documenting. A few commented-out
  dead lines in `distance_utils.py` (e.g. 218, 345, 444) can be removed.

---

## P3 — Hygiene & polish

### P3.1 — `.gitignore` gaps
- **Where:** current `.gitignore` lists `.ipynb_checkpoints/ *.pdf *.h5ad *.parquet *.pt
  *.json *.geojson *.gz`.
- **Issues:** no Python artifacts (`__pycache__/`, `*.pyc`, `.ipynb_checkpoints` is there but
  `*.egg-info/`, `.venv/` absent). Also `*.pkl`/`*.joblib` were removed from ignore in a recent
  commit — confirm that's intended given `inputs/` legitimately ships pickled palettes (those
  are committed deliberately, so blanket-ignoring `*.pkl` would need a `!inputs/*.pkl` exception).
- **Fix:** add `__pycache__/`, `*.pyc`, `.venv/`, `*.egg-info/`; decide pickle policy explicitly.

### P3.2 — Add a `CITATION.cff`
- **Where:** none present. README cites Bangs et al., 2025 and Zenodo DOI
  `10.5281/zenodo.15858877`.
- **Fix:** add `CITATION.cff` with authors, title, and the Zenodo DOI so GitHub shows a
  "Cite this repository" button.

### P3.3 — Notebook cell outputs embed paths and bloat diffs (your decision)
- **Where:** in-scope notebooks carry executed outputs containing `/home/workspace` strings.
- **Options (pick at implementation):**
  - **Strip all outputs** (nbstripout / "Clear All Outputs") — cleanest, smallest repo,
    removes embedded paths; readers re-run to see results.
  - **Keep outputs** for illustration — then scrub path/bucket strings in the leaking cells.
  - **Strip only leaking cells** — keep result plots/tables, clear the few cells printing paths.
- **Recommendation:** strip all outputs for a clean public release; results are reproducible
  from the data + notebooks.

### P3.4 — `config/paths.*` typo
- **Where:** `config/paths.py:18` and `config/paths.R:10` — comment reads
  "Lodal directory with scRNAseq ref" → "Local".

### P3.5 — Update `CLAUDE.md` to match reality
- **Where:** CLAUDE.md describes a `02b_distance...`/`03_distance_tls_radial` downstream
  naming that no longer matches the files (you confirmed CLAUDE.md is the outdated side here).
- **Fix:** update the `downstream_analysis/` description to the current
  `02 / 03_..._vis / 04_..._vis` names. (CLAUDE.md is a dev doc — optional for the public repo,
  but worth keeping accurate if it ships.)

---

## Suggested implementation order

1. **P0** leaks + checkpoint untracking (blockers).
2. **P1.1** refactor `config/paths.*` to relative/env-var; then **P1.2–P1.5** README/env fixes.
3. **P2** `utils/` correctness (P2.1/P2.2/P2.3 first — they're real bugs; then P2.4/P2.5).
4. **P3** gitignore, CITATION, output-stripping decision, typos, CLAUDE.md.

## Verification checklist (after you implement)

- `git grep -n 'wftissueimmaifi'` → no hits in published files.
- `git ls-files | grep -c '\.ipynb_checkpoints/'` → `0`.
- `git grep -n '/home/workspace'` → only intended placeholders (ideally none outside docs).
- In a fresh checkout, edit only the documented env var / `BASE_DIR`, then run a setup cell
  (`from config.paths import ...; import plotting_utils as pu; import de_utils`) — imports
  succeed with no global side effects and no `colorcet`/unbound-name errors.
- `python -c "import ast,sys; [ast.parse(open(f).read()) for f in ['utils/plotting_utils.py','utils/de_utils.py','utils/distance_utils.py']]"` parses clean;
  optionally run a linter (`ruff`/`pyflakes`) to confirm unused imports are gone.
- README env-create commands run as written on a clean machine.
