---
name: notebook-cleanup
description: Clean a manuscript-analysis Jupyter notebook (Python or R) for publication readiness, following this repo's codereview.md. Use when preparing any notebook under initial_processing/, downstream_analysis/ (incl. 05_cellchat/), or geo/ for the public GitHub release — read the notebook, CHECK CELL IDs, propose a per-notebook edit plan, then implement on approval (NotebookEdit if IDs exist, else programmatic rebuild). This is a living checklist; update it as new patterns are found.
---

# Notebook cleanup (publication readiness)

Guided cleanup of one notebook at a time, grounded in [`codereview.md`](../../../codereview.md)
at the repo root. The standing tiers there (P0 leaks → P3 polish) are the source of truth;
this skill is the per-notebook execution procedure plus repo-specific rules.

## Workflow

1. **Read the whole notebook** first (Read tool on the `.ipynb`). Do not edit before reading —
   `NotebookEdit` requires it and you need full context (later cells reveal import usage,
   hardcoded paths, etc.).
2. **Propose a per-notebook edit plan** to the user: a short, numbered list of concrete edits
   keyed to the checklist below, naming the cell(s) and the change. Flag (do not silently make)
   any change that alters results or behavior.
3. **Wait for approval**, then implement (see "Editing mechanics" — check cell IDs first!).
   The user may ask for a subset, or ask you to plan-first and approve before *any* edit.
4. **Report** what changed and give the `nbdiff-web` command (see Verification).

## Editing mechanics — CHECK CELL IDs FIRST

Before any edit, check whether cells have real `id` fields:
`python3 -c "import json; print([c.get('id') for c in json.load(open(PATH))['cells']])"`

- **Real UUID-ish IDs** → edit normally with `NotebookEdit` (targets by `id`). Most committed
  notebooks (`initial_processing/`, top-level `downstream_analysis/`, `geo/`) have these.
- **IDs are `None`** (the Read tool then shows synthetic `cell-0`, `cell-1`, … labels) →
  `NotebookEdit`'s id targeting **misfires and scrambles/clobbers cells** (this happened on the
  first cellchat notebook and silently deleted a load cell). Do NOT use `NotebookEdit`. Instead
  **rebuild the notebook with a small `python3` script**: build the cells list in order
  (new + modified + preserved), copy `outputs`/`execution_count` from the originals to honor
  "keep outputs," assign unique `id`s (e.g. `f'cc01-{i:02d}'`), `json.dump(..., indent=1)`, then
  verify order + preserved outputs. The R cellchat notebooks had no IDs.
- Inserting two markdown cells "at the beginning" with `NotebookEdit` prepends, so insert the
  *second* cell (pinned-env) first, then the title, to end up with title → pinned-env.

## Standing directions

- **HARD RULE — never edit code that changes function or outputs without explicit approval.**
  You may implement *non-functional* edits directly (see "safe" below). Anything that could
  change what the code computes, produces, saves, or how it runs MUST be proposed and approved
  first — propose it, wait for an explicit "yes," then implement. When unsure which side an edit
  falls on, treat it as functional and ask.
  - **R notebooks:** treat removing `library(...)` calls as functional (ask first) — unlike
    Python imports, R package loads can have load-time side effects (e.g. S3 method registration).
    Also watch for functions used without their package loaded (e.g. `pal_d3` needs `library(ggsci)`) —
    flag the missing `library()` as a reproducibility fix.
  - **Safe to implement directly (non-functional):** removing genuinely-unused imports; fixing
    comments/markdown/typos; whitespace/formatting (e.g. `os.path.exists (x)` → `os.path.exists(x)`);
    deleting trailing empty cells; adding the Pinned Environment block; turning a discarded
    expression into a `print(...)` of the same value.
  - **Requires approval first (functional):** changing any computation, parameter, threshold,
    seed, or algorithm; changing what/where data is read or written (incl. swapping a data source
    like the HISE→local load); replacing a redefined function with an imported one; converting a
    no-op into an `assert`; algorithm-version switches (e.g. leiden `flavor="igraph"`); anything
    that alters control flow.
- **Do NOT strip cell outputs.** Leave all executed outputs in place for now. Do not propose
  output-stripping. (The user may do this separately later.)
- **Add a "Pinned Environment" markdown block** to each notebook if one does not already exist
  (see below).
- **Never run git commands.** The user handles all commits.

## Pinned Environment block

If the notebook lacks a pinned-environment cell, insert a markdown cell **immediately after the
title cell**. The relative-link depth depends on how deep the notebook lives:
- one level deep (`initial_processing/`, top-level `downstream_analysis/`, `geo/`): `../conda_envs/<env>.yml`
- two levels deep (`downstream_analysis/05_cellchat/`): `../../conda_envs/<env>.yml`

If a pinned-env cell exists but sits *before* the title cell, offer to reorder it after the title
(cosmetic, non-functional):

```
**Pinned Environment:** [`conda_envs/<env>.yml`](../conda_envs/<env>.yml)
```

Pick `<env>` from the mapping below (verify against the notebook's own `session_info` output /
`os.path.basename(sys.prefix)` if present; ask the user if unclear):

| Notebook | Env file |
|---|---|
| `initial_processing/00_scRNAseq_reference` | `space2_20250604.yml` |
| `initial_processing/01_label_prediction_scanvi` | `cellcharter2_20250604.yml` (GPU) |
| `initial_processing/02_label_refine` | `space2_20250604.yml` |
| `initial_processing/03_cd4_cite_integration` | `r_seurat_20250604.yml` |
| `initial_processing/04_cellcharter` | `cellcharter2_20250604.yml` (GPU) |
| `initial_processing/05_merge_cell_labels` … `08_plot_labels_zones` | `space2_20250604.yml` |
| `downstream_analysis/00`–`04`, `save_final_adata`, `geo/` | `space2_20250604.yml` |
| `downstream_analysis/05_cellchat/00_*` (Python prep) | `space2_20250604.yml` |
| `downstream_analysis/05_cellchat/01_*`, `02_*` (R CellChat) | `cellchat.yml` |

## Per-notebook cleanup checklist

Apply each item that's present; skip what doesn't apply.

1. **Secrets / leaks (P0).** Remove any leaked credentials or internal cloud identifiers
   (e.g. `gs://wftissueimmaifi-...` bucket IDs, tokens). Replace with a placeholder or delete
   the line. These can hide in markdown "terminal command" cells too.
2. **Hardcoded `/home/workspace/...` paths (P1.1).** Replace literal absolute paths in **code
   cells** with values from `config.paths` (`DATA_DIR`, `BASE_OUTDIR`, `REF_DIR`, `INPUTS_DIR`,
   `FUNCTIONS_DIR`) joined via `os.path.join`. Watch for inline `pd.read_csv('/home/workspace/...')`
   calls — these are common and easy to miss.
   - **Bootstrap depth:** the `sys.path.append` to reach the repo root depends on notebook depth —
     `parents[0]` one level deep, `parents[1]` two levels deep (`05_cellchat/`). R uses
     `source(file.path(dirname(getwd()), "config", "paths.R"))` one level / `dirname(dirname(getwd()))`
     two levels; after sourcing, `BASE_OUTDIR`/`INPUTS_DIR`/etc. are available in the R global env.
   - A notebook that has **no config bootstrap at all** (e.g. fresh `geo/`/`cellchat/` notebooks)
     needs one added (depth-appropriate) — this is the "add config path loading" task.
2b. **Internal-platform / broken data fetches.** Replace data loads that depend on Allen-internal
   tooling — e.g. HISE `hp.cache_files(...)`, a `download_files.extend(...)` referencing an
   undefined list, or hard-coded internal file UUIDs — with a local read of the artifact produced
   by an earlier notebook (e.g. `sc.read_h5ad(os.path.join(BASE_OUTDIR, ...))`). These blocks
   are usually broken for external users (undefined names) and leak internal identifiers. Flag
   for approval since it changes the data source; confirm the local path with the user.
3. **Setup-markdown accuracy.** "Local file info"–style cells should tell users to set the
   relevant variable **in `config/paths.py`** (e.g. `DATA_DIR`, `REF_DIR`), not a local notebook
   variable. Fix grammar/typos while there. **R notebooks** source `config/paths.R` instead —
   point users there (`config/paths.R`), not the `.py`.
4. **Unused imports (P2.4).** Trim the import cell to what the notebook actually uses (confirm
   each by scanning ALL cells — and note that matches inside docstrings/comments are NOT usage).
   Exception: keep the *conventional* scanpy-stack imports even when unused — `numpy`, `pandas`,
   `anndata`, `squidpy` — a reader expects these; the user's standing choice ("trim, keep
   conventional") is to keep these four and drop everything else unused (incl. `matplotlib`/`seaborn`
   when there's no plotting). When in doubt about a conventional one, confirm rather than strip
   silently. Also drop the inline global
   `warnings.simplefilter(action='ignore', category=Warning)` if present (it silences warnings
   process-wide) and a stray bare `sys.version_info` statement.
5. **Duplicated utility functions (reuse).** If the notebook redefines functions that already
   live in `utils/` (`de_utils.py`, `plotting_utils.py`, `distance_utils.py`) — e.g.
   `DE_test`, `DE_volcano`, `filter_adata_expressed_in_n_cells` — flag them and propose importing
   from the module instead (`from de_utils import ...` — the repo convention, matching
   `distance_utils.py`; or `import plotting_utils as pu`). **Diff the notebook's definition against
   the module's before replacing** (a quick `difflib` script). If identical → swap freely. If it has
   *diverged* (e.g. `01_DE`'s `run_zone_DE_analysis` used `dotplot` while `de_utils` used
   `matrixplot`), flag and offer: (a) source as-is (changes the figure/result), (b) update the module
   to match the notebook so the output is preserved — safe when that module function is used only by
   this notebook — then source, or (c) keep it inline. Get approval; this is functional.
5b. **Docstring coverage.** Every function (`def` in Python, `name <- function(...)` in R) should
   have a doc block stating purpose, params, and returns. Fill gaps: Python → numpy-style docstring;
   R notebooks → a plain comment-header block (the repo's notebooks aren't packages, so don't use
   roxygen unless asked). Documenting is non-functional/safe, but describe ONLY what the code
   actually does — if a function's intent is unclear, flag it rather than guessing. Also fix
   incomplete docstrings (e.g. an empty `param :` type annotation).
5c. **Undefined names used but never defined (flag — functional).** Scan for variables/functions
   referenced but defined nowhere in the notebook — they error on a clean run. Recurring cases seen:
   palette vars (`zone_palette_mapped`, `zone12_palette_mapped`), output dirs (`dist_out_dir` used in
   `08` but only `structures_dir` defined), and functions called but undefined (`netVisual_bubble_targetsOnly`).
   Flag with options: define it, rename to the intended existing name, or (for a whole cell that can't
   run) comment the cell out with a NOTE. Don't guess silently.
5d. **Customized copies of library functions (esp. R).** Some notebooks inline-adapt library
   internals (e.g. CellChat's `netAnalysis_signalingRole_heatmap`). Review these *carefully* for
   subtle bugs (e.g. `is.nan(x)` where the empties are `NA` → use `is.na`), document them with
   markdown explaining which stock function is adapted and what was customized (and add a doc header
   if it's a real function def), and watch for `library(Pkg())` stray-`()` typos.
6. **Silent no-op checks → assertions (correctness, flag).** A line like
   `adata.obs_names.equals(other.index)` whose return value is discarded should become
   `assert <expr>, "<message>"`. Related typo to catch: `adata.obs[col] = np` (assigns the numpy
   *module*) → almost always meant `np.nan`.
7. **Minor cleanups.** Remove stray spaces (`os.path.exists (x)` → `os.path.exists(x)`), trailing
   empty cells, commented-out dead code, and narrow bare `except:` clauses where obvious.
8. **Cross-notebook path / filename consistency (flag — functional).** Inputs one notebook reads
   should match what an earlier notebook writes. Recurring drift seen: `adata_distance*.h5ad` vs
   `adata_labeled*.h5ad`; `cite_expressed_genes_by_annotation.csv` read from
   `distance_analysis/.../R/` in one notebook but `cell_labeling/.../cite_integration/` in another;
   palettes loaded from `inputs/` vs `inputs/palettes/`. Flag mismatches and ask which path is canonical.
9. **Unused boilerplate cells.** A "Zissou colormap" definition cell is copy-pasted into many
   downstream notebooks; remove it ONLY where `colormap`/`colormap_r` are unused (they ARE used in
   the distance-vis notebooks `03`/`04`). Same idea for any defined-but-unused setup cell.

## Markdown documentation (often the explicit ask)

Improving markdown is non-functional/safe, and frequently the main request. Patterns that helped:
- **Bare notebooks** (no title / few markdown cells — e.g. fresh cellchat notebooks): add a title +
  one-paragraph overview, a `## Local file info` header, and a short `## <step>` header before each
  logical code block describing what it does and why.
- **Manual / visual-inspection workflows** (`02_label_refine`, `06_label_zones`): explicitly state
  that labels/zones are assigned by *eyeballing the inline outputs*, and describe the loop
  (subcluster/cluster → inspect markers in the plots above → hand-map clusters to labels via a
  `*_mapping`/`zone_map` dict). Keep the existing biological rationale.
- **Customized library functions**: explain which stock function is adapted and what changed (see 5d).
- Keep additions concise; don't restate code line-by-line.

## Per-folder READMEs

Each analysis folder gets a brief `README.md` (see `initial_processing/README.md` for the template):
one-paragraph intro (what it produces, run-in-order, what it consumes/feeds) + a
`Notebook | Purpose | Environment` table + a closing note pointing to `conda_envs/` and
`config/paths.py`(`.R`). Sub-pipelines (e.g. `05_cellchat/`) get their own small table; note any
non-canonical `*_variations/` subfolders.

## Flag-only (do not change without explicit OK)

These alter results/behavior — list them in the plan as flagged items, don't edit:
- Algorithm-version switches (e.g. leiden `flavor="igraph"`) — changes clustering output.
- `os.makedirs(DATA_DIR)` on a missing **input** dir — silently proceeds against empty data;
  a clear error would be better, but it's a behavior change.

## Verification (after edits)

- Show the diff: `nbdiff-web <path/to/notebook.ipynb>` (working tree vs HEAD; add `--no-outputs`
  to focus on source). Requires `nbdime`.
- Confirm no remaining leaks/paths in the edited notebook (the user can run
  `git grep -n '/home/workspace' <notebook>` and `git grep -n 'wftissueimmaifi' <notebook>`).
- Note: editing without re-running leaves execution counts out of sync until the notebook is
  next run top-to-bottom — expected, since we don't execute.
