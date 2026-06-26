---
name: notebook-cleanup
description: Clean a manuscript-analysis Jupyter notebook for publication readiness, following this repo's codereview.md. Use when preparing an initial_processing/ or downstream_analysis/ notebook for the public GitHub release — read the notebook, propose a per-notebook edit plan, then implement on approval. This is a living checklist; update it as new patterns are found.
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
3. **Wait for approval**, then implement with `NotebookEdit`. The user may ask for a subset.
4. **Report** what changed and give the `nbdiff-web` command (see Verification).

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
title cell** in this exact format (notebooks live one dir deep, so the link is `../conda_envs/...`).
If a pinned-env cell exists but sits *before* the title cell, offer to reorder it after the title
(cosmetic, non-functional):

```
**Pinned Environment:** [`conda_envs/<env>.yml`](../conda_envs/<env>.yml)
```

Pick `<env>` from the README/CLAUDE notebook→environment mapping (verify against the notebook's
own `session_info` output / `os.path.basename(sys.prefix)` if present; ask the user if unclear):

| Notebook | Env file |
|---|---|
| `initial_processing/00_scRNAseq_reference` | `space2_20250604.yml` |
| `initial_processing/01_label_prediction_scanvi` | `cellcharter2_20250604.yml` (GPU) |
| `initial_processing/02_label_refine` | `space2_20250604.yml` |
| `initial_processing/03_cd4_cite_integration` | `r_seurat_20250604.yml` |
| `initial_processing/04_cellcharter` | `cellcharter2_20250604.yml` (GPU) |
| `initial_processing/05_merge_cell_labels` … `08_plot_labels_zones` | `space2_20250604.yml` |
| `downstream_analysis/` (in-scope notebooks) | `space2_20250604.yml` |

## Per-notebook cleanup checklist

Apply each item that's present; skip what doesn't apply.

1. **Secrets / leaks (P0).** Remove any leaked credentials or internal cloud identifiers
   (e.g. `gs://wftissueimmaifi-...` bucket IDs, tokens). Replace with a placeholder or delete
   the line. These can hide in markdown "terminal command" cells too.
2. **Hardcoded `/home/workspace/...` paths (P1.1).** Replace literal absolute paths in **code
   cells** with values from `config.paths` (`DATA_DIR`, `BASE_OUTDIR`, `REF_DIR`, `INPUTS_DIR`,
   `FUNCTIONS_DIR`) joined via `os.path.join`. Watch for inline `pd.read_csv('/home/workspace/...')`
   calls — these are common and easy to miss.
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
   Exception: a few imports are *conventional* in this stack (e.g. `anndata`/`squidpy` in a
   scanpy notebook); if an unused import is one a reader would expect to see, confirm with the
   user before removing rather than stripping it silently. Also drop the inline global
   `warnings.simplefilter(action='ignore', category=Warning)` if present (it silences warnings
   process-wide) and a stray bare `sys.version_info` statement.
5. **Duplicated utility functions (reuse).** If the notebook redefines functions that already
   live in `utils/` (`de_utils.py`, `plotting_utils.py`, `distance_utils.py`) — e.g.
   `DE_test`, `DE_volcano`, `filter_adata_expressed_in_n_cells` — flag them and propose importing
   from the module instead (`import de_utils as du` / `import plotting_utils as pu`). Confirm the
   notebook's definition matches the module's before replacing; if it has diverged, flag rather
   than silently swap.
5b. **Docstring coverage.** Every function (`def` in Python, `name <- function(...)` in R) should
   have a doc block stating purpose, params, and returns. Fill gaps: Python → numpy-style docstring;
   R notebooks → a plain comment-header block (the repo's notebooks aren't packages, so don't use
   roxygen unless asked). Documenting is non-functional/safe, but describe ONLY what the code
   actually does — if a function's intent is unclear, flag it rather than guessing. Also fix
   incomplete docstrings (e.g. an empty `param :` type annotation).
6. **Silent no-op checks → assertions (correctness).** A line like
   `adata.obs_names.equals(other.index)` whose return value is discarded should become
   `assert <expr>, "<message>"`.
7. **Minor cleanups.** Remove stray spaces (`os.path.exists (x)` → `os.path.exists(x)`), trailing
   empty cells, commented-out dead code, and narrow bare `except:` clauses where obvious.

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
