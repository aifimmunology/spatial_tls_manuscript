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

- **Do NOT strip cell outputs.** Leave all executed outputs in place for now. Do not propose
  output-stripping. (The user may do this separately later.)
- **Add a "Pinned Environment" markdown block** to each notebook if one does not already exist
  (see below).
- **Never run git commands.** The user handles all commits.

## Pinned Environment block

If the notebook lacks a pinned-environment cell, insert a markdown cell **immediately after the
title cell** in this exact format (notebooks live one dir deep, so the link is `../conda_envs/...`):

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
3. **Setup-markdown accuracy.** "Local file info"–style cells should tell users to set the
   relevant variable **in `config/paths.py`** (e.g. `DATA_DIR`, `REF_DIR`), not a local notebook
   variable. Fix grammar/typos while there.
4. **Unused imports (P2.4).** Trim the import cell to what the notebook actually uses (confirm
   each by scanning all cells). Also drop the inline global
   `warnings.simplefilter(action='ignore', category=Warning)` if present (it silences warnings
   process-wide) and a stray bare `sys.version_info` statement.
5. **Duplicated utility functions (reuse).** If the notebook redefines functions that already
   live in `utils/` (`de_utils.py`, `plotting_utils.py`, `distance_utils.py`) — e.g.
   `DE_test`, `DE_volcano`, `filter_adata_expressed_in_n_cells` — flag them and propose importing
   from the module instead (`import de_utils as du` / `import plotting_utils as pu`). Confirm the
   notebook's definition matches the module's before replacing; if it has diverged, flag rather
   than silently swap.
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
