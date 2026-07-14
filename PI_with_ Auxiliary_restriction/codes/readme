# `code/` — Main Analysis Notebook

`PI_CODE.ipynb` implements the paper's numerical illustration block by
block. Each block corresponds to a specific assumption, proposition, or
corollary in Section 2 of the paper (and its appendix extensions), in
the same order they appear there. Running the notebook top to bottom
regenerates every figure and every reported number in the paper.

## Setup

Put `cps2020_v1.csv` (from `data/`) in this same folder before running,
or edit the `pd.read_csv(...)` path in Block 01 to point at
`data/cps2020_v1.csv`. Then open `PI_CODE.ipynb` in Jupyter and run all
cells (Kernel → Restart & Run All).

## Block-by-block correspondence to the paper

| Cell(s) | Content | Paper reference | Output |
|---|---|---|---|
| 0 | Imports and plotting style | — | — |
| 1 | Data setup, Assumption 1 (`Q` positive definite) | Section 2 opening | — |
| 2 | Interval construction (interval-privacy mechanism) | Appendix A | — |
| 3–4 | Proposition 1 (`basic`) | Section 2 | `fig01_unconstrained_region` |
| 5–6 | Proposition 2 (`mean_hyperplane`) | Section 2.1 | `fig02_mean_known` |
| 7–8 | Proposition 3 (`width`) | Section 2.2 | `figB1_coord_bounds` |
| 9 | Proposition (`directional_contraction_mean`) | Section 2.3 | numbers only |
| 10–11 | Corollary (`local_mean_contraction`) | Section 2.3 | `figB2_local_quadratic_check` |
| 12–13 | Proposition (`transformed_convexity`), retargeting to $\theta_f$ | Section 2.4 | `fig03_transformed_region` |
| 14–15 | Proposition (`moment_aux_convexity`), no retargeting | Section 2.5 | `fig04_moment_only_band` |
| 16–17 | Proposition (`conditional_fwl`) | Section 2.6 | `fig05_conditional_race` |
| 18–20 | Supplementary 3-D detail for `conditional_fwl` (own-covariate tangency) | Appendix B | `figB3_case_c_3d`, `figB4_case_c_projection` |
| 21–22 | Propositions (`external_sharp`, `external_contraction`) | Section 2.7 | `fig06_external_v` |
| 23 | Corollary (`between_cell`) | Section 2.7 | numbers only |
| 24–25 | Proposition (`appendix_quantile`), known median | Appendix C | `figC1_median_quantile` |
| 26 | Results recorder | — | `results_summary.csv` |

## Two things worth knowing before you read the output

1. **Blocks 10–11 report an honest negative finding, not a bug.** The
   local-quadratic corollary requires a continuous score; on this data
   `educ` takes only 12 distinct values, so the corollary's ratio does
   not converge to 1. The same cell also reruns the identical
   construction on a synthetic continuous regressor to confirm this is
   a property of the data, not an error in the formula or the code.
2. **Block 18–20's figure needs two analytically-derived directions to
   visually show tangency.** The conditional-on-race segment's
   endpoints touch the boundary of the full 3-D region exactly, but a
   generic random sample of directions will almost never land on the
   exact direction where this happens. The code explicitly solves for
   these two directions (they are exactly reproducible, because race is
   itself a component of the covariate vector) before building the
   plotted polytope.

## Requirements

```
numpy
pandas
scipy
matplotlib
```
