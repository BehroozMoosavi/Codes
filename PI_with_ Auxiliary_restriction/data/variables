# Data Dictionary: `cps2020_v1.csv`

Every column in the CPS ASEC extract, built by `CPSdata.ipynb` from the
Census Bureau API and used by `code/PI_CODE.ipynb`. Raw extract: 24,905
rows, 13 columns, before the further filtering applied inside
`PI_CODE.ipynb` (employed respondents only, top income percentile
dropped), which leaves $n=22{,}397$.

## 1. Variables used directly in the analysis

| Column | Type | Range / categories | Used as | Notes |
|---|---|---|---|---|
| `WSAL_VAL` | numeric | positive, annual $ | `ystar` (income, rescaled to $1000s) | Wage and salary income — the point-valued outcome treated as latent and coarsened into an interval by the interval-privacy mechanism. |
| `educ_numeric` | numeric | 0–20 (years), 12 distinct values | `educ` | Recoded from `A_HGA` via `educ_years_dict` (below). Enters $\tilde x=(1,\mathrm{educ})$. Only 12 distinct values because several `A_HGA` categories collapse to the same year count — this is the direct cause of the discreteness issue documented for Blocks 10–11 in `code/README.md`. |
| `PRDTRACE` | numeric code | 1–26 | `race` (recoded to binary: 1 if `PRDTRACE==1`, else 2) | Detailed race code. Code 1 is "White Only" (see `race_dict` below). Collapsed to a binary own-covariate $x_1$ for the conditional-mean result (Blocks 16–20). |
| `A_AGE` | numeric | 20–50 (years) | `age` | Used as the external variable $v$ (Blocks 21–23), binned into 5 quantile-based cells. Restricted to 20–50 **by construction** in `CPSdata.ipynb` (see below), not filtered later. |
| `wor_status` | categorical (derived) | `Employed`, `Unemployed`, `Not in experienced labor force` | sample filter | `PI_CODE.ipynb` keeps only `Employed`. Derived from `A_EXPLF` via `work_dict` (below). |

## 2. Variables present in the file but not used

| Column | Type | Notes |
|---|---|---|
| `WS_VAL` | numeric | A separate wage-value field from the API; not referenced in the analysis. |
| `TCERNVAL` | 0/1 flag | Top-code/allocation flag on earnings value, from the raw API. Not used. |
| `TCWSVAL` | 0/1 flag | Top-code/allocation flag on wage-salary value, from the raw API. Not used. |
| `A_SEX` | numeric | **Constant = 1 for every row.** This is a deliberate sample restriction, not a data artifact: `CPSdata.ipynb` filters to `A_SEX == 1` (male respondents) at extraction time, alongside the income and age filters. Because there is no variation left, `A_SEX` cannot be used as a covariate (Assumption 1 requires $Q$ nonsingular), which is why it is dropped from the analysis. |
| `A_HGA` | numeric code | 31–46. Raw CPS education-attainment code; source for `educ_numeric`/`educ_description`. Not used directly. |
| `A_EXPLF` | numeric code | 0/1/2. Raw labor-force status code; source for `wor_status`. Not used directly. |
| `PRDTRACE.1` | numeric | Exact duplicate of `PRDTRACE` (the variable list passed to the API repeats `PRDTRACE`). Not used. |
| `educ_description` | text (derived) | Human-readable label matching `educ_numeric`, e.g. "Bachelor's degree (BA,AB,BS)" for 16. Kept for reference/auditing. |

## 3. Sample restrictions applied in `CPSdata.ipynb`

Applied together, at extraction time, before any derived variable is
constructed:

```python
df_filtered = df[(df["A_SEX"] == 1) & (df["WSAL_VAL"] > 0) & (df["A_AGE"].between(20, 50))].copy()
```

- `A_SEX == 1`: male respondents only.
- `WSAL_VAL > 0`: positive wage and salary income.
- `A_AGE` between 20 and 50 inclusive.

`code/PI_CODE.ipynb` applies one further restriction on top of this:
`wor_status == "Employed"`, and drops the top 1% of the income
distribution.

## 4. Code-to-label crosswalks (from `CPSdata.ipynb`)

### Education (`A_HGA` → `educ_numeric` / `educ_description`)

| `A_HGA` code | `educ_numeric` | `educ_description` |
|---|---|---|
| 0 | 0 | Children |
| 31 | 0 | Less Than 1st Grade |
| 32 | 4 | 1st,2nd,3rd,or 4th grade |
| 33 | 6 | 5th Or 6th Grade |
| 34 | 8 | 7th and 8th grade |
| 35 | 9 | 9th Grade |
| 36 | 10 | 10th Grade |
| 37 | 11 | 11th Grade |
| 38 | 12 | 12th Grade No Diploma |
| 39 | 12 | High school graduate-high school diploma |
| 40 | 14 | Some College But No Degree |
| 41 | 14 | Assc degree-occupation/vocation |
| 42 | 14 | Assc degree-academic program |
| 43 | 16 | Bachelor's degree (BA,AB,BS) |
| 44 | 18 | Master's degree (MA,MS,MENG,MED,MSW,MBA) |
| 45 | 20 | Professional school degree (MD,DDS,DVM,L...) |
| 46 | 20 | Doctorate degree (PHD,EDD) |

Note that three distinct codes (40, 41, 42) all map to `educ_numeric=14`,
and two (45, 46) both map to 20 — this is why `educ_numeric` has only 12
distinct values despite 17 raw `A_HGA` categories.

### Employment status (`A_EXPLF` → `wor_status`)

| `A_EXPLF` code | `wor_status` |
|---|---|
| 0 | Not in experienced labor force |
| 1 | Employed |
| 2 | Unemployed |

### Race (`PRDTRACE`, reference labels — full detail, not all used)

| `PRDTRACE` code | Label |
|---|---|
| 1 | White Only |
| 2 | Black Only |
| 3 | American Indian, Alaskan Native Only |
| 4 | Asian Only |
| 5 | Hawaiian/Pacific Islander Only |
| 6–26 | Multiple-race combinations (see `CPSdata.ipynb` for the full mapping) |

`code/PI_CODE.ipynb` collapses this to a binary `race` variable: 1 if
`PRDTRACE==1` (White Only), 2 otherwise.
