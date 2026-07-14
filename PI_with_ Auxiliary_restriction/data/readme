# `data/` — CPS Extract and Construction Code

## Contents

```
CPSdata.ipynb    -- pulls the raw extract from the Census Bureau API and cleans it
cps2020_v1.csv   -- the resulting extract, used by code/PI_CODE.ipynb
README.md        -- this file
VARIABLES.md     -- full data dictionary for every column in cps2020_v1.csv
```

## Source

The data comes from the Census Bureau's CPS ASEC (Annual Social and
Economic Supplement) API. `CPSdata.ipynb` requests the following
variables directly from the API:

```
WS_VAL, WSAL_VAL, TCERNVAL, TCWSVAL, A_SEX, PRDTRACE, A_AGE, A_HGA, A_EXPLF
```

This requires a Census Bureau API key. The notebook expects it to be
available as the environment variable `Census_API_key`, loaded via a
local script (edit the `filename = "..."` path at the top of the
notebook to point at your own key file, or set the environment variable
directly before running).

## Sample restrictions applied

`CPSdata.ipynb` restricts the raw API response to:

- `A_SEX == 1` (male respondents only)
- `WSAL_VAL > 0` (positive wage and salary income)
- `A_AGE` between 20 and 50 inclusive

These three filters are applied together in one step early in the
notebook, before any of the derived variables (`educ_numeric`,
`educ_description`, `wor_status`) are constructed. `code/PI_CODE.ipynb`
applies one further restriction on top of this extract: it keeps only
respondents with `wor_status == "Employed"` and drops the top 1% of the
income distribution.

## Derived variables

`CPSdata.ipynb` constructs three variables not present in the raw API
response, all as direct recodes of raw CPS codes:

- `educ_description`, `educ_numeric` — from `A_HGA` (highest grade
  attained)
- `wor_status` — from `A_EXPLF` (experienced labor force status)

The exact code-to-label mappings are in `VARIABLES.md`.

## Regenerating the extract

Running `CPSdata.ipynb` end to end overwrites `cps2020_v1.csv`. This is
only necessary if you want to pull a fresh extract (e.g., a different
year or a different set of sample restrictions) — the version already
in this folder is the one used to produce every result in the paper.

## See also

`VARIABLES.md` in this folder documents every column in
`cps2020_v1.csv`, including which are used directly in the analysis,
which are present but unused, and the full code-to-label crosswalk for
education and race.
