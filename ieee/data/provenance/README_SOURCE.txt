CAISO OASIS ACTUAL-DEMAND SNAPSHOT: STEP 3
Prepared 27 October 2026

PURPOSE AND SOURCE CHANGE

This package supplies one consistent hourly system-demand series with
provider-stamped UTC interval endpoints. It replaces the EMS spreadsheet
series for the corrected experiment. It is not an exact reproduction of the
submitted results and does not recover separate observations by splitting
the averaged cells in the 2019 EMS workbook. All model estimates, accuracy
tables, diagnostic summaries, and statistical tests must be regenerated.

DATA IDENTITY

Publisher: California Independent System Operator (CAISO).
Service: https://oasis.caiso.com/oasisapi/SingleZip
Report: SLD_FCST
Market run: ACTUAL
Area: CA ISO-TAC
Data item: SYS_FCST_ACT_MW
Value column in the source report: MW
Source interval columns: INTERVALSTARTTIME_GMT and INTERVALENDTIME_GMT.

Although the report name contains FCST, this package selects only its ACTUAL
hourly-demand records. Day-ahead and other forecast records are excluded.
The target is system demand, not a solar-subtracted net-load series.

PERIOD AND DATA PROCESSING

The intended source period is local 1 January 2019 through 31 December 2023,
using America/Los_Angeles for calendar covariates. The UTC interval-start
range is 2019-01-01T08:00:00Z through 2024-01-01T07:00:00Z, with the final
interval ending at 2024-01-01T08:00:00Z. It contains 43,824 physical hours.

The acquisition partitions this range into 66 nonoverlapping requests of
at most 28 elapsed days. It filters the exact report/market/area/data item,
retains intervals in each requested half-open interval, and sorts by the
provider's UTC interval start. Every expected hour must occur exactly once.
Any gap, duplicate retained interval, nonpositive/nonfinite load, or interval
length other than one hour causes acquisition to fail. No value is filled,
averaged, interpolated, or adjusted by the acquisition code. ACTUAL describes
the provider's reporting status; it does not certify that the provider never
estimated or subsequently revised an underlying measurement.

FILE CONTENTS AND REPRODUCIBILITY

caiso_oasis_actual_2019_2023.csv   Canonical selected series and row lineage.
source_manifest.csv              Query URLs, retrieval times, source hashes,
                                 byte counts, and per-request record counts.
source_metadata.json             Coverage, units, checks, and protocol limits.
raw/oasis_actual_NNN.zip          Original unmodified OASIS ZIP responses.
raw/oasis_actual_NNN.json         Each response's URL/retrieval/hash record.
acquire_oasis_actual.py           Acquisition code, supplied for provenance.
independent_source_validation.json Independent raw-to-canonical checks.

Every canonical value points to a SourceArchive, SourceMember, and SourceRow.
SourceRow counts the header as row 1 of the original source CSV. The raw
archives include other areas returned by OASIS; only CA ISO-TAC records enter
the canonical study series. Keep the supplied raw responses because future
requests to the same URL may return a revised historical vintage.

MATLAB STEP

Put 451_step3_oasis_source.zip and prepare_caiso_oasis_inputs.m in the same
working folder. Leave this source ZIP intact, then run:

  [step3DataFile, step3Folder] = prepare_caiso_oasis_inputs(pwd);

The MATLAB helper works offline, extracts into a new directory, validates
the delivered records and manifest, and saves timezone-preserving MAT inputs.
It does not estimate forecasts. Return step3_validation.zip for review.

FORECASTING PROTOCOL TO USE AFTER INPUT VALIDATION

Volatility calibration: local year 2019, complete past-only 48-change windows.
Evaluation era: local years 2020-2023.
Annual lag: 8,760 elapsed hours, not a same-calendar-date operation.
Initial static-model calibration: first 720 complete model-ready rows.
Default rolling training window: 720 hours.
Expected aligned forecast origins after initial calibration: 34,343.
First expected forecast interval start: 2020-01-31T09:00:00Z
                                        (01:00 Pacific on 31 January 2020).

The interval-start time is the forecast origin for that hour's demand. The
backtest assumes the preceding completed hour's realized load is available
at that origin. These published historical actuals were retrieved in October
2026 and can contain later provider revisions. This is a retrospective
evaluation, not a reconstruction of the original real-time data vintages.

REFERENCE

CAISO OASIS Interface Specification, section on SLD_FCST and hourly ACTUAL:
https://www.caiso.com/documents/oasis-interfacespecification_v5_1_1clean_fall2017release.pdf
Exact executed requests and retrieval times are in source_manifest.csv.
