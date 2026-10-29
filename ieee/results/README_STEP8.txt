VA-REN REVISION: STEP 8 SAVED-FORECAST ANALYSIS
Revision 1.0

PURPOSE
This step analyzes the completed Step 7 forecasts. It fits no forecasting
models and changes no forecasts, predictors, grids, calibration choices,
rolling windows, or tuning schedules. The Step 7 folder is read only.

The MATLAB implementation has been source-reviewed. Its statistical
formulas have an independent Python arithmetic reference. MATLAB was not
available in the preparation environment, so the wrapper first executes
deterministic numerical self-tests in your MATLAB session. A successful
analysis run validates execution; it does not establish model superiority.

FILES
  run_varen_step8_analysis.m       Run the complete analysis and package it.
  varen_forecast_inference.m       Numerical inference without file writes.
  test_varen_forecast_inference.m  Deterministic numerical self-tests.
  README_STEP8.txt                 These instructions and the protocol.

INSTALL AND RUN
Extract the three .m files beside your earlier scripts in:
  /Users/behroozmoosavi/Desktop/ieee

Keep your complete Step 7 output folder. The Step 7 results ZIP contains
the input/results/source files needed if that folder must be restored.

In MATLAB:

  cd('/Users/behroozmoosavi/Desktop/ieee');
  clear run_varen_step8_analysis varen_forecast_inference
  clear test_varen_forecast_inference

  [step8Folder, step8Zip] = run_varen_step8_analysis( ...
      '/Users/behroozmoosavi/Desktop/ieee/step7_full_tp2b7a8ed2_2b71_4adc_961d_617f7ae720a7');

Alternatively call run_varen_step8_analysis with no argument and select
step7_full_results.mat in the completed Step 7 folder. You may also pass
that MAT file's full path directly.

The wrapper writes to a new step8_analysis_* folder beside itself. It does
not overwrite the Step 7 results. The inference core uses base MATLAB
arithmetic, erfc, and gammainc; the earlier Statistics Toolbox installation
can remain unchanged. No new parallel or GPU configuration is needed.

This is postprocessing rather than the many-hour forecasting computation.
Its actual runtime is printed. Wait for the ZIP path and the MATLAB prompt.
Return step8_analysis_results.zip. If it fails, return the complete error
and step8_failure.zip if generated. Preserve the output folder.

FIXED ANALYSIS CHOICES
1. Sample: all 34,343 common, contiguous physical hourly forecast origins.
   Target: reported CAISO system demand in MW, one hour ahead. Calendar
   indicators use America/Los_Angeles; time dependence uses elapsed hours.
2. Three paired comparisons, always with VA-REN as model A:
      B = calibrated StaticEN
      B = FixedEN05
      B = LSBoost24h
3. Primary reported loss: squared forecast error. Absolute error is a
   secondary loss. MAPE is descriptive.
4. Bartlett/Newey-West HAC maximum lag: 168 elapsed hours for the main
   report; 24, 48, and 720 hours are declared sensitivity checks. These
   choices are not selected by their p-values. The 720-hour check reflects
   the length of the rolling training window; it is not an assertion that
   all dependence ends at 720 hours.
5. Mean-loss tests are two-sided. Holm adjustments cover the three planned
   comparison pairs separately for each loss, lag, and test type. The main
   inferential family is squared-loss mean comparisons at lag 168. Other
   families are secondary/sensitivity analyses, not alternative ways to
   obtain a favorable primary result.
6. Effect intervals are asymptotic, paired, pointwise 95% confidence
   intervals. They are not simultaneous familywise intervals, prediction
   intervals, or an equivalence test.
7. The supplementary joint conditional-moment test uses four instruments:
   constant, the existing causal normalized volatility, midday indicator
   (10:00--14:00 Pacific inclusive), and evening indicator (16:00--21:00
   inclusive). These quantities are available at the forecast origin.
8. These are revision-stage exploratory analyses. Forecast errors and
   aggregate metrics have already been examined during this project; this
   protocol is not claimed to be an untouched preregistered experiment.

EFFECT DIRECTION AND FORMULAS
Let e_A,t and e_B,t be target minus forecast for the same hourly target.
For squared or absolute loss, d_t = loss_B,t - loss_A,t. Positive values
favor VA-REN; negative values favor the comparator.

For a centered scalar or vector series u_t with n rows, the HAC estimate is

  Gamma_j = sum(t=j+1:n, u_t * u_(t-j)') / n,
  Omega_L = Gamma_0 + sum(j=1:L,
                (1-j/(L+1)) * (Gamma_j + Gamma_j')).

The normalization is n at every lag, not n-j. No gaps are compressed into
consecutive hours. There is no automatic lag-zero rule merely because the
forecast horizon is one hour.

The mean-loss statistic is z = sqrt(n)*mean(d)/sqrt(Omega_L). The two-sided
normal p-value is erfc(abs(z)/sqrt(2)); the mean-loss interval uses
mean(d) +/- 1.959963984540054*sqrt(Omega_L/n).

Let a=mean(e_A^2), b=mean(e_B^2). The reported RMSE improvement is
sqrt(b)-sqrt(a), in MW. Its centered influence series is

  (e_B,t^2-b)/(2*sqrt(b)) - (e_A,t^2-a)/(2*sqrt(a)).

The percentage improvement is 100*(1-sqrt(a/b)). Its centered influence is

  -50*(e_A,t^2-a)/sqrt(a*b)
       + 50*sqrt(a)*(e_B,t^2-b)/(b^(3/2)).

Applying the same HAC calculation to each influence series gives its
standard error sqrt(Omega_L/n). This includes the pairing/covariance
between the forecasts. Separate unpaired standard errors are not used.
MAE improvement uses mean(abs(e_B)-abs(e_A)) and its paired HAC variance.

JOINT CONDITIONAL-MOMENT TEST
With h_t the four-vector of instruments, form g_t=h_t*d_t and calculate

  W = n * mean(g)' * inverse(HAC(g)) * mean(g).

The implementation uses a Cholesky solve, not an explicit matrix inverse.
The reference distribution is chi-square with four degrees of freedom
when the long-run covariance of the four moment series is nonsingular.
A positive column scaling
is used only for numerical conditioning; it does not change the moment
conditions or the Wald statistic. No ridge or pseudoinverse is silently
substituted for a singular covariance.

The null tested is joint orthogonality to these four instruments. It is
implied by equal conditional predictive accuracy, but a finite instrument
set does not test every possible conditioning variable. Its p-value is
an upper-tail chi-square probability. Rejection is nondirectional and
does not by itself show that VA-REN is better overall or in every regime.
Giacomini and White (2006), Theorem 1 and Comment 5, allow a HAC covariance
for this construction; their Section 3.4 discusses unconditional loss
comparisons for methods using rolling estimation.

UNDEFINED CASES AND MULTIPLICITY
The calibrated StaticEN selected alpha=1 and produces the same forecasts
as LASSO. They remain separate labeled output columns for continuity, but
represent one forecasting rule. No inferential comparison between those
identical outputs is included in the three planned pairs.

If a requested test has zero variance or a singular covariance, the core
records an explicit undefined status and NaN p-value. It does not invent
p=1 or evidence of equivalence. For Holm ordering only, an undefined test
occupies its planned family slot with p=1; reported p-values remain NaN.
Thus an undefined test cannot shrink the planned three-comparison family.

INTERPRETATION AND LIMITS
The tests and confidence intervals use large-sample approximations under
appropriate finite-moment, dependence, and covariance-estimation
conditions. Source/alignment checks do not prove those stochastic
conditions. Persistent changes in the mean loss or other nonstationarity
can affect the interpretation of pooled inference. The yearly and
clock-window tables are descriptive diagnostics, not a claim that a
single pooled model is adequate for every part of the evaluation period.
The inference is conditional on the frozen data, calibration, and method
choices; it does not quantify the uncertainty of the full research and
model-selection process. The 8760-hour predictor lag also extends raw-data
memory beyond the 720-hour regression fitting window.

Sensitivity results must be considered together. A failure to reject
equal accuracy does not establish equivalence. Statistical significance
does not establish an operationally material improvement or cost saving.
The nonlinear comparison describes this LSBoost configuration and budget,
not all boosting implementations or all modern forecasting models.

The Step 7 grid and tuning policy remain frozen. Upper-bound tree-count
selections are disclosed. This analysis does not expand the grid using
evaluation errors. LSBoost retunes count every 24 elapsed origins and
refits hourly; Step 5 retunes lambda hourly. Runtime measurements come
from separate runs, not a controlled repeated timing benchmark.

The target remains system demand, not constructed solar-subtracted net
load. No solar/weather predictors have been added. Published actuals may
contain historical revisions; completed-hour load is assumed available
at the next origin. The annual lag remains 8760 elapsed hours.

REFERENCES
Diebold, F.X., and Mariano, R.S. (1995). Comparing Predictive Accuracy.
  Journal of Business & Economic Statistics 13(3), 253--263.
  https://doi.org/10.1080/07350015.1995.10524599
Newey, W.K., and West, K.D. (1987). A Simple, Positive Semi-Definite,
  Heteroskedasticity and Autocorrelation Consistent Covariance Matrix.
  Econometrica 55(3), 703--708. https://doi.org/10.2307/1913610
Giacomini, R., and White, H. (2006). Tests of Conditional Predictive Ability.
  Econometrica 74(6), 1545--1578.
  https://doi.org/10.1111/j.1468-0262.2006.00718.x
Holm, S. (1979). A Simple Sequentially Rejective Multiple Test Procedure.
  Scandinavian Journal of Statistics 6(2), 65--70.
  https://www.jstor.org/stable/4615733
