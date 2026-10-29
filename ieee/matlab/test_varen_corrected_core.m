function test_varen_corrected_core()
%TEST_VAREN_CORRECTED_CORE User-run MATLAB causality and alignment smoke test.
% Prepared in an environment without MATLAB; not executed during the audit.
% Uses synthetic observed data only. It is not an empirical forecasting result.
% The small tuning grid deliberately reduces test cost; do not use its scores
% as paper results. Requires varen_corrected_core.m on the MATLAB path.

    n=10000;
    timeUTC=datetime(2019,1,1,0,0,0,'TimeZone','UTC')+hours((0:n-1)');
    h=(0:n-1)';
    loadMW=25000+2500*sin(2*pi*h/24)+1000*cos(2*pi*h/168) ...
        +700*sin(2*pi*h/(24*37))+200*sin(sqrt(2)*h);
    opts=struct('IsObserved',true(n,1),'MaxForecasts',3,'K',2, ...
        'AlphaGrid',[.3 .7],'LambdaGrid',[.01 .1], ...
        'OutputDir','','Verbose',false);
    [R,A]=varen_corrected_core(timeUTC,loadMW,opts);
    assert(size(R.PredictionMW,1)==3 && size(R.PredictionMW,2)==7);
    assert(all(R.Metrics.N==3));
    assert(all(isfinite(R.PredictionMW(:))));
    assert(all(abs(R.PredictionMW(:,1)-loadMW(R.RawIndex-1))<1e-10));
    assert(all(abs(R.PredictionMW(:,2)-loadMW(R.RawIndex-168))<1e-10));

    % Critical negative-control test: change the current target and all future
    % targets. Forecasts through that target's origin must be unchanged.
    cutoff=R.RawIndex(2);
    perturbed=loadMW;
    perturbed(cutoff:end)=1.17*perturbed(cutoff:end)+300;
    [S,B]=varen_corrected_core(timeUTC,perturbed,opts);
    before=R.RawIndex<=cutoff;
    predictionDifference=abs(R.PredictionMW(before,:)-S.PredictionMW(before,:));
    assert(max(predictionDifference(:))<1e-8, ...
        'Future target perturbation changed an earlier/current-origin forecast.');
    assert(max(abs(A.ForecastAlpha(before)-B.ForecastAlpha(before)))<1e-12, ...
        'Current/future target entered alpha.');
    lambdaDifference=abs(R.Lambda(before,3:end)-S.Lambda(before,3:end));
    assert(max(lambdaDifference(:))<1e-12, ...
        'Current/future target entered penalty selection.');
    assert(R.StaticAlpha==S.StaticAlpha);

    % Independent, direct formula check for the past-only 48-hour variance.
    z=[NaN;100*diff(log(loadMW))];
    for j=1:numel(R.RawIndex)
        t=R.RawIndex(j);
        expected=var(z(t-A.Config.Wv:t-1),0);
        assert(abs(expected-A.ForecastVolatility(j))<1e-10*max(1,expected));
    end

    % Invalid clock input must fail instead of silently repairing chronology.
    bad=timeUTC;
    bad(100)=bad(99);
    rejected=false;
    try
        varen_corrected_core(bad,loadMW,opts);
    catch ME
        rejected=strcmp(ME.identifier,'VAREN:Hourly');
    end
    assert(rejected,'Repeated physical timestamps must be rejected.');
    badOpts=opts;
    badOpts.IsObserved(100)=false;
    rejected=false;
    try
        varen_corrected_core(timeUTC,loadMW,badOpts);
    catch ME
        rejected=strcmp(ME.identifier,'VAREN:Observed');
    end
    assert(rejected,'Artificial/unobserved targets must be rejected.');
    fprintf('PASS: alignment, causal variance, future perturbation, and input checks.\n');
end
