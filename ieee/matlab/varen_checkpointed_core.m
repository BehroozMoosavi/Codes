function [R, A, S] = varen_checkpointed_core(timeUTC, loadMW, opts, control)
%VAREN_CHECKPOINTED_CORE Resumable execution of the tested corrected core.
% [R,A,S] = varen_checkpointed_core(timeUTC,loadMW,opts,control)
% opts.OutputDir is required and must be dedicated to this one run.
% control.CheckpointEvery defaults to 250 completed origins.
% control.StopAfter defaults to Inf and is an ABSOLUTE completed-origin cap.
% A paused call returns R=[], A=[], S.Status='paused'. A completed call
% returns the usual R/A and S.Status='complete'. StopAfter can change on
% resume; all mathematical settings (including MaxForecasts) cannot.
%
% Numerical helper functions and the forecast calculations below are copied
% unchanged from varen_corrected_core.m tested in the user's MATLAB R2025a.
% Checkpoint/resume behavior requires its own 12/24-origin equivalence check.
% The original core file is not modified by this implementation.
%
% Resume requires exact normalized input arrays, this source file's text,
% mathematical configuration, MATLAB/Statistics Toolbox version and platform.
% A checkpoint stores complete origins only and preserves initial alpha
% selection. No forecast prefix or calibration search is rerun on resume.
% Latest/previous checkpoints protect the last valid state across an
% interrupted save. Only unreadable/invalid snapshots may fall back; an
% identity/configuration mismatch always stops. Use ONE MATLAB process per
% output directory. Do not edit or exchange checkpoint files between runs.
%
% R.TotalSeconds accumulates active compute represented by saved states and
% the successful call, excluding paused time and checkpoint I/O. Aborted work
% after the last committed checkpoint cannot be timed after a process crash.
% R.RuntimeSeconds retains each original per-model timing unchanged. S also
% reports current-call checkpoint/read/output overhead. Final R/A are cached
% before publishing files, so completed calls can rebuild outputs idempotently.
% No MATLAB execution of this checkpoint extension is claimed here.
% Revision 1.0, 2026-10-27.
%
% Original mathematical protocol documentation follows.
%
% STATUS: prepared and source-reviewed; NOT executed in MATLAB in the audit
% environment. The 12/24-origin resume comparison must pass before a full run.
%
% REQUIRED INPUT CONTRACT
%   timeUTC: timezone-aware datetime column, one true physical hour per row.
%   loadMW: positive finite observed system load, in MW, same length.
%   opts.IsObserved: logical column, same length, true for every input row.
%   The caller must verify source timestamps and observation provenance.
%   This function DOES NOT repair the CAISO Excel files' DST conventions.
%   In particular, assigning UTC to naive local-clock labels is not a repair.
%   No interpolation or fabricated target value is accepted by this contract.
%
% EXAMPLE (after an independently verified UTC data preparation):
%   opts = struct('IsObserved', isObserved, ...
%                 'OutputDir', fullfile(pwd,'varen_checkpointed_results'));
%   control = struct('CheckpointEvery',250,'StopAfter',Inf);
%   [R,A,S] = varen_checkpointed_core(timeUTC,loadMW,opts,control);
%
% CORRECTED PROTOCOL (different from the submitted experiment)
%   * Predict y(t) using only data through t-1 and known calendar covariates.
%   * z(t)=100*diff(log(y)); annual lag is exactly 8760 elapsed hours.
%   * V(t)=sample variance of z(t-Wv:t-1), requiring a full finite window.
%   * Fixed min/max calibrated on full past-only windows in local year 2019.
%   * AlphaVA(t)=max(0.001, 1-clip((V(t)-Vmin)/(Vmax-Vmin),0,1)).
%   * Each fit uses training-only predictor/response scaling; each internal
%     validation fold estimates its own scales. Zero-SD columns get scale1.
%   * Initial 720 model-ready rows select the static alpha once.
%   * Five deterministic expanding training blocks validate on disjoint
%     chronological blocks, using the same split/grid for all regressions.
%     A block fit remains fixed; each validation row uses actual lags known
%     at its one-step forecast origin. VA alpha is computed at block start
%     from strictly prior z and held fixed within that validation block.
%     This is block-forward tuning, not a complete nested rolling rerun.
%   * Lambda minimizes validation squared LEVEL error in MW^2.
%   * All penalized models solve, in their training-standardized units,
%       SSE/(2*n) + lambda*alpha*||beta||_1
%                  + lambda*(1-alpha)*||beta||_2^2/2,
%     with an unpenalized intercept. Ridge uses the exact alpha=0 solver.
%   * Seven models share identical forecast origins: persistence, weekly
%     naive, Ridge, LASSO, initially tuned StaticEN, FixedEN(alpha=.5), VA-REN.
%   * FixedEN(.5) is a prespecified timing-control; it is not chosen from
%     evaluation-sample mean alpha. No nonlinear benchmark is implemented.
%
% OUTPUTS
%   R: aligned hourly predictions, errors, chosen alpha/lambda, effective
%      penalty coefficients, active-set counts, timings, summary metrics.
%   A: raw volatility/alpha diagnostics, configuration, calibration metadata.
%   S: completion/resume status, checkpoint locations and current-call timing.
%   Required opts.OutputDir receives checkpoints and the final MAT/CSV files.
%   File writing cannot be suppressed in this checkpointed implementation.
%   No DM/GW p-values or bootstrap selection intervals are generated here.
%   For window-length robustness, set EvaluationWarmupRows to the largest
%   W being compared in EVERY run, so the forecast origins remain identical.
%
% MATLAB Statistics and Machine Learning Toolbox is required for lasso.
% No random cross-validation, parallel RNG, or silent fitting catches occur.
% Revised empirical results must be generated before camera-ready submission.

    callStart = tic;
    if nargin < 3, opts = struct(); end
    if nargin < 4, control = struct(); end
    control = checkpoint_controls(control);
    if ~isfield(opts,'OutputDir') || isempty(opts.OutputDir)
        error('VAREN:CheckpointOutput','Supply a dedicated nonempty opts.OutputDir.');
    end
    if ~isfield(opts,'IsObserved')
        error('VAREN:Provenance', ...
            'Supply opts.IsObserved after verifying the physical-hour data.');
    end
    cfg = default_options(opts);

    if ~isdatetime(timeUTC) || isempty(timeUTC.TimeZone)
        error('VAREN:TimeZone','timeUTC must be a timezone-aware datetime.');
    end
    timeUTC = timeUTC(:);
    timeUTC.TimeZone = 'UTC'; % Converts instants; never localizes naive labels.
    loadMW = double(loadMW(:));
    observed = opts.IsObserved(:);
    if ~islogical(observed) || numel(observed) ~= numel(loadMW) || ~all(observed)
        error('VAREN:Observed', ...
            'IsObserved must be logical, same length, and true for every row.');
    end
    if numel(timeUTC) ~= numel(loadMW) || any(isnat(timeUTC)) || ...
            any(~isfinite(loadMW)) || any(loadMW<=0)
        error('VAREN:Data','Require equal lengths, finite positive load, and no NaT.');
    end
    if any(diff(timeUTC) ~= hours(1))
        error('VAREN:Hourly', ...
            'Require strictly ordered, unique, uninterrupted physical hourly instants.');
    end
    if exist('lasso','file') ~= 2
        error('VAREN:Toolbox','MATLAB lasso (Statistics Toolbox) is required.');
    end
    timeLocal = timeUTC;
    timeLocal.TimeZone = 'America/Los_Angeles';
    N = numel(loadMW);
    if N <= 8761 + max(cfg.W,cfg.CalibrationRows)
        error('VAREN:Length','Insufficient history for annual lags and calibration.');
    end

    paths = checkpoint_paths(cfg.OutputDir);
    identity = checkpoint_identity(timeUTC,loadMW,observed,cfg);
    readClock = tic;
    [checkpoint,fallbackUsed,latestUsable] = read_checkpoint(paths,identity);
    resumeReadSeconds = toc(readClock);
    checkpointSeconds = 0;
    finalOutputSeconds = 0;
    resumed = ~isempty(checkpoint);
    initialCompleted = 0;
    computeBase = 0;
    if resumed
        initialCompleted = checkpoint.CompletedOrigins;
        computeBase = checkpoint.CumulativeComputeSeconds;
        fprintf('Checkpoint resume: %d/%d complete origins.\n', ...
            checkpoint.CompletedOrigins,checkpoint.TotalOrigins);
        if checkpoint.HasFinalResults
            R = checkpoint.FinalR; A = checkpoint.FinalA;
            if ~checkpoint.Finalized || ~all_final_outputs_exist(paths.OutputDir)
                outputClock = tic;
                publish_outputs(R,A,paths.OutputDir);
                finalOutputSeconds = toc(outputClock);
                checkpoint.Finalized = true;
                writeClock = tic;
                latestUsable = commit_checkpoint(checkpoint,paths,latestUsable);
                checkpointSeconds = checkpointSeconds+toc(writeClock);
            end
            S = run_status('complete',checkpoint,paths,resumed,initialCompleted, ...
                fallbackUsed,callStart,checkpointSeconds,resumeReadSeconds, ...
                finalOutputSeconds,true);
            return;
        end
        if checkpoint.CompletedOrigins>=control.StopAfter && ...
                checkpoint.CompletedOrigins<checkpoint.TotalOrigins
            R = []; A = [];
            S = run_status('paused',checkpoint,paths,resumed,initialCompleted, ...
                fallbackUsed,callStart,checkpointSeconds,resumeReadSeconds,0,false);
            return;
        end
    end

    clockStart = tic;
    logY = log(loadMW);
    z = [NaN; 100*diff(logY)];
    [X, predictorNames] = make_predictors(z,logY,timeLocal);

    % Shift BEFORE taking the moving window: row t uses z(t-1), never z(t).
    zPast = [NaN; z(1:end-1)];
    V = movvar(zPast,[cfg.Wv-1,0],0,'omitnan');
    completeCount = movsum(double(isfinite(zPast)),[cfg.Wv-1,0]);
    V(completeCount ~= cfg.Wv) = NaN;
    calibrationMask = year(timeLocal)==cfg.VolCalibrationYear & isfinite(V);
    if sum(calibrationMask)<cfg.Wv
        error('VAREN:VolCal','Too few complete volatility-calibration windows.');
    end
    vMin = min(V(calibrationMask));
    vMax = max(V(calibrationMask));
    if ~(vMax>vMin)
        error('VAREN:VolRange','Calibration volatility range must be positive.');
    end
    unclipped = (V-vMin)/(vMax-vMin);
    volNorm = min(1,max(0,unclipped));
    volNorm(~isfinite(V)) = NaN;
    alphaVA = max(cfg.AlphaFloor,1-volNorm);
    alphaVA(~isfinite(V)) = NaN;

    ready = all(isfinite(X),2) & isfinite(z) & isfinite(alphaVA) & ...
        year(timeLocal)>=cfg.FirstEvaluationYear & ...
        year(timeLocal)<=cfg.LastEvaluationYear;
    rawRows = find(ready);
    if numel(rawRows)<=max(cfg.W,cfg.CalibrationRows)
        error('VAREN:ModelReady','Insufficient complete evaluation-era rows.');
    end
    if any(diff(rawRows)~=1)
        error('VAREN:Rows','Model-ready input must be a contiguous hourly sequence.');
    end
    if timeUTC(find(calibrationMask,1,'last')) >= timeUTC(rawRows(1))
        error('VAREN:CalibrationOrder','Volatility calibration must precede model-ready data.');
    end
    Xr = X(rawRows,:);
    zr = z(rawRows);
    yr = loadMW(rawRows);
    previousY = loadMW(rawRows-1);
    alphaR = alphaVA(rawRows);
    Tr = numel(rawRows);

    % No evaluation observations are used to select the fixed static alpha.
    if resumed
        staticAlpha = checkpoint.StaticAlpha;
        alphaScores = checkpoint.AlphaScores;
        alphaCalibrationSeconds = checkpoint.AlphaCalibrationSeconds;
    else
        alphaClock = tic;
        initial = (1:cfg.CalibrationRows)';
        alphaScores = NaN(numel(cfg.AlphaGrid),1);
        for j=1:numel(cfg.AlphaGrid)
            scores = forward_cv(Xr(initial,:),zr(initial),yr(initial), ...
                previousY(initial),alphaR(initial),cfg.AlphaGrid(j),false,cfg);
            alphaScores(j) = min(scores);
        end
        % A fixed tie policy: favor the larger alpha among numerical exact ties.
        candidate = find(alphaScores==min(alphaScores),1,'last');
        staticAlpha = cfg.AlphaGrid(candidate);
        alphaCalibrationSeconds = toc(alphaClock);
    end

    firstModelRow = max([cfg.W,cfg.CalibrationRows,cfg.EvaluationWarmupRows])+1;
    if firstModelRow>Tr
        error('VAREN:Warmup','No evaluation observations remain after the requested warmup.');
    end
    evalRows = (firstModelRow:Tr)';
    if isfinite(cfg.MaxForecasts)
        evalRows = evalRows(1:min(numel(evalRows),cfg.MaxForecasts));
    end
    nForecast = numel(evalRows);
    names = {'Persistence','SeasonalNaive','Ridge','LASSO', ...
             'StaticEN','FixedEN05','VAREN'};
    M = numel(names);
    predicted = NaN(nForecast,M);
    chosenAlpha = NaN(nForecast,M);
    chosenLambda = NaN(nForecast,M);
    active = NaN(nForecast,M);
    seconds = zeros(nForecast,M);
    boundary = false(nForecast,M);
    if resumed
        if checkpoint.TotalOrigins~=nForecast || ...
                ~isequal(checkpoint.EvalRawRows,rawRows(evalRows)) || ...
                checkpoint.VolatilityMin~=vMin || checkpoint.VolatilityMax~=vMax
            error('VAREN:CheckpointAlignment','Recomputed input alignment/calibration differs.');
        end
        completed = checkpoint.CompletedOrigins;
        predicted(1:completed,:) = checkpoint.Predicted;
        chosenAlpha(1:completed,:) = checkpoint.ChosenAlpha;
        chosenLambda(1:completed,:) = checkpoint.ChosenLambda;
        active(1:completed,:) = checkpoint.Active;
        seconds(1:completed,:) = checkpoint.Seconds;
        boundary(1:completed,:) = checkpoint.Boundary;
    else
        completed = 0;
        checkpoint = struct('SchemaVersion',1,'Identity',identity, ...
            'CompletedOrigins',0,'TotalOrigins',nForecast, ...
            'EvalRawRows',rawRows(evalRows),'VolatilityMin',vMin,'VolatilityMax',vMax, ...
            'StaticAlpha',staticAlpha,'AlphaScores',alphaScores, ...
            'AlphaCalibrationSeconds',alphaCalibrationSeconds, ...
            'CumulativeComputeSeconds',0,'HasFinalResults',false, ...
            'Finalized',false,'FinalR',[],'FinalA',[]);
        computeNow = computeBase+toc(clockStart)-checkpointSeconds;
        writeClock = tic;
        checkpoint = snapshot(checkpoint,completed,predicted,chosenAlpha, ...
            chosenLambda,active,seconds,boundary,computeNow);
        latestUsable = commit_checkpoint(checkpoint,paths,latestUsable);
        checkpointSeconds = checkpointSeconds+toc(writeClock);
    end
    lastSavedCompleted = checkpoint.CompletedOrigins;
    stopAt = min(nForecast,control.StopAfter);
    try
    for i=(completed+1):stopAt
        t = evalRows(i);
        r = rawRows(t);
        fitRows = (t-cfg.W):(t-1);
        Xtr = Xr(fitRows,:);
        ztr = zr(fitRows);
        ytr = yr(fitRows);
        lagtr = previousY(fitRows);
        atr = alphaR(fitRows);
        fixedAlphas = [NaN NaN 0 1 staticAlpha cfg.FixedAlpha alphaR(t)];

        tick = tic;
        predicted(i,1) = loadMW(r-1);
        seconds(i,1) = toc(tick);
        tick = tic;
        predicted(i,2) = loadMW(r-168);
        seconds(i,2) = toc(tick);
        for m=3:M
            tick = tic;
            scores = forward_cv(Xtr,ztr,ytr,lagtr,atr, ...
                fixedAlphas(m),m==M,cfg);
            % Lambda is ascending; favor stronger shrinkage on exact ties.
            best = find(scores==min(scores),1,'last');
            lambda = cfg.LambdaGrid(best);
            [zhat,beta] = predict_path(Xtr,ztr,Xr(t,:),fixedAlphas(m),lambda,cfg);
            forecast = previousY(t)*exp(zhat/100);
            if ~isscalar(forecast) || ~isfinite(forecast) || forecast<=0
                error('VAREN:Forecast','Invalid forecast at raw row %d, model %s.',r,names{m});
            end
            predicted(i,m) = forecast;
            chosenAlpha(i,m) = fixedAlphas(m);
            chosenLambda(i,m) = lambda;
            active(i,m) = sum(beta~=0);
            boundary(i,m) = best==1 || best==numel(cfg.LambdaGrid);
            seconds(i,m) = toc(tick);
        end
        if cfg.Verbose && (i==1 || mod(i,100)==0 || i==nForecast)
            fprintf('Corrected VA-REN: %d/%d origins; UTC %s\n', ...
                i,nForecast,char(timeUTC(r)));
        end
        % Commit only after every model at this origin has succeeded.
        completed = i;
        if mod(completed,control.CheckpointEvery)==0 || completed==stopAt
            computeNow = computeBase+toc(clockStart)-checkpointSeconds;
            writeClock = tic;
            checkpoint = snapshot(checkpoint,completed,predicted,chosenAlpha, ...
                chosenLambda,active,seconds,boundary,computeNow);
            latestUsable = commit_checkpoint(checkpoint,paths,latestUsable);
            checkpointSeconds = checkpointSeconds+toc(writeClock);
            lastSavedCompleted = completed;
            fprintf('Checkpoint committed: %d/%d complete origins.\n',completed,nForecast);
        end
    end
    catch ME
        % Preserve successful origins preceding a solver error; never the
        % partially filled current row. A failed checkpoint write remains fatal.
        if completed>lastSavedCompleted
            try
                computeNow = computeBase+toc(clockStart)-checkpointSeconds;
                checkpoint = snapshot(checkpoint,completed,predicted,chosenAlpha, ...
                    chosenLambda,active,seconds,boundary,computeNow);
                commit_checkpoint(checkpoint,paths,latestUsable);
            catch checkpointError
                fprintf('Could not preserve extra completed origins: %s\n',checkpointError.message);
            end
        end
        rethrow(ME);
    end
    if completed<nForecast
        R = []; A = [];
        S = run_status('paused',checkpoint,paths,resumed,initialCompleted, ...
            fallbackUsed,callStart,checkpointSeconds,resumeReadSeconds,0,false);
        fprintf('Paused after %d origins; call again with a larger StopAfter to resume.\n',completed);
        return;
    end

    target = yr(evalRows);
    errorMW = target-predicted;
    usedRaw = rawRows(evalRows);
    R = struct();
    R.Protocol = 'Causal block-forward validation with training-only rolling scaling';
    R.TimeUTC = timeUTC(usedRaw);
    R.TimeLocal = timeLocal(usedRaw);
    R.RawIndex = usedRaw;
    R.TargetMW = target;
    R.ModelNames = names;
    R.PredictionMW = predicted;
    R.ErrorMW = errorMW;
    R.Alpha = chosenAlpha;
    R.Lambda = chosenLambda;
    R.EffectiveL1 = chosenLambda.*chosenAlpha;
    % Coefficient multiplying ||beta||_2^2 in the objective (includes 1/2).
    R.EffectiveL2 = chosenLambda.*(1-chosenAlpha)/2;
    R.ActiveSetSize = active;
    R.RuntimeSeconds = seconds;
    R.LambdaAtGridBoundary = boundary;
    R.StaticAlpha = staticAlpha;
    R.StaticAlphaCalibrationSeconds = alphaCalibrationSeconds;
    R.AlphaCalibrationTable = table(cfg.AlphaGrid(:),alphaScores, ...
        'VariableNames',{'Alpha','MinCVLevelMSE'});
    R.Metrics = accuracy_table(target,predicted,names,seconds,boundary);
    localHour = hour(R.TimeLocal);
    regimeLabels = {'All','Midday10to14','Evening16to21'};
    regimes = {true(nForecast,1),localHour>=10 & localHour<=14, ...
                                   localHour>=16 & localHour<=21};
    R.RegimeMetrics = table();
    for j=1:numel(regimes)
        ix = regimes{j};
        if ~any(ix), continue; end
        q = accuracy_table(target(ix),predicted(ix,:),names,seconds(ix,:),boundary(ix,:));
        q.Regime = repmat(regimeLabels(j),M,1);
        R.RegimeMetrics = [R.RegimeMetrics; q]; %#ok<AGROW>
    end
    R.TotalSeconds = computeBase+toc(clockStart)-checkpointSeconds;

    A = struct();
    A.Config = cfg;
    A.MATLABVersion = version;
    A.ToolboxVersions = ver;
    A.ComputerArchitecture = computer;
    A.InputTimeUTC = timeUTC;
    A.InputLoadMW = loadMW;
    A.InputIsObserved = observed;
    A.PredictorNames = predictorNames;
    A.VolatilityCalibrationMin = vMin;
    A.VolatilityCalibrationMax = vMax;
    A.VolatilityCalibrationCount = sum(calibrationMask);
    A.VolatilityCalibrationStartUTC = timeUTC(find(calibrationMask,1,'first'));
    A.VolatilityCalibrationEndUTC = timeUTC(find(calibrationMask,1,'last'));
    A.StaticCalibrationStartUTC = timeUTC(rawRows(1));
    A.StaticCalibrationEndUTC = timeUTC(rawRows(cfg.CalibrationRows));
    A.ForecastVolatility = V(usedRaw);
    A.ForecastNormalizedVolatility = volNorm(usedRaw);
    A.ForecastUnclippedVolatility = unclipped(usedRaw);
    A.ForecastAlpha = alphaVA(usedRaw);
    A.Diagnostics = volatility_table(A,regimes,regimeLabels);
    A.Scope = ['System load only; no solar/weather predictors or nonlinear model. ' ...
        'Input DST and observation provenance must be verified separately.'];
    if any(mean(boundary(:,3:end),1)>.05)
        warning('VAREN:LambdaBoundary', ...
            ['Some models select a lambda-grid endpoint in more than 5%% of origins. ' ...
             'Disclose it; change any grid using training/calibration evidence only.']);
    end
    % Store immutable completed results BEFORE publishing any final files.
    % A crash during publication can therefore be repaired without refitting
    % or changing the completed results' timing fields.
    checkpoint.FinalR = R; checkpoint.FinalA = A;
    checkpoint.HasFinalResults = true;
    checkpoint.CumulativeComputeSeconds = R.TotalSeconds;
    writeClock = tic;
    latestUsable = commit_checkpoint(checkpoint,paths,latestUsable);
    checkpointSeconds = checkpointSeconds+toc(writeClock);
    outputClock = tic;
    publish_outputs(R,A,paths.OutputDir);
    finalOutputSeconds = toc(outputClock);
    checkpoint.Finalized = true;
    writeClock = tic;
    latestUsable = commit_checkpoint(checkpoint,paths,latestUsable); %#ok<NASGU>
    checkpointSeconds = checkpointSeconds+toc(writeClock);
    S = run_status('complete',checkpoint,paths,resumed,initialCompleted, ...
        fallbackUsed,callStart,checkpointSeconds,resumeReadSeconds,finalOutputSeconds,false);
end

function cfg = default_options(opts)
    cfg = struct('W',720,'Wv',48,'K',5,'CalibrationRows',720,'EvaluationWarmupRows',720, ...
        'InitialTrainFraction',.5,'AlphaFloor',.001,'FixedAlpha',.5, ...
        'AlphaGrid',[.05 .1 .2 .3 .4 .5 .6 .7 .8 .9 1], ...
        'LambdaGrid',logspace(-5,2,60),'VolCalibrationYear',2019, ...
        'FirstEvaluationYear',2020,'LastEvaluationYear',2023, ...
        'MaxForecasts',Inf,'OutputDir',fullfile(pwd,'varen_corrected_results'), ...
        'Verbose',true,'RelTol',1e-6,'MaxIter',100000);
    f = fieldnames(opts);
    for j=1:numel(f)
        if strcmp(f{j},'IsObserved'), continue; end
        if ~isfield(cfg,f{j})
            error('VAREN:Option','Unknown option: %s',f{j});
        end
        cfg.(f{j}) = opts.(f{j});
    end
    validateattributes(cfg.W,{'numeric'},{'scalar','integer','>=',20});
    validateattributes(cfg.Wv,{'numeric'},{'scalar','integer','>=',2});
    validateattributes(cfg.K,{'numeric'},{'scalar','integer','>=',2});
    validateattributes(cfg.CalibrationRows,{'numeric'},{'scalar','integer','>=',20});
    validateattributes(cfg.EvaluationWarmupRows,{'numeric'},{'scalar','integer','nonnegative'});
    validateattributes(cfg.InitialTrainFraction,{'numeric'},{'scalar','>',0,'<',1});
    validateattributes(cfg.AlphaFloor,{'numeric'},{'scalar','>',0,'<',1});
    validateattributes(cfg.FixedAlpha,{'numeric'},{'scalar','>',0,'<=',1});
    validateattributes(cfg.AlphaGrid,{'numeric'},{'vector','finite','>',0,'<=',1});
    validateattributes(cfg.LambdaGrid,{'numeric'},{'vector','finite','positive'});
    if ~(isscalar(cfg.MaxForecasts) && ...
            (cfg.MaxForecasts==Inf || (cfg.MaxForecasts>=1 && mod(cfg.MaxForecasts,1)==0)))
        error('VAREN:MaxForecasts','MaxForecasts must be a positive integer or Inf.');
    end
    cfg.AlphaGrid = unique(cfg.AlphaGrid(:)','sorted');
    cfg.LambdaGrid = unique(cfg.LambdaGrid(:)','sorted');
    cfg.OutputDir = char(cfg.OutputDir);
    if cfg.VolCalibrationYear>=cfg.FirstEvaluationYear || ...
            cfg.LastEvaluationYear<cfg.FirstEvaluationYear
        error('VAREN:Years','Calibration year must precede the evaluation years.');
    end
end

function [X,names] = make_predictors(z,logY,timeLocal)
    changeLags = [1:24 168 8760];
    levelLags = [1 24 168 8760];
    X = NaN(numel(z),70);
    names = cell(1,70);
    j=0;
    for L=changeLags
        j=j+1; X(:,j)=lag_vector(z,L);
        names{j}=sprintf('LogDeltaLag_%d',L);
    end
    for L=levelLags
        j=j+1; X(:,j)=lag_vector(logY,L);
        names{j}=sprintf('LogLoadLag_%d',L);
    end
    for h=1:23
        j=j+1; X(:,j)=double(hour(timeLocal)==h);
        names{j}=sprintf('HourDummy_%02d',h);
    end
    for d=2:7
        j=j+1; X(:,j)=double(weekday(timeLocal)==d);
        names{j}=sprintf('WeekdayDummy_%d',d);
    end
    for m=2:12
        j=j+1; X(:,j)=double(month(timeLocal)==m);
        names{j}=sprintf('MonthDummy_%02d',m);
    end
    assert(j==70);
end

function v = lag_vector(x,L)
    v = [NaN(L,1); x(1:end-L)];
end

function scores = forward_cv(X,z,y,previousY,alphaPast,fixedAlpha,isVA,cfg)
    n = size(X,1);
    firstTrain = max(10,floor(cfg.InitialTrainFraction*n));
    edges = round(linspace(firstTrain+1,n+1,cfg.K+1));
    if any(diff(edges)<1)
        error('VAREN:CVSize','Too few observations for nonempty chronological folds.');
    end
    squaredSum = zeros(1,numel(cfg.LambdaGrid));
    count=0;
    for k=1:cfg.K
        validation = edges(k):(edges(k+1)-1);
        training = 1:(edges(k)-1);
        a = fixedAlpha;
        if isVA, a=alphaPast(validation(1)); end
        [zhat,~] = predict_path(X(training,:),z(training),X(validation,:), ...
            a,cfg.LambdaGrid,cfg);
        predicted = previousY(validation).*exp(zhat/100);
        if any(~isfinite(predicted(:)))
            error('VAREN:CVOverflow','Non-finite internal validation prediction.');
        end
        squaredSum = squaredSum + sum((y(validation)-predicted).^2,1);
        count = count + numel(validation);
    end
    scores = squaredSum/count;
    if any(~isfinite(scores))
        error('VAREN:CVLoss','Non-finite validation loss.');
    end
end

function [zhat,B] = predict_path(X,z,Xtest,alpha,lambdas,cfg)
    n = size(X,1);
    muX = mean(X,1);
    sdX = std(X,0,1);
    sdX(sdX==0)=1;
    muZ = mean(z);
    sdZ = std(z,0);
    if sdZ==0, sdZ=1; end
    Xs = (X-muX)./sdX;
    Xt = (Xtest-muX)./sdX;
    zs = (z-muZ)/sdZ;
    lambdas = lambdas(:)';
    if alpha==0
        % Explicit centering includes a genuinely unpenalized intercept.
        centerX = mean(Xs,1);
        centerZ = mean(zs);
        C = Xs-centerX;
        w = zs-centerZ;
        gram = C'*C/n;
        rhs = C'*w/n;
        B = NaN(size(X,2),numel(lambdas));
        for j=1:numel(lambdas)
            B(:,j) = (gram+lambdas(j)*eye(size(gram)))\rhs;
        end
        intercept = centerZ-centerX*B;
    else
        lastwarn('');
        [B,FI] = lasso(Xs,zs,'Alpha',alpha,'Lambda',lambdas, ...
            'Standardize',false,'Intercept',true,'CV','resubstitution', ...
            'RelTol',cfg.RelTol,'MaxIter',cfg.MaxIter);
        [fitWarning,fitWarningID] = lastwarn;
        if ~isempty(fitWarning)
            error('VAREN:FitWarning', ...
                'lasso returned warning %s: %s. Resolve before using forecasts.', ...
                fitWarningID,fitWarning);
        end
        % lasso may reorder its grid; explicitly align every column.
        [distance,order] = min(abs(FI.Lambda(:)-lambdas),[],1);
        if any(distance>1e-12*max(1,max(lambdas))) || numel(unique(order))~=numel(lambdas)
            error('VAREN:LambdaPath','Returned lasso path does not match requested grid.');
        end
        B = B(:,order);
        intercept = FI.Intercept(order);
    end
    zhat = muZ+sdZ*(Xt*B+intercept);
    if any(~isfinite(B(:))) || any(~isfinite(zhat(:)))
        error('VAREN:Fit','Non-finite coefficients or predictions.');
    end
end

function T = accuracy_table(y,pred,names,seconds,boundary)
    e = y-pred;
    rmse = sqrt(mean(e.^2,1));
    mae = mean(abs(e),1);
    mape = 100*mean(abs(e./y),1);
    reference = find(strcmp(names,'StaticEN'),1);
    improvement = 100*(rmse(reference)-rmse)/rmse(reference);
    T = table(names(:),repmat(numel(y),numel(names),1),rmse(:),mae(:), ...
        mape(:),improvement(:),sum(seconds,1)',mean(boundary,1)', ...
        'VariableNames',{'Model','N','RMSE_MW','MAE_MW','MAPE_pct', ...
            'RMSE_Improvement_pct','ForecastLoopSeconds','LambdaBoundaryFraction'});
end

function T = volatility_table(A,regimes,labels)
    T=table();
    for j=1:numel(regimes)
        ix=regimes{j};
        if ~any(ix), continue; end
        a=A.ForecastAlpha(ix);
        v=A.ForecastNormalizedVolatility(ix);
        u=A.ForecastUnclippedVolatility(ix);
        q=prctile(a,[0 5 25 50 75 95 100]);
        row=table(labels(j),sum(ix),mean(v),std(v),mean(a),std(a), ...
            q(1),q(2),q(3),q(4),q(5),q(6),q(7),mean(u<=0),mean(u>=1), ...
            mean(a==A.Config.AlphaFloor), ...
            'VariableNames',{'Regime','N','MeanVol','SDVol','MeanAlpha','SDAlpha', ...
                'MinAlpha','P05Alpha','P25Alpha','MedianAlpha','P75Alpha', ...
                'P95Alpha','MaxAlpha','LowerClipFraction','UpperClipFraction', ...
                'AlphaFloorFraction'});
        T=[T;row]; %#ok<AGROW>
    end
end

function write_outputs(R,A,outDir)
    if ~exist(outDir,'dir'), mkdir(outDir); end
    matPath=fullfile(outDir,'varen_corrected_results.mat');
    if exist(matPath,'file')
        error('VAREN:ExistingOutput', ...
            'Output already exists. Use a new OutputDir to preserve the prior run.');
    end
    save(matPath,'R','A','-v7.3');
    prediction=table(R.TimeUTC,R.TimeLocal,R.RawIndex,R.TargetMW, ...
        'VariableNames',{'TimestampUTC','TimestampPT','RawIndex','ObservedLoadMW'});
    parameters=table(R.TimeUTC,'VariableNames',{'TimestampUTC'});
    for m=1:numel(R.ModelNames)
        name=R.ModelNames{m};
        prediction.([name '_PredictionMW'])=R.PredictionMW(:,m);
        prediction.([name '_ErrorMW'])=R.ErrorMW(:,m);
        parameters.([name '_Alpha'])=R.Alpha(:,m);
        parameters.([name '_Lambda'])=R.Lambda(:,m);
        parameters.([name '_EffectiveL1'])=R.EffectiveL1(:,m);
        parameters.([name '_EffectiveL2'])=R.EffectiveL2(:,m);
        parameters.([name '_ActiveSet'])=R.ActiveSetSize(:,m);
        parameters.([name '_Seconds'])=R.RuntimeSeconds(:,m);
    end
    parameters.Variance=A.ForecastVolatility;
    parameters.NormalizedVolatility=A.ForecastNormalizedVolatility;
    parameters.UnclippedVolatility=A.ForecastUnclippedVolatility;
    writetable(prediction,fullfile(outDir,'aligned_predictions.csv'));
    writetable(parameters,fullfile(outDir,'parameters_and_runtime.csv'));
    writetable(R.Metrics,fullfile(outDir,'accuracy_common_sample.csv'));
    writetable(R.RegimeMetrics,fullfile(outDir,'accuracy_by_clock_window.csv'));
    writetable(A.Diagnostics,fullfile(outDir,'volatility_diagnostics.csv'));
    writetable(R.AlphaCalibrationTable,fullfile(outDir,'static_alpha_calibration.csv'));
end


function c = checkpoint_controls(control)
    c = struct('CheckpointEvery',250,'StopAfter',Inf);
    if ~isstruct(control) || ~isscalar(control)
        error('VAREN:CheckpointControl','control must be one structure.');
    end
    names = fieldnames(control);
    for j=1:numel(names)
        if ~isfield(c,names{j})
            error('VAREN:CheckpointControl','Unknown control field: %s',names{j});
        end
        c.(names{j})=control.(names{j});
    end
    validateattributes(c.CheckpointEvery,{'numeric'},{'scalar','finite','integer','positive'});
    if ~(isnumeric(c.StopAfter) && isreal(c.StopAfter) && isscalar(c.StopAfter) && ...
            (c.StopAfter==Inf || (isfinite(c.StopAfter) && c.StopAfter>=0 && mod(c.StopAfter,1)==0)))
        error('VAREN:CheckpointStop','StopAfter must be a nonnegative integer or Inf.');
    end
end

function p = checkpoint_paths(outDir)
    if ~exist(outDir,'dir')
        [ok,msg]=mkdir(outDir);
        if ~ok, error('VAREN:CheckpointFolder','%s',msg); end
    end
    p = struct('OutputDir',outDir, ...
        'Latest',fullfile(outDir,'varen_checkpoint.mat'), ...
        'Previous',fullfile(outDir,'varen_checkpoint_previous.mat'));
    if ~isfile(p.Latest) && ~isfile(p.Previous)
        names=final_output_names();
        if any(isfile(fullfile(string(outDir),string(names))))
            error('VAREN:ExistingOutput','Final output exists without this run''s checkpoint. Use a fresh OutputDir.');
        end
    end
end

function id = checkpoint_identity(timeUTC,loadMW,observed,cfg)
    selfPath=[mfilename('fullpath') '.m'];
    id = struct('InputTimeUTC',timeUTC,'InputLoadMW',loadMW, ...
        'InputIsObserved',observed,'MathConfig',rmfield(cfg,{'OutputDir','Verbose'}), ...
        'SourceText',fileread(selfPath),'MATLABVersion',version, ...
        'StatisticsToolbox',ver('stats'),'ComputerArchitecture',computer);
end

function [c,fallback,latestUsable] = read_checkpoint(p,identity)
    c=[]; fallback=false; latestUsable=false;
    latestExists=isfile(p.Latest); previousExists=isfile(p.Previous);
    if latestExists
        try
            candidate=read_snapshot(p.Latest);
        catch ME
            fprintf('Latest checkpoint is unreadable/invalid: %s\n',ME.message);
            candidate=[];
        end
        if ~isempty(candidate)
            check_identity(candidate.Identity,identity);
            c=candidate; latestUsable=true; return;
        end
    end
    if previousExists
        try
            candidate=read_snapshot(p.Previous);
        catch ME
            error('VAREN:CheckpointUnreadable','Previous checkpoint is also unreadable/invalid: %s',ME.message);
        end
        check_identity(candidate.Identity,identity);
        c=candidate; fallback=true;
        warning('VAREN:CheckpointFallback', ...
            'Using the previous valid checkpoint (%d completed origins). Latest was missing/unreadable.',c.CompletedOrigins);
    elseif latestExists
        error('VAREN:CheckpointUnreadable','Latest checkpoint is unreadable and no previous checkpoint is available.');
    end
end

function check_identity(saved,current)
    names=fieldnames(current);
    for j=1:numel(names)
        f=names{j};
        if ~isfield(saved,f) || ~isequaln(saved.(f),current.(f))
            error('VAREN:CheckpointIdentity', ...
                'Resume refused: %s differs from the saved run. Use the original inputs/code/settings or a fresh OutputDir.',f);
        end
    end
end

function c = read_snapshot(path)
    loaded=load(path,'checkpoint');
    if ~isfield(loaded,'checkpoint') || ~isstruct(loaded.checkpoint) || ~isscalar(loaded.checkpoint)
        error('VAREN:CheckpointInvalid','Missing scalar checkpoint structure.');
    end
    c=loaded.checkpoint;
    required={'SchemaVersion','Identity','CompletedOrigins','TotalOrigins','EvalRawRows', ...
        'VolatilityMin','VolatilityMax','StaticAlpha','AlphaScores','AlphaCalibrationSeconds', ...
        'CumulativeComputeSeconds','HasFinalResults','Finalized','FinalR','FinalA', ...
        'Predicted','ChosenAlpha','ChosenLambda','Active','Seconds','Boundary'};
    if ~all(isfield(c,required)) || c.SchemaVersion~=1 || ~isstruct(c.Identity)
        error('VAREN:CheckpointInvalid','Unsupported/incomplete checkpoint schema.');
    end
    n=c.CompletedOrigins; total=c.TotalOrigins;
    if ~(isnumeric(n)&&isscalar(n)&&isfinite(n)&&n>=0&&n==fix(n) && ...
            isnumeric(total)&&isscalar(total)&&isfinite(total)&&total>=1&&total==fix(total)&&n<=total) || ...
            ~iscolumn(c.EvalRawRows) || numel(c.EvalRawRows)~=total || ...
            any(diff(c.EvalRawRows)~=1)
        error('VAREN:CheckpointInvalid','Invalid completed-origin counts or row index.');
    end
    fields={'Predicted','ChosenAlpha','ChosenLambda','Active','Seconds','Boundary'};
    for j=1:numel(fields)
        if ~isequal(size(c.(fields{j})),[n 7])
            error('VAREN:CheckpointInvalid','Invalid completed-prefix shape: %s',fields{j});
        end
    end
    a=c.ChosenAlpha(:,3:7); lam=c.ChosenLambda(:,3:7); active=c.Active(:,3:7);
    if any(~isfinite(c.Predicted(:))) || any(c.Predicted(:)<=0) || ...
            any(~isfinite(a(:))) || any(a(:)<0 | a(:)>1) || ...
            any(~isfinite(lam(:))) || any(lam(:)<=0) || ...
            any(~isfinite(active(:))) || any(active(:)<0 | active(:)>70) || ...
            any(active(:)~=fix(active(:))) || any(~isfinite(c.Seconds(:))) || ...
            any(c.Seconds(:)<0) || ~islogical(c.Boundary)
        error('VAREN:CheckpointInvalid','A completed row contains invalid predictions/parameters/timing.');
    end
    times=[c.AlphaCalibrationSeconds c.CumulativeComputeSeconds];
    if numel(times)~=2 || any(~isfinite(times)) || any(times<0) || ...
            ~isscalar(c.StaticAlpha) || ~isfinite(c.StaticAlpha) || ...
            c.StaticAlpha<=0 || c.StaticAlpha>1 || ...
            ~isvector(c.AlphaScores) || any(~isfinite(c.AlphaScores)) || ...
            ~isscalar(c.VolatilityMin) || ~isscalar(c.VolatilityMax) || ...
            ~isfinite(c.VolatilityMin) || ~isfinite(c.VolatilityMax) || ...
            c.VolatilityMax<=c.VolatilityMin || ...
            ~islogical(c.HasFinalResults) || ~isscalar(c.HasFinalResults) || ...
            ~islogical(c.Finalized) || ~isscalar(c.Finalized) || ...
            (c.Finalized && ~c.HasFinalResults)
        error('VAREN:CheckpointInvalid','Invalid calibration/completion metadata.');
    end
    if c.HasFinalResults
        if n~=total || ~isstruct(c.FinalR) || ~isstruct(c.FinalA) || ...
                ~isequal(c.FinalR.RawIndex,c.EvalRawRows) || ...
                ~isequal(c.FinalR.PredictionMW,c.Predicted) || ...
                ~isequaln(c.FinalR.Alpha,c.ChosenAlpha) || ...
                ~isequaln(c.FinalR.Lambda,c.ChosenLambda) || ...
                ~isequal(c.FinalR.RuntimeSeconds,c.Seconds) || ...
                c.FinalR.TotalSeconds~=c.CumulativeComputeSeconds
            error('VAREN:CheckpointInvalid','Cached final results do not match the completed prefix.');
        end
    end
end

function c = snapshot(c,n,predicted,alpha,lambda,active,seconds,boundary,computeSeconds)
    c.CompletedOrigins=n;
    c.Predicted=predicted(1:n,:); c.ChosenAlpha=alpha(1:n,:);
    c.ChosenLambda=lambda(1:n,:); c.Active=active(1:n,:);
    c.Seconds=seconds(1:n,:); c.Boundary=boundary(1:n,:);
    c.CumulativeComputeSeconds=computeSeconds;
end

function latestUsable = commit_checkpoint(checkpoint,p,latestUsable)
    % Write and verify a separate file before replacing the current one.
    % Keep the previous valid generation until the new generation is usable.
    tempPath=[tempname(p.OutputDir) '.mat'];
    save(tempPath,'checkpoint','-v7.3');
    read_snapshot(tempPath);
    if latestUsable && isfile(p.Latest)
        [ok,msg]=copyfile(p.Latest,p.Previous,'f');
        if ~ok, error('VAREN:CheckpointBackup','%s',msg); end
    end
    [ok,msg]=movefile(tempPath,p.Latest,'f');
    if ~ok, error('VAREN:CheckpointCommit','%s',msg); end
    latestUsable=true;
end

function names = final_output_names()
    names={'varen_corrected_results.mat','aligned_predictions.csv', ...
        'parameters_and_runtime.csv','accuracy_common_sample.csv', ...
        'accuracy_by_clock_window.csv','volatility_diagnostics.csv', ...
        'static_alpha_calibration.csv'};
end

function yes = all_final_outputs_exist(outDir)
    names=final_output_names();
    yes=all(isfile(fullfile(string(outDir),string(names))));
    if yes
        for j=1:numel(names)
            d=dir(fullfile(outDir,names{j}));
            if isempty(d) || d.bytes==0, yes=false; return; end
        end
    end
end

function publish_outputs(R,A,outDir)
    % The original output writer receives a fresh staging directory. Until
    % every file is published the checkpoint remains Finalized=false, so an
    % interrupted publication can be reconstructed from cached exact R/A.
    stage=tempname(outDir);
    [ok,msg]=mkdir(stage);
    if ~ok, error('VAREN:FinalStage','%s',msg); end
    write_outputs(R,A,stage);
    names=final_output_names();
    for j=1:numel(names)
        [ok,msg]=movefile(fullfile(stage,names{j}),fullfile(outDir,names{j}),'f');
        if ~ok, error('VAREN:FinalPublish','%s',msg); end
    end
    if ~all_final_outputs_exist(outDir)
        error('VAREN:FinalPublish','A final output file is missing or empty.');
    end
    rmdir(stage);
end

function S = run_status(status,c,p,resumed,initialCompleted,fallback,callStart, ...
        checkpointSeconds,resumeReadSeconds,finalOutputSeconds,reused)
    S=struct('Status',status,'CompletedOrigins',c.CompletedOrigins, ...
        'TotalOrigins',c.TotalOrigins,'Resumed',resumed, ...
        'StartedFromCompletedOrigins',initialCompleted, ...
        'OriginsComputedThisCall',c.CompletedOrigins-initialCompleted, ...
        'CheckpointFile',p.Latest,'PreviousCheckpointFile',p.Previous, ...
        'CheckpointFallbackUsed',fallback, ...
        'CumulativeComputeSeconds',c.CumulativeComputeSeconds, ...
        'CurrentCallSeconds',toc(callStart), ...
        'CurrentCallCheckpointSeconds',checkpointSeconds, ...
        'CurrentCallResumeReadSeconds',resumeReadSeconds, ...
        'CurrentCallFinalOutputSeconds',finalOutputSeconds, ...
        'FinalOutputsReused',reused);
end
