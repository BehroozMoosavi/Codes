function [B,A,S] = varen_lsboost_checkpointed_core(timeUTC,loadMW,calibration,opts,control)
%VAREN_LSBOOST_CHECKPOINTED_CORE Hourly LSBoost fits with 24-hour retuning.
%
% [B,A,S] = varen_lsboost_checkpointed_core(timeUTC,loadMW,calibration,opts,control)
% calibration.C and calibration.A are the validated Step 6b calibration.
% opts.IsObserved is required. opts.OutputDir is required for checkpoints.
% opts.MaxForecasts defaults to Inf. opts.Verbose defaults to true.
% control.CheckpointEvery=250; StopAfter=Inf (absolute completed-origin cap);
% WriteCheckpoints=true. With WriteCheckpoints=false, no state is read/saved.
% A paused call returns B=[], A=[], S.Status='paused'. A complete call
% returns all aligned B/A arrays and S.Status='complete'. No CSV is written.
%
% FROZEN COMPUTATIONAL PROTOCOL, revision 1.0:
% - The Step 6b-selected structure is fixed: max splits 3, min leaf 20.
% - LearnRate=.1, CycleGrid=[25 50 100 200 400 800 1600], seed 451.
% - Tune tree count at global forecast origins 1,25,49,... using the previous
%   720 hours and the same five expanding chronological folds in MW loss.
% - Hold that selected count for 24 physical-hour origins, but REFIT the
%   chosen ensemble on the most recent 720 hours at EVERY forecast origin.
% - This explicit computation-budget policy differs from the Step 6 pilot's
%   hourly tree-count tuning and the regressions' hourly lambda tuning.
%   It is not claimed to reproduce those different tuning policies.
% - The 70 predictors and z-to-MW transformation are unchanged. Initial
%   structure calibration is reused, not refitted or charged a second time.
% - First-origin CV must reproduce the Step 6b selected-structure scores and
%   fold losses. Later hyperparameters use only data available at retuning.
%
% CHECKPOINT CONTRACT:
% Only completed origins are committed. Cache updates from a retuning event
% become committed only after that origin's final fit/prediction succeeds.
% Source, normalized input arrays, calibration evidence, math settings,
% MATLAB/Statistics versions and platform must match exactly on resume.
% A valid checkpoint with a different identity stops; only unreadable or
% invalid latest snapshots can fall back to the previous valid generation.
% One MATLAB process must own an output folder. Preserve sources while a run
% is active. Latest/previous snapshots mitigate interrupted saves; no claim
% of universal filesystem atomicity is made.
%
% Per-origin fitting/tuning times survive resume unchanged. Cumulative
% compute excludes checkpoint I/O, resume-read time and pauses. A crash may
% leave abandoned work after the last saved origin unrecorded. Completed
% checkpoints rebuild outputs without refitting. The wrapper owns outputs.
%
% Source-reviewed only; MATLAB execution is validated by the wrapper's
% uninterrupted-versus-resumed preflight. The fitting/predictor helpers are
% copied unchanged from the MATLAB-tested Step 6 core.

    callStart=tic;
    B=[]; A=[];
    if nargin<4, opts=struct(); end
    if nargin<5, control=struct(); end
    cfg=fixed_options(opts);
    ctl=checkpoint_controls(control);
    if ctl.WriteCheckpoints && isempty(cfg.OutputDir)
        error('VAREN:BoostCheckpointOutput','Supply a dedicated opts.OutputDir.');
    end
    [timeUTC,loadMW,observed]=validate_input(timeUTC,loadMW,opts.IsObserved);
    validate_calibration(calibration,timeUTC,loadMW,observed,cfg);
    if exist('fitrensemble','file')~=2 || exist('templateTree','file')~=2
        error('VAREN:BoostToolbox','MATLAB fitrensemble and templateTree are required.');
    end
    callerRng=rng;
    rngCleanup=onCleanup(@()rng(callerRng)); %#ok<NASGU>
    [oldWarning,oldWarningID]=lastwarn;
    warningCleanup=onCleanup(@()lastwarn(oldWarning,oldWarningID)); %#ok<NASGU>
    sourcePath=[mfilename('fullpath') '.m'];
    sourceText=fileread(sourcePath);
    identity=checkpoint_identity(timeUTC,loadMW,observed,calibration,cfg,sourceText);
    paths=checkpoint_paths(cfg.OutputDir,ctl.WriteCheckpoints);
    checkpointSeconds=0; resumeReadSeconds=0;
    checkpoint=[]; fallback=false; latestUsable=false;
    if ctl.WriteCheckpoints
        readClock=tic;
        [checkpoint,fallback,latestUsable]=read_checkpoint(paths,identity);
        resumeReadSeconds=toc(readClock);
    end
    resumed=~isempty(checkpoint);
    timeLocal=timeUTC;
    timeLocal.TimeZone='America/Los_Angeles';
    logY=log(loadMW);
    z=[NaN;100*diff(logY)];
    [X,predictorNames]=make_predictors(z,logY,timeLocal);
    ready=all(isfinite(X),2) & isfinite(z) & year(timeLocal)>=2020 & year(timeLocal)<=2023;
    rawRows=find(ready);
    if ~isequal(rawRows,(8762:43824)')
        error('VAREN:BoostRows','Expected the frozen 35063 model-ready physical hours.');
    end
    Xr=X(rawRows,:); zr=z(rawRows); yr=loadMW(rawRows);
    previousY=loadMW(rawRows-1);
    evalRows=(721:numel(rawRows))';
    if isfinite(cfg.MaxForecasts)
        evalRows=evalRows(1:min(numel(evalRows),cfg.MaxForecasts));
    end
    usedRaw=rawRows(evalRows);
    total=numel(evalRows);
    data=empty_arrays(total);
    if resumed
        if checkpoint.TotalOrigins~=total || ~isequal(checkpoint.EvalRawRows,usedRaw)
            error('VAREN:BoostCheckpointGeometry','Saved forecast geometry differs.');
        end
        completed=checkpoint.CompletedOrigins;
        data=restore_prefix(data,checkpoint,completed);
        cache=cache_from_checkpoint(checkpoint);
        baseCompute=checkpoint.CumulativeComputeSeconds;
        reproMSE=checkpoint.CalibrationReproductionMaxMSE;
        reproFoldMSE=checkpoint.CalibrationReproductionMaxFoldMSE;
        if cfg.Verbose
            fprintf('LSBoost checkpoint resume: %d/%d complete origins.\n',completed,total);
        end
    else
        completed=0;
        cache=struct('OriginIndex',0,'Cycles',NaN,'Scores',NaN(1,7),'FoldMSE',NaN(5,7));
        baseCompute=0; reproMSE=NaN; reproFoldMSE=NaN;
        checkpoint=struct('SchemaVersion',1,'Identity',identity, ...
            'CompletedOrigins',0,'TotalOrigins',total,'EvalRawRows',usedRaw, ...
            'CumulativeComputeSeconds',0);
    end
    initialCompleted=completed;
    lastSaved=completed;
    retainedCompute=baseCompute;
    if ~resumed
        retainedCompute=max(0,toc(callStart)-resumeReadSeconds);
    end
    if ~resumed && ctl.WriteCheckpoints
        checkpoint=snapshot(checkpoint,data,completed,cache,reproMSE,reproFoldMSE, ...
            retainedCompute);
        check_source(sourcePath,sourceText);
        saveClock=tic;
        latestUsable=commit_checkpoint(checkpoint,paths,latestUsable);
        checkpointSeconds=checkpointSeconds+toc(saveClock);
    end
    stopAt=min(total,ctl.StopAfter);
    structure=cfg.StructureGrid(calibration.C.SelectedStructureIndex,:);
    try
        for i=(completed+1):stopAt
            originClock=tic;
            t=evalRows(i); r=rawRows(t);
            training=(t-cfg.W):(t-1);
            candidateCache=cache;
            cvTime=0;
            retune=mod(i-1,cfg.RetuneEvery)==0;
            candidateReproMSE=reproMSE; candidateReproFold=reproFoldMSE;
            if retune
                cvClock=tic;
                [scores,foldMSE]=forward_cv(Xr(training,:),zr(training), ...
                    yr(training),previousY(training),structure,cfg,false);
                cvTime=toc(cvClock);
                chosen=find(scores==min(scores),1,'first');
                candidateCache=struct('OriginIndex',i,'Cycles',cfg.CycleGrid(chosen), ...
                    'Scores',scores,'FoldMSE',foldMSE);
                if i==1
                    j=calibration.C.SelectedStructureIndex;
                    referenceScores=calibration.C.Scores(j,:);
                    referenceFold=squeeze(calibration.C.FoldMSE(j,:,:));
                    candidateReproMSE=max(abs(scores-referenceScores));
                    candidateReproFold=max(abs(foldMSE-referenceFold),[],'all');
                    if candidateReproMSE>1e-7 || candidateReproFold>1e-7
                        error('VAREN:BoostCalibrationReproduction', ...
                            'First-origin CV did not reproduce Step 6b calibration.');
                    end
                end
            elseif cache.OriginIndex~=1+cfg.RetuneEvery*floor((i-1)/cfg.RetuneEvery)
                error('VAREN:BoostTuningCache','Missing or misaligned held tuning state.');
            end
            fitClock=tic;
            mdl=fit_ensemble(Xr(training,:),zr(training),structure,candidateCache.Cycles,cfg);
            zhat=predict_checked(mdl,Xr(t,:),candidateCache.Cycles);
            forecast=previousY(t)*exp(zhat/100);
            if ~isscalar(forecast) || ~isfinite(forecast) || forecast<=0
                error('VAREN:BoostForecast','Invalid forecast at raw row %d.',r);
            end
            fitTime=toc(fitClock);
            % Advance every retained field only after this full origin succeeds.
            data.PredictionMW(i)=forecast;
            data.NumLearningCycles(i)=candidateCache.Cycles;
            data.CVSeconds(i)=cvTime;
            data.FitPredictSeconds(i)=fitTime;
            data.Retuned(i)=retune;
            data.TuningOriginRawIndex(i)=usedRaw(candidateCache.OriginIndex);
            data.CycleCVScores(i,:)=candidateCache.Scores;
            data.CycleCVFoldMSE(i,:,:)=reshape(candidateCache.FoldMSE,[1 5 7]);
            data.RuntimeSeconds(i)=toc(originClock);
            cache=candidateCache;
            reproMSE=candidateReproMSE; reproFoldMSE=candidateReproFold;
            completed=i;
            retainedCompute=baseCompute+max(0,toc(callStart)-resumeReadSeconds-checkpointSeconds);
            if cfg.Verbose && (i==1 || mod(i,100)==0 || i==total)
                fprintf('LSBoost24h: %d/%d origins; UTC %s; trees=%d; retuned=%d\n', ...
                    i,total,char(timeUTC(r)),cache.Cycles,retune);
            end
            if ctl.WriteCheckpoints && mod(completed,ctl.CheckpointEvery)==0
                check_source(sourcePath,sourceText);
                checkpoint=snapshot(checkpoint,data,completed,cache,reproMSE,reproFoldMSE, ...
                    retainedCompute);
                saveClock=tic;
                latestUsable=commit_checkpoint(checkpoint,paths,latestUsable);
                checkpointSeconds=checkpointSeconds+toc(saveClock);
                lastSaved=completed;
                if cfg.Verbose
                    fprintf('LSBoost checkpoint committed: %d/%d complete origins.\n',completed,total);
                end
            end
        end
        check_source(sourcePath,sourceText);
        if ~ctl.WriteCheckpoints || completed>lastSaved
            checkpoint=snapshot(checkpoint,data,completed,cache,reproMSE,reproFoldMSE,retainedCompute);
        end
        if ctl.WriteCheckpoints && completed>lastSaved
            saveClock=tic;
            latestUsable=commit_checkpoint(checkpoint,paths,latestUsable); %#ok<NASGU>
            checkpointSeconds=checkpointSeconds+toc(saveClock);
            lastSaved=completed;
            if cfg.Verbose
                fprintf('LSBoost checkpoint committed: %d/%d complete origins.\n',completed,total);
            end
        end
    catch ME
        if ctl.WriteCheckpoints && completed>lastSaved
            try
                check_source(sourcePath,sourceText);
                checkpoint=snapshot(checkpoint,data,completed,cache,reproMSE,reproFoldMSE, ...
                    retainedCompute);
                commit_checkpoint(checkpoint,paths,latestUsable);
                fprintf('Preserved %d complete origins after fitting error.\n',completed);
            catch saveError
                fprintf('Could not save additional completed origins: %s\n',saveError.message);
            end
        end
        rethrow(ME);
    end
    if completed<total
        S=run_status('paused',checkpoint,paths,resumed,initialCompleted,fallback, ...
            callStart,checkpointSeconds,resumeReadSeconds,ctl.WriteCheckpoints);
        if cfg.Verbose
            fprintf('LSBoost paused at %d/%d; call again with a larger StopAfter to resume.\n',completed,total);
        end
        return;
    end

    B=struct('Protocol','LSBoost hourly refit with 24-elapsed-hour retuning v1.0', ...
        'ModelName','LSBoost24h','TimeUTC',timeUTC(usedRaw),'TimeLocal',timeLocal(usedRaw), ...
        'RawIndex',usedRaw,'TargetMW',loadMW(usedRaw));
    fields=array_fields();
    for j=1:numel(fields), B.(fields{j})=checkpoint.(fields{j}); end
    B.ErrorMW=B.TargetMW-B.PredictionMW;
    B.CycleAtLowerBoundary=B.NumLearningCycles==cfg.CycleGrid(1);
    B.CycleAtUpperBoundary=B.NumLearningCycles==cfg.CycleGrid(end);
    B.SelectedStructureIndex=calibration.C.SelectedStructureIndex;
    B.SelectedMaxNumSplits=structure(1); B.SelectedMinLeafSize=structure(2);
    B.LearnRate=cfg.LearnRate; B.RetuneEvery=cfg.RetuneEvery;
    B.InitialCalibrationNumLearningCycles=calibration.C.SelectedNumLearningCycles;
    B.CalibrationSeconds=0;
    B.PriorCalibrationSeconds=calibration.C.CalibrationSeconds;
    B.CalibrationTable=calibration.C.CalibrationTable;
    B.CalibrationFoldMSE=calibration.C.FoldMSE;
    B.TotalSeconds=checkpoint.CumulativeComputeSeconds;
    A=struct();
    A.Config=cfg; A.MATLABVersion=version; A.ToolboxVersions=ver('stats');
    A.ComputerArchitecture=computer;
    A.InputTimeUTC=timeUTC; A.InputLoadMW=loadMW; A.InputIsObserved=observed;
    A.PredictorNames=predictorNames; A.ModelReadyRawIndex=rawRows;
    A.StaticCalibrationRawIndex=(8762:9481)';
    A.StaticCalibrationStartUTC=timeUTC(8762); A.StaticCalibrationEndUTC=timeUTC(9481);
    A.ForecastTrainingFirstRawIndex=usedRaw-cfg.W;
    A.ForecastTrainingLastRawIndex=usedRaw-1;
    A.ChronologicalFolds=chronological_folds(cfg.W,cfg);
    A.CalibrationReproductionMaxMSE=reproMSE;
    A.CalibrationReproductionMaxFoldMSE=reproFoldMSE;
    A.Scope=['System demand; same 70 lag/calendar predictors and 720-hour histories. ' ...
        'Tree count retuned every 24 elapsed forecast origins; ensemble refitted hourly. ' ...
        'This disclosed cost-driven tuning frequency differs from hourly EN lambda tuning. ' ...
        'CycleCVScores/FoldMSE are held from TuningOriginRawIndex between retunes. ' ...
        'No solar/weather variables. Historical published actuals may be revised; ' ...
        'completed-hour load is assumed available at the next origin.'];
    S=run_status('complete',checkpoint,paths,resumed,initialCompleted,fallback, ...
        callStart,checkpointSeconds,resumeReadSeconds,ctl.WriteCheckpoints);
end

function cfg=fixed_options(opts)
    if ~isstruct(opts) || ~isscalar(opts) || ~isfield(opts,'IsObserved')
        error('VAREN:BoostProvenance','Supply scalar opts with IsObserved.');
    end
    cfg=struct('W',720,'K',5,'CalibrationRows',720,'EvaluationWarmupRows',720, ...
        'InitialTrainFraction',.5,'FirstEvaluationYear',2020,'LastEvaluationYear',2023, ...
        'StructureGrid',[3 20;3 5;15 20;15 5],'LearnRate',.1, ...
        'CycleGrid',[25 50 100 200 400 800 1600],'RandomSeed',451, ...
        'RetuneEvery',24,'MaxForecasts',Inf,'OutputDir','','Verbose',true);
    allowed={'IsObserved','MaxForecasts','OutputDir','Verbose'};
    names=fieldnames(opts);
    for j=1:numel(names)
        if ~ismember(names{j},allowed)
            error('VAREN:BoostOption','Unknown or frozen option: %s.',names{j});
        end
        if ~strcmp(names{j},'IsObserved'), cfg.(names{j})=opts.(names{j}); end
    end
    if ~(isnumeric(cfg.MaxForecasts) && isreal(cfg.MaxForecasts) && isscalar(cfg.MaxForecasts) && ...
            (cfg.MaxForecasts==Inf || (isfinite(cfg.MaxForecasts) && cfg.MaxForecasts>=1 && ...
            cfg.MaxForecasts==fix(cfg.MaxForecasts))))
        error('VAREN:BoostMaxForecasts','MaxForecasts must be a positive integer or Inf.');
    end
    if ~((ischar(cfg.OutputDir) && (isrow(cfg.OutputDir)||isempty(cfg.OutputDir))) || ...
            (isstring(cfg.OutputDir) && isscalar(cfg.OutputDir) && ~ismissing(cfg.OutputDir)))
        error('VAREN:BoostOutputDir','OutputDir must be one path or empty.');
    end
    cfg.OutputDir=char(cfg.OutputDir);
    if ~(isscalar(cfg.Verbose) && (islogical(cfg.Verbose) || ...
            (isnumeric(cfg.Verbose)&&ismember(cfg.Verbose,[0 1]))))
        error('VAREN:BoostVerbose','Verbose must be logical or 0/1.');
    end
    cfg.Verbose=logical(cfg.Verbose);
end

function ctl=checkpoint_controls(control)
    ctl=struct('CheckpointEvery',250,'StopAfter',Inf,'WriteCheckpoints',true);
    if ~isstruct(control) || ~isscalar(control)
        error('VAREN:BoostCheckpointControl','control must be a scalar structure.');
    end
    names=fieldnames(control);
    for j=1:numel(names)
        if ~isfield(ctl,names{j}), error('VAREN:BoostCheckpointControl','Unknown control: %s',names{j}); end
        ctl.(names{j})=control.(names{j});
    end
    validateattributes(ctl.CheckpointEvery,{'numeric'},{'scalar','finite','integer','positive'});
    if ~(isnumeric(ctl.StopAfter) && isreal(ctl.StopAfter) && isscalar(ctl.StopAfter) && ...
            (ctl.StopAfter==Inf || (isfinite(ctl.StopAfter)&&ctl.StopAfter>=0&&ctl.StopAfter==fix(ctl.StopAfter))))
        error('VAREN:BoostCheckpointStop','StopAfter must be a nonnegative integer or Inf.');
    end
    if ~islogical(ctl.WriteCheckpoints) || ~isscalar(ctl.WriteCheckpoints)
        error('VAREN:BoostCheckpointControl','WriteCheckpoints must be a logical scalar.');
    end
end

function [t,y,o]=validate_input(t,y,o)
    if ~isdatetime(t) || isempty(t.TimeZone) || ~isnumeric(y) || ~isreal(y)
        error('VAREN:BoostInput','Require timezone-aware datetimes and real numeric loads.');
    end
    t=t(:); t.TimeZone='UTC'; y=double(y(:)); o=o(:);
    if numel(t)~=43824 || numel(y)~=43824 || numel(o)~=43824 || ...
            ~islogical(o) || ~all(o) || any(isnat(t)) || any(~isfinite(y)) || any(y<=0) || ...
            any(diff(t)~=hours(1)) || ...
            t(1)~=datetime(2019,1,1,8,0,0,'TimeZone','UTC') || ...
            t(end)~=datetime(2024,1,1,7,0,0,'TimeZone','UTC')
        error('VAREN:BoostInput','Require the complete contiguous observed 43824-hour source.');
    end
end

function validate_calibration(cal,t,y,o,cfg)
    if ~isstruct(cal) || ~isscalar(cal) || ~all(isfield(cal,{'C','A'}))
        error('VAREN:BoostCalibration','Supply the validated Step 6b C/A calibration.');
    end
    C=cal.C; A=cal.A;
    neededC={'ForecastOrigins','Scores','FoldMSE','SelectedStructureIndex', ...
        'SelectedMaxNumSplits','SelectedMinLeafSize','SelectedNumLearningCycles', ...
        'LearnRate','CalibrationSeconds','CalibrationTable'};
    neededA={'Config','InputTimeUTC','InputLoadMW','InputIsObserved', ...
        'CalibrationRawIndex','MATLABVersion','PrefixCheck'};
    if ~all(isfield(C,neededC)) || ~all(isfield(A,neededA)) || C.ForecastOrigins~=0 || ...
            ~isequal(A.InputTimeUTC,t(1:9481)) || ~isequal(A.InputLoadMW,y(1:9481)) || ...
            ~isequal(A.InputIsObserved,o(1:9481)) || ~isequal(A.CalibrationRawIndex,(8762:9481)') || ...
            ~strcmp(A.MATLABVersion,version) || ~A.PrefixCheck.Pass
        error('VAREN:BoostCalibration','Calibration prefix, status or MATLAB version differs.');
    end
    if ~isequal(A.Config.CycleGrid,cfg.CycleGrid) || ...
            ~isequal(A.Config.StructureGrid,cfg.StructureGrid) || ...
            A.Config.LearnRate~=cfg.LearnRate || A.Config.RandomSeed~=cfg.RandomSeed || ...
            A.Config.K~=cfg.K || A.Config.CalibrationRows~=720 || ...
            A.Config.InitialTrainFraction~=.5 || ...
            ~isequal(size(C.Scores),[4 7]) || ~isequal(size(C.FoldMSE),[4 5 7]) || ...
            any(~isfinite(C.Scores),'all') || any(C.Scores<0,'all') || ...
            any(~isfinite(C.FoldMSE),'all') || any(C.FoldMSE<0,'all')
        error('VAREN:BoostCalibration','Calibration grid or loss arrays are invalid.');
    end
    [~,j]=min(min(C.Scores,[],2)); [~,k]=min(C.Scores(j,:));
    pooled=squeeze(mean(C.FoldMSE,2));
    if max(abs(pooled-C.Scores),[],'all')>1e-7 || j~=1 || k~=6 || ...
            C.SelectedStructureIndex~=j || C.SelectedMaxNumSplits~=3 || ...
            C.SelectedMinLeafSize~=20 || C.SelectedNumLearningCycles~=800 || ...
            C.LearnRate~=.1 || ~isfinite(C.CalibrationSeconds) || C.CalibrationSeconds<0
        error('VAREN:BoostCalibration','Expected reviewed Step 6b structure 3/20 and 800 trees.');
    end
end

function id=checkpoint_identity(t,y,o,cal,cfg,source)
    calCfg=cal.A.Config;
    if isfield(calCfg,'Verbose'), calCfg=rmfield(calCfg,'Verbose'); end
    calIdentity=struct('Scores',cal.C.Scores,'FoldMSE',cal.C.FoldMSE, ...
        'SelectedStructureIndex',cal.C.SelectedStructureIndex, ...
        'SelectedNumLearningCycles',cal.C.SelectedNumLearningCycles, ...
        'Config',calCfg,'InputTimeUTC',cal.A.InputTimeUTC, ...
        'InputLoadMW',cal.A.InputLoadMW,'InputIsObserved',cal.A.InputIsObserved);
    id=struct('InputTimeUTC',t,'InputLoadMW',y,'InputIsObserved',o, ...
        'Calibration',calIdentity,'MathConfig',rmfield(cfg,{'OutputDir','Verbose'}), ...
        'SourceText',source,'MATLABVersion',version,'StatisticsToolbox',ver('stats'), ...
        'ComputerArchitecture',computer);
end

function p=checkpoint_paths(folder,enabled)
    p=struct('OutputDir',folder,'Latest','','Previous','');
    if ~enabled, return; end
    if ~isfolder(folder)
        [ok,msg]=mkdir(folder);
        if ~ok, error('VAREN:BoostCheckpointFolder','%s',msg); end
    end
    p.Latest=fullfile(folder,'lsboost_checkpoint.mat');
    p.Previous=fullfile(folder,'lsboost_checkpoint_previous.mat');
end

function fields=array_fields()
    fields={'PredictionMW','NumLearningCycles','RuntimeSeconds','CVSeconds', ...
        'FitPredictSeconds','Retuned','TuningOriginRawIndex','CycleCVScores','CycleCVFoldMSE'};
end

function d=empty_arrays(n)
    d=struct('PredictionMW',NaN(n,1),'NumLearningCycles',NaN(n,1), ...
        'RuntimeSeconds',NaN(n,1),'CVSeconds',NaN(n,1), ...
        'FitPredictSeconds',NaN(n,1),'Retuned',false(n,1), ...
        'TuningOriginRawIndex',NaN(n,1),'CycleCVScores',NaN(n,7), ...
        'CycleCVFoldMSE',NaN(n,5,7));
end

function d=restore_prefix(d,c,n)
    fields=array_fields();
    for j=1:numel(fields)
        f=fields{j};
        if strcmp(f,'CycleCVFoldMSE'), d.(f)(1:n,:,:)=c.(f);
        else, d.(f)(1:n,:)=c.(f); end
    end
end

function cache=cache_from_checkpoint(c)
    cache=struct('OriginIndex',c.LastTuningOriginIndex,'Cycles',c.LastChosenCycles, ...
        'Scores',c.LastCVScores,'FoldMSE',c.LastCVFoldMSE);
end

function c=snapshot(c,d,n,cache,repro,reproFold,compute)
    fields=array_fields();
    for j=1:numel(fields)
        f=fields{j};
        if strcmp(f,'CycleCVFoldMSE'), c.(f)=d.(f)(1:n,:,:);
        else, c.(f)=d.(f)(1:n,:); end
    end
    c.CompletedOrigins=n; c.CumulativeComputeSeconds=compute;
    c.LastTuningOriginIndex=cache.OriginIndex; c.LastChosenCycles=cache.Cycles;
    c.LastCVScores=cache.Scores; c.LastCVFoldMSE=cache.FoldMSE;
    c.CalibrationReproductionMaxMSE=repro;
    c.CalibrationReproductionMaxFoldMSE=reproFold;
end

function [c,fallback,latestUsable]=read_checkpoint(p,id)
    c=[]; fallback=false; latestUsable=false;
    latestExists=isfile(p.Latest); previousExists=isfile(p.Previous);
    if latestExists
        try
            candidate=read_snapshot(p.Latest);
        catch ME
            fprintf('Latest LSBoost checkpoint unreadable/invalid: %s\n',ME.message);
            candidate=[];
        end
        if ~isempty(candidate)
            check_identity(candidate.Identity,id);
            c=candidate; latestUsable=true; return;
        end
    end
    if previousExists
        try
            candidate=read_snapshot(p.Previous);
        catch ME
            error('VAREN:BoostCheckpointUnreadable','Previous checkpoint invalid: %s',ME.message);
        end
        check_identity(candidate.Identity,id);
        c=candidate; fallback=true;
        warning('VAREN:BoostCheckpointFallback', ...
            'Using previous LSBoost checkpoint with %d complete origins.',c.CompletedOrigins);
    elseif latestExists
        error('VAREN:BoostCheckpointUnreadable','Latest invalid and no previous snapshot exists.');
    end
end

function check_identity(saved,current)
    fields=fieldnames(current);
    for j=1:numel(fields)
        f=fields{j};
        if ~isfield(saved,f) || ~isequaln(saved.(f),current.(f))
            error('VAREN:BoostCheckpointIdentity', ...
                'Resume refused: %s differs. Restore original inputs/code/settings or use a new folder.',f);
        end
    end
end

function c=read_snapshot(path)
    q=load(path,'checkpoint');
    if ~isfield(q,'checkpoint') || ~isstruct(q.checkpoint) || ~isscalar(q.checkpoint)
        error('VAREN:BoostCheckpointInvalid','Missing scalar checkpoint.');
    end
    c=q.checkpoint;
    required=[{'SchemaVersion','Identity','CompletedOrigins','TotalOrigins','EvalRawRows', ...
        'CumulativeComputeSeconds','LastTuningOriginIndex','LastChosenCycles','LastCVScores', ...
        'LastCVFoldMSE','CalibrationReproductionMaxMSE','CalibrationReproductionMaxFoldMSE'},array_fields()];
    if ~all(isfield(c,required)) || c.SchemaVersion~=1 || ~isstruct(c.Identity) || ~isscalar(c.Identity)
        error('VAREN:BoostCheckpointInvalid','Unsupported checkpoint schema.');
    end
    n=c.CompletedOrigins; total=c.TotalOrigins;
    if ~(isnumeric(n)&&isreal(n)&&isscalar(n)&&isfinite(n)&&n>=0&&n==fix(n) && ...
            isnumeric(total)&&isreal(total)&&isscalar(total)&&isfinite(total)&&total>=1&& ...
            total<=34343&&total==fix(total)&&n<=total) || ...
            ~isequal(c.EvalRawRows,(9482:9481+total)') || ...
            ~isscalar(c.CumulativeComputeSeconds) || ~isfinite(c.CumulativeComputeSeconds) || ...
            c.CumulativeComputeSeconds<0
        error('VAREN:BoostCheckpointInvalid','Invalid origin geometry or compute time.');
    end
    fields=array_fields();
    for j=1:numel(fields)
        f=fields{j}; expected=[n 1];
        if strcmp(f,'CycleCVScores'), expected=[n 7]; end
        if strcmp(f,'CycleCVFoldMSE'), expected=[n 5 7]; end
        if ~isequal(size(c.(f)),expected)
            error('VAREN:BoostCheckpointInvalid','Invalid completed-prefix shape: %s',f);
        end
    end
    if any(~isfinite(c.PredictionMW)) || any(c.PredictionMW<=0) || ...
            any(~ismember(c.NumLearningCycles,[25 50 100 200 400 800 1600])) || ...
            any(~isfinite(c.CycleCVScores),'all') || any(c.CycleCVScores<0,'all') || ...
            any(~isfinite(c.CycleCVFoldMSE),'all') || any(c.CycleCVFoldMSE<0,'all') || ...
            any(~isfinite(c.RuntimeSeconds)) || any(c.RuntimeSeconds<0) || ...
            any(~isfinite(c.CVSeconds)) || any(c.CVSeconds<0) || ...
            any(~isfinite(c.FitPredictSeconds)) || any(c.FitPredictSeconds<0) || ...
            ~islogical(c.Retuned) || ~isequal(c.Retuned,mod((1:n)'-1,24)==0) || ...
            any(c.CVSeconds(~c.Retuned)~=0)
        error('VAREN:BoostCheckpointInvalid','Invalid predictions, tuning schedule or timings.');
    end
    if n==0
        if c.LastTuningOriginIndex~=0 || ~isscalar(c.LastChosenCycles) || ~isnan(c.LastChosenCycles) || ...
                ~isequal(size(c.LastCVScores),[1 7]) || ~all(isnan(c.LastCVScores)) || ...
                ~isequal(size(c.LastCVFoldMSE),[5 7]) || ~all(isnan(c.LastCVFoldMSE),'all') || ...
                ~isscalar(c.CalibrationReproductionMaxMSE) || ~isnan(c.CalibrationReproductionMaxMSE) || ...
                ~isscalar(c.CalibrationReproductionMaxFoldMSE) || ~isnan(c.CalibrationReproductionMaxFoldMSE)
            error('VAREN:BoostCheckpointInvalid','Invalid empty tuning cache.');
        end
    else
        tuningIndex=1+24*floor(((1:n)'-1)/24);
        [~,choice]=min(c.CycleCVScores,[],2); grid=[25 50 100 200 400 800 1600]';
        pooled=reshape(mean(c.CycleCVFoldMSE,2),n,7);
        if ~isequal(c.TuningOriginRawIndex,c.EvalRawRows(tuningIndex)) || ...
                ~isequal(c.NumLearningCycles,grid(choice)) || ...
                ~isequal(c.NumLearningCycles,c.NumLearningCycles(tuningIndex)) || ...
                ~isequal(c.CycleCVScores,c.CycleCVScores(tuningIndex,:)) || ...
                ~isequal(c.CycleCVFoldMSE,c.CycleCVFoldMSE(tuningIndex,:,:)) || ...
                max(abs(pooled-c.CycleCVScores),[],'all')>1e-7 || ...
                c.LastTuningOriginIndex~=tuningIndex(end) || ...
                c.LastChosenCycles~=c.NumLearningCycles(end) || ...
                ~isequal(c.LastCVScores,c.CycleCVScores(end,:)) || ...
                ~isequal(c.LastCVFoldMSE,squeeze(c.CycleCVFoldMSE(end,:,:)))
            error('VAREN:BoostCheckpointInvalid','Held tuning cache does not match the completed prefix.');
        end
        repro=[c.CalibrationReproductionMaxMSE c.CalibrationReproductionMaxFoldMSE];
        if numel(repro)~=2 || any(~isfinite(repro)) || any(repro<0) || any(repro>1e-7)
            error('VAREN:BoostCheckpointInvalid','Initial calibration reproduction is invalid.');
        end
    end
end

function latestUsable=commit_checkpoint(checkpoint,p,latestUsable)
    tempPath=[tempname(p.OutputDir) '.mat'];
    cleanup=onCleanup(@()remove_temp(tempPath)); %#ok<NASGU>
    save(tempPath,'checkpoint','-v7.3');
    verified=read_snapshot(tempPath);
    if ~isequaln(verified,checkpoint)
        error('VAREN:BoostCheckpointVerify','Saved checkpoint did not reproduce in memory.');
    end
    if latestUsable && isfile(p.Latest)
        [ok,msg]=copyfile(p.Latest,p.Previous,'f');
        if ~ok, error('VAREN:BoostCheckpointBackup','%s',msg); end
    end
    [ok,msg]=movefile(tempPath,p.Latest,'f');
    if ~ok, error('VAREN:BoostCheckpointCommit','%s',msg); end
    latestUsable=true;
end

function remove_temp(path)
    if isfile(path), delete(path); end
end

function check_source(path,expected)
    if ~strcmp(fileread(path),expected)
        error('VAREN:BoostCheckpointSource','Source changed during execution. Restore it before resuming.');
    end
end

function S=run_status(status,c,p,resumed,initialCompleted,fallback,callStart,cpSeconds,readSeconds,enabled)
    S=struct('Status',status,'CompletedOrigins',c.CompletedOrigins,'TotalOrigins',c.TotalOrigins, ...
        'Resumed',resumed,'StartedFromCompletedOrigins',initialCompleted, ...
        'OriginsComputedThisCall',c.CompletedOrigins-initialCompleted, ...
        'CheckpointFile',p.Latest,'PreviousCheckpointFile',p.Previous, ...
        'CheckpointFallbackUsed',fallback,'CumulativeComputeSeconds',c.CumulativeComputeSeconds, ...
        'CurrentCallSeconds',toc(callStart),'CurrentCallCheckpointSeconds',cpSeconds, ...
        'CurrentCallResumeReadSeconds',readSeconds,'WriteCheckpoints',enabled, ...
        'LastTuningOriginIndex',c.LastTuningOriginIndex,'LastNumLearningCycles',c.LastChosenCycles);
end


function T=chronological_folds(n,cfg)
    firstTrain=max(10,floor(cfg.InitialTrainFraction*n));
    edges=round(linspace(firstTrain+1,n+1,cfg.K+1));
    if any(diff(edges)<1)
        error('VAREN:BoostCVSize','Too few observations for nonempty chronological folds.');
    end
    fold=(1:cfg.K)';
    trainStart=ones(cfg.K,1);
    trainEnd=edges(1:end-1)'-1;
    validationStart=edges(1:end-1)';
    validationEnd=edges(2:end)'-1;
    trainN=trainEnd;
    validationN=validationEnd-validationStart+1;
    T=table(fold,trainStart,trainEnd,validationStart,validationEnd,trainN,validationN, ...
        'VariableNames',{'Fold','TrainStart','TrainEnd','ValidationStart', ...
        'ValidationEnd','TrainN','ValidationN'});
end

function [scores,foldMSE,check]=forward_cv(X,z,y,previousY,structure,cfg,doPrefixCheck)
    folds=chronological_folds(size(X,1),cfg);
    nCycles=numel(cfg.CycleGrid);
    squaredSum=zeros(1,nCycles);
    count=0;
    foldMSE=NaN(cfg.K,nCycles);
    check=struct();
    for k=1:cfg.K
        training=folds.TrainStart(k):folds.TrainEnd(k);
        validation=folds.ValidationStart(k):folds.ValidationEnd(k);
        if max(training)>=min(validation)
            error('VAREN:BoostCVCausality','Chronological training must precede validation.');
        end
        mdl=fit_ensemble(X(training,:),z(training),structure,max(cfg.CycleGrid),cfg);
        zhat=NaN(numel(validation),nCycles);
        for j=1:nCycles
            zhat(:,j)=predict_checked(mdl,X(validation,:),cfg.CycleGrid(j));
        end
        if doPrefixCheck && k==1
            checkClock=tic;
            small=fit_ensemble(X(training,:),z(training),structure,cfg.CycleGrid(1),cfg);
            direct=predict_checked(small,X(validation,:),cfg.CycleGrid(1));
            discrepancy=max(abs(direct-zhat(:,1)));
            tolerance=1e-10*max(1,max(abs(direct)));
            check=struct('Pass',discrepancy<=tolerance, ...
                'MaxAbsDifferenceZ',discrepancy,'ToleranceZ',tolerance, ...
                'Seconds',toc(checkClock),'PrefixCycles',cfg.CycleGrid(1), ...
                'FullCycles',max(cfg.CycleGrid),'StructureIndex',1, ...
                'TrainingRows',training(:),'ValidationRows',validation(:));
            if ~check.Pass
                error('VAREN:BoostPrefix', ...
                    'Staged prediction differs from separately fitted prefix: %.17g > %.17g.', ...
                    discrepancy,tolerance);
            end
        end
        predicted=previousY(validation).*exp(zhat/100);
        if any(~isfinite(predicted(:))) || any(predicted(:)<=0)
            error('VAREN:BoostCVForecast','Non-finite or nonpositive validation prediction.');
        end
        losses=(y(validation)-predicted).^2;
        foldMSE(k,:)=mean(losses,1);
        squaredSum=squaredSum+sum(losses,1);
        count=count+numel(validation);
    end
    scores=squaredSum/count;
    if any(~isfinite(scores)) || any(~isfinite(foldMSE(:)))
        error('VAREN:BoostCVLoss','Non-finite chronological validation loss.');
    end
end

function mdl=fit_ensemble(X,z,structure,nCycles,cfg)
    if any(~isfinite(X(:))) || any(~isfinite(z)) || size(X,1)~=numel(z)
        error('VAREN:BoostFitInput','Fit inputs must be finite and have equal row counts.');
    end
    % Every fit starts from the same recorded RNG state, independent of
    % other fits or any future checkpoint/resume partition.
    rng(cfg.RandomSeed,'twister');
    lastwarn('');
    learner=templateTree('Type','regression', ...
        'MaxNumSplits',structure(1),'MinLeafSize',structure(2), ...
        'NumVariablesToSample','all','PredictorSelection','allsplits', ...
        'Reproducible',true,'Surrogate','off','MergeLeaves','off');
    mdl=fitrensemble(X,z,'Method','LSBoost','Learners',learner, ...
        'NumLearningCycles',nCycles,'LearnRate',cfg.LearnRate, ...
        'Resample','off','CrossVal','off','ResponseTransform','none', ...
        'Options',statset('UseParallel',false));
    [message,identifier]=lastwarn;
    if ~isempty(message)
        error('VAREN:BoostFitWarning', ...
            'LSBoost returned warning %s: %s. Resolve before using forecasts.', ...
            identifier,message);
    end
    if mdl.NumTrained~=nCycles
        error('VAREN:BoostTreeCount', ...
            'Requested %d trees but LSBoost trained %d. No silent truncation is allowed.', ...
            nCycles,mdl.NumTrained);
    end
end

function zhat=predict_checked(mdl,X,m)
    if any(~isfinite(X(:))) || m>mdl.NumTrained
        error('VAREN:BoostPredictInput','Invalid prediction inputs or requested prefix.');
    end
    lastwarn('');
    zhat=predict(mdl,X,'Learners',1:m,'UseParallel',false);
    [message,identifier]=lastwarn;
    if ~isempty(message)
        error('VAREN:BoostPredictWarning', ...
            'LSBoost prediction returned warning %s: %s.',identifier,message);
    end
    zhat=zhat(:);
    if numel(zhat)~=size(X,1) || any(~isfinite(zhat))
        error('VAREN:BoostPrediction','Invalid transformed-response prediction.');
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
