function [B,A] = varen_lsboost_core(timeUTC,loadMW,opts)
%VAREN_LSBOOST_CORE Chronologically tuned nonlinear comparator for VA-REN.
%
% [B,A] = varen_lsboost_core(timeUTC,loadMW,opts)
%
% Revision 1.0. Prepared and source-reviewed; MATLAB execution must be
% checked by the Step 6 pilot. This function does not write any files.
% opts.IsObserved is required. MaxForecasts defaults to 24, not a full run.
%
% Protocol:
%   - Identical 70 lag/calendar predictors and z=100*diff(log(load)) target
%     transformation as varen_corrected_core. Forecast load(t) from data
%     through t-1. No solar or weather predictors are introduced here.
%   - Initial 720 model-ready rows choose one tree structure using five
%     expanding chronological folds: 360 training rows, then five disjoint
%     validation blocks of 72 rows. Selection minimizes pooled MW^2 loss.
%   - Structures [MaxNumSplits,MinLeafSize], in exact-tie preference order:
%     [3,20;3,5;15,20;15,5]. LearnRate=.1 is fixed in advance.
%   - In each rolling 720-hour window, the same five chronological folds
%     choose the number of trees from [25,50,100,200]. Exact ties select
%     fewer trees. The tree structure and learning rate remain fixed.
%   - A block fit stays fixed during its validation block. Each validation
%     row uses actual lags known at that one-step forecast origin. This is
%     block-forward tuning, not a complete nested rolling rerun.
%   - Fit 200 trees per CV fold once, and use predict(...,'Learners',1:m)
%     to score each prefix. The first calibration fold also compares the
%     25-tree prefix against a separately fitted 25-tree ensemble. Failure
%     stops the run; this check's extra time is reported separately.
%   - At each forecast origin, refit the chosen number of trees on all 720
%     available training rows. Convert predicted z to load using the last
%     observed load. No feature/response scaling is required by trees.
%   - No random CV, row resampling, predictor subsampling, default Bayesian
%     optimization, early-stopping on evaluation outcomes, or silent fit
%     fallbacks. RNG is reset to a fixed seed before each fit and the
%     caller's RNG state is restored on exit, including errors.
%
% This is MATLAB least-squares boosting (LSBoost), not XGBoost or LightGBM.
% It is a nonlinear benchmark under the stated input/tuning protocol, not
% a claim to exhaust contemporary forecasting methods. Source and target
% provenance must be established before calling this function.
%
% Official MATLAB R2025a API references:
% https://www.mathworks.com/help/releases/R2025a/stats/fitrensemble.html
% https://www.mathworks.com/help/releases/R2025a/stats/templatetree.html
% https://www.mathworks.com/help/releases/R2025a/stats/classreg.learning.regr.compactregressionensemble.predict.html

    if nargin<3, opts=struct(); end
    if ~isstruct(opts) || ~isscalar(opts) || ~isfield(opts,'IsObserved')
        error('VAREN:BoostProvenance', ...
            'Supply scalar opts with IsObserved after verifying physical-hour data.');
    end
    cfg=default_options(opts);
    if ~isdatetime(timeUTC) || isempty(timeUTC.TimeZone)
        error('VAREN:BoostTimeZone','timeUTC must be a timezone-aware datetime.');
    end
    timeUTC=timeUTC(:);
    timeUTC.TimeZone='UTC';
    loadMW=double(loadMW(:));
    observed=opts.IsObserved(:);
    if ~islogical(observed) || numel(observed)~=numel(loadMW) || ~all(observed)
        error('VAREN:BoostObserved', ...
            'IsObserved must be logical, same length, and true for every row.');
    end
    if numel(timeUTC)~=numel(loadMW) || any(isnat(timeUTC)) || ...
            any(~isfinite(loadMW)) || any(loadMW<=0)
        error('VAREN:BoostData','Require equal lengths, finite positive load, and no NaT.');
    end
    if any(diff(timeUTC)~=hours(1))
        error('VAREN:BoostHourly', ...
            'Require strictly ordered, unique, uninterrupted physical hourly instants.');
    end
    if numel(loadMW)<=8761+max([cfg.W,cfg.CalibrationRows,cfg.EvaluationWarmupRows])
        error('VAREN:BoostLength','Insufficient history for annual lags and calibration.');
    end
    if exist('fitrensemble','file')~=2 || exist('templateTree','file')~=2
        error('VAREN:BoostToolbox', ...
            'MATLAB fitrensemble and templateTree (Statistics Toolbox) are required.');
    end

    callerRng=rng;
    rngCleanup=onCleanup(@()rng(callerRng)); %#ok<NASGU>
    [oldWarning,oldWarningID]=lastwarn;
    warningCleanup=onCleanup(@()lastwarn(oldWarning,oldWarningID)); %#ok<NASGU>
    coreClock=tic;
    timeLocal=timeUTC;
    timeLocal.TimeZone='America/Los_Angeles';
    logY=log(loadMW);
    z=[NaN;100*diff(logY)];
    [X,predictorNames]=make_predictors(z,logY,timeLocal);
    ready=all(isfinite(X),2) & isfinite(z) & ...
        year(timeLocal)>=cfg.FirstEvaluationYear & ...
        year(timeLocal)<=cfg.LastEvaluationYear;
    rawRows=find(ready);
    if numel(rawRows)<=max([cfg.W,cfg.CalibrationRows,cfg.EvaluationWarmupRows]) || ...
            any(diff(rawRows)~=1)
        error('VAREN:BoostModelReady', ...
            'Require a sufficiently long contiguous sequence of model-ready hours.');
    end
    Xr=X(rawRows,:);
    zr=z(rawRows);
    yr=loadMW(rawRows);
    previousY=loadMW(rawRows-1);
    initial=(1:cfg.CalibrationRows)';
    calibrationFolds=chronological_folds(cfg.CalibrationRows,cfg);
    rollingFolds=chronological_folds(cfg.W,cfg);
    nStructures=size(cfg.StructureGrid,1);
    nCycles=numel(cfg.CycleGrid);
    calibrationScores=NaN(nStructures,nCycles);
    calibrationFoldMSE=NaN(nStructures,cfg.K,nCycles);
    prefixCheck=struct();
    calibrationClock=tic;
    for j=1:nStructures
        [scores,foldMSE,check]=forward_cv(Xr(initial,:),zr(initial), ...
            yr(initial),previousY(initial),cfg.StructureGrid(j,:),cfg,j==1);
        calibrationScores(j,:)=scores;
        calibrationFoldMSE(j,:,:)=reshape(foldMSE,[1,cfg.K,nCycles]);
        if j==1, prefixCheck=check; end
        if cfg.Verbose
            fprintf('LSBoost structure calibration: %d/%d; splits=%d, leaf=%d\n', ...
                j,nStructures,cfg.StructureGrid(j,1),cfg.StructureGrid(j,2));
        end
    end
    calibrationSeconds=toc(calibrationClock)-prefixCheck.Seconds;
    bestByStructure=min(calibrationScores,[],2);
    selectedStructure=find(bestByStructure==min(bestByStructure),1,'first');
    selected=cfg.StructureGrid(selectedStructure,:);
    initialCycleIndex=find(calibrationScores(selectedStructure,:)== ...
        min(calibrationScores(selectedStructure,:)),1,'first');
    calibrationTable=make_calibration_table(cfg,calibrationScores, ...
        selectedStructure,initialCycleIndex);

    firstModelRow=max([cfg.W,cfg.CalibrationRows,cfg.EvaluationWarmupRows])+1;
    evalRows=(firstModelRow:numel(rawRows))';
    if isfinite(cfg.MaxForecasts)
        evalRows=evalRows(1:min(numel(evalRows),cfg.MaxForecasts));
    end
    nForecast=numel(evalRows);
    predicted=NaN(nForecast,1);
    selectedCycles=NaN(nForecast,1);
    totalSeconds=NaN(nForecast,1);
    cvSeconds=NaN(nForecast,1);
    fitSeconds=NaN(nForecast,1);
    cvScores=NaN(nForecast,nCycles);
    cvFoldMSE=NaN(nForecast,cfg.K,nCycles);
    trainingFirstRaw=NaN(nForecast,1);
    trainingLastRaw=NaN(nForecast,1);
    if cfg.Verbose
        fprintf('LSBoost chosen structure: splits=%d, leaf=%d, LearnRate=%g\n', ...
            selected(1),selected(2),cfg.LearnRate);
    end
    for i=1:nForecast
        originClock=tic;
        t=evalRows(i);
        r=rawRows(t);
        training=(t-cfg.W):(t-1);
        trainingFirstRaw(i)=rawRows(training(1));
        trainingLastRaw(i)=rawRows(training(end));
        if trainingLastRaw(i)~=r-1
            error('VAREN:BoostCausality','Training must end immediately before the forecast.');
        end
        cvClock=tic;
        [scores,foldMSE]=forward_cv(Xr(training,:),zr(training), ...
            yr(training),previousY(training),selected,cfg,false);
        cvSeconds(i)=toc(cvClock);
        chosen=find(scores==min(scores),1,'first');
        m=cfg.CycleGrid(chosen);
        fitClock=tic;
        mdl=fit_ensemble(Xr(training,:),zr(training),selected,m,cfg);
        zhat=predict_checked(mdl,Xr(t,:),m);
        forecast=previousY(t)*exp(zhat/100);
        if ~isscalar(forecast) || ~isfinite(forecast) || forecast<=0
            error('VAREN:BoostForecast','Invalid forecast at raw row %d.',r);
        end
        fitSeconds(i)=toc(fitClock);
        predicted(i)=forecast;
        selectedCycles(i)=m;
        cvScores(i,:)=scores;
        cvFoldMSE(i,:,:)=reshape(foldMSE,[1,cfg.K,nCycles]);
        totalSeconds(i)=toc(originClock);
        if cfg.Verbose && (i==1 || mod(i,100)==0 || i==nForecast)
            fprintf('LSBoost: %d/%d origins; UTC %s; trees=%d\n', ...
                i,nForecast,char(timeUTC(r)),m);
        end
    end

    usedRaw=rawRows(evalRows);
    B=struct();
    B.Protocol='LSBoost chronological rolling benchmark v1.0';
    B.ModelName='LSBoost';
    B.TimeUTC=timeUTC(usedRaw);
    B.TimeLocal=timeLocal(usedRaw);
    B.RawIndex=usedRaw;
    B.TargetMW=loadMW(usedRaw);
    B.PredictionMW=predicted;
    B.ErrorMW=B.TargetMW-predicted;
    B.NumLearningCycles=selectedCycles;
    B.RuntimeSeconds=totalSeconds;
    B.CVSeconds=cvSeconds;
    B.FitPredictSeconds=fitSeconds;
    B.CycleCVScores=cvScores;
    B.CycleCVFoldMSE=cvFoldMSE;
    B.CycleAtLowerBoundary=selectedCycles==cfg.CycleGrid(1);
    B.CycleAtUpperBoundary=selectedCycles==cfg.CycleGrid(end);
    B.SelectedStructureIndex=selectedStructure;
    B.SelectedMaxNumSplits=selected(1);
    B.SelectedMinLeafSize=selected(2);
    B.LearnRate=cfg.LearnRate;
    B.CalibrationSeconds=calibrationSeconds;
    B.CalibrationTable=calibrationTable;
    B.CalibrationFoldMSE=calibrationFoldMSE;
    B.InitialCalibrationNumLearningCycles=cfg.CycleGrid(initialCycleIndex);

    A=struct();
    A.Config=cfg;
    A.MATLABVersion=version;
    A.ToolboxVersions=ver('stats');
    A.ComputerArchitecture=computer;
    A.InputTimeUTC=timeUTC;
    A.InputLoadMW=loadMW;
    A.InputIsObserved=observed;
    A.PredictorNames=predictorNames;
    A.ModelReadyRawIndex=rawRows;
    A.StaticCalibrationRawIndex=rawRows(initial);
    A.StaticCalibrationStartUTC=timeUTC(rawRows(initial(1)));
    A.StaticCalibrationEndUTC=timeUTC(rawRows(initial(end)));
    A.ForecastTrainingFirstRawIndex=trainingFirstRaw;
    A.ForecastTrainingLastRawIndex=trainingLastRaw;
    A.ChronologicalFolds=rollingFolds;
    A.CalibrationChronologicalFolds=calibrationFolds;
    A.PrefixCheck=prefixCheck;
    A.Scope=['System demand; one-step forecasts from completed-hour load. ' ...
        '8760 elapsed-hour annual lag; Pacific calendar dummies. ' ...
        'Training-only chronological tuning; no solar/weather inputs. ' ...
        'Historical published actuals may be revised; availability is assumed. ' ...
        'LSBoost comparator, not XGBoost or LightGBM.'];
    B.TotalSeconds=toc(coreClock);
end

function cfg=default_options(opts)
    cfg=struct('W',720,'K',5,'CalibrationRows',720, ...
        'EvaluationWarmupRows',720,'InitialTrainFraction',.5, ...
        'FirstEvaluationYear',2020,'LastEvaluationYear',2023, ...
        'StructureGrid',[3 20;3 5;15 20;15 5], ...
        'LearnRate',.1,'CycleGrid',[25 50 100 200], ...
        'RandomSeed',451,'MaxForecasts',24,'Verbose',true);
    names=fieldnames(opts);
    for j=1:numel(names)
        if strcmp(names{j},'IsObserved'), continue; end
        if ~isfield(cfg,names{j})
            error('VAREN:BoostOption','Unknown option: %s',names{j});
        end
        cfg.(names{j})=opts.(names{j});
    end
    validateattributes(cfg.W,{'numeric'},{'scalar','integer','>=',20,'finite'});
    validateattributes(cfg.K,{'numeric'},{'scalar','integer','>=',2,'finite'});
    validateattributes(cfg.CalibrationRows,{'numeric'},{'scalar','integer','>=',20,'finite'});
    validateattributes(cfg.EvaluationWarmupRows,{'numeric'}, ...
        {'scalar','integer','nonnegative','finite'});
    validateattributes(cfg.InitialTrainFraction,{'numeric'}, ...
        {'scalar','>',0,'<',1,'finite'});
    validateattributes(cfg.FirstEvaluationYear,{'numeric'},{'scalar','integer','finite'});
    validateattributes(cfg.LastEvaluationYear,{'numeric'},{'scalar','integer','finite'});
    if cfg.LastEvaluationYear<cfg.FirstEvaluationYear
        error('VAREN:BoostYears','LastEvaluationYear must not precede FirstEvaluationYear.');
    end
    validateattributes(cfg.StructureGrid,{'numeric'}, ...
        {'2d','nonempty','integer','positive','finite','ncols',2});
    if size(unique(cfg.StructureGrid,'rows'),1)~=size(cfg.StructureGrid,1)
        error('VAREN:BoostGrid','StructureGrid must not contain duplicate rows.');
    end
    validateattributes(cfg.LearnRate,{'numeric'},{'scalar','>',0,'<=',1,'finite'});
    validateattributes(cfg.CycleGrid,{'numeric'}, ...
        {'vector','nonempty','integer','positive','finite'});
    cfg.CycleGrid=unique(cfg.CycleGrid(:)','sorted');
    if numel(cfg.CycleGrid)<2
        error('VAREN:BoostGrid','CycleGrid must contain at least two distinct values.');
    end
    validateattributes(cfg.RandomSeed,{'numeric'}, ...
        {'scalar','integer','nonnegative','<=',2^32-1,'finite'});
    if ~(isnumeric(cfg.MaxForecasts) && isscalar(cfg.MaxForecasts) && ...
            (cfg.MaxForecasts==Inf || ...
            (cfg.MaxForecasts>=1 && mod(cfg.MaxForecasts,1)==0)))
        error('VAREN:BoostMaxForecasts','MaxForecasts must be a positive integer or Inf.');
    end
    if ~(isscalar(cfg.Verbose) && ...
            (islogical(cfg.Verbose) || (isnumeric(cfg.Verbose) && ismember(cfg.Verbose,[0 1]))))
        error('VAREN:BoostVerbose','Verbose must be a logical scalar or 0/1.');
    end
    cfg.Verbose=logical(cfg.Verbose);
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

function T=make_calibration_table(cfg,scores,selectedStructure,selectedCycle)
    nStructures=size(cfg.StructureGrid,1);
    nCycles=numel(cfg.CycleGrid);
    structureIndex=repelem((1:nStructures)',nCycles);
    maxSplits=cfg.StructureGrid(structureIndex,1);
    minLeaf=cfg.StructureGrid(structureIndex,2);
    learnRate=repmat(cfg.LearnRate,nStructures*nCycles,1);
    cycles=repmat(cfg.CycleGrid(:),nStructures,1);
    mse=reshape(scores',[],1);
    chosenStructure=structureIndex==selectedStructure;
    chosenCycle=chosenStructure & cycles==cfg.CycleGrid(selectedCycle);
    T=table(structureIndex,maxSplits,minLeaf,learnRate,cycles,mse, ...
        chosenStructure,chosenCycle,'VariableNames',{'StructureIndex', ...
        'MaxNumSplits','MinLeafSize','LearnRate','NumLearningCycles', ...
        'CV_MSE_MW2','SelectedStructure','SelectedInitialCycles'});
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

