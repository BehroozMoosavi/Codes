function [C,A] = varen_lsboost_calibration(timeUTC,loadMW,opts)
%VAREN_LSBOOST_CALIBRATION Bounded calibration-only LSBoost tree-count audit.
%
% [C,A] = varen_lsboost_calibration(timeUTC,loadMW,opts)
%
% Revision 1.0. Source-reviewed; MATLAB execution remains the user's check.
% This function NEVER forecasts an evaluation origin and writes no files.
% It accepts exactly raw rows 1:9481 of the frozen CAISO input, ending at
% 2020-01-31T08:00:00Z. The first evaluation target at 09:00Z and all later
% observations must be removed BEFORE calling this function.
%
% Required opts.IsObserved: true logical flags for the supplied prefix.
% Only the optional opts.Verbose can be changed. All mathematical settings
% are fixed here so this is one bounded calibration audit, not an open-ended
% automatic search or a response to evaluation forecast errors.
%
% The 720 model-ready calibration rows are raw rows 8762:9481. The same
% 70 predictors, log-change response, MW loss, five chronological folds,
% four tree structures and fixed learning rate from the Step 6 pilot are
% retained. The sole search extension is the tree-count grid:
%   [25 50 100 200 400 800 1600].
% There is no automatic extension beyond 1600 or shrinking of this grid.
%
% The 1600-tree path is fitted once per structure/fold. Candidate prefixes
% are scored using the verified predict(...,'Learners',1:m) convention.
% The first 25-tree prefix is also checked against a separate 25-tree fit.
% The wrapper must reconcile all original 25/50/100/200 calibration losses
% with the saved Step 6 calibration before accepting this audit.
%
% C retains all candidate/fold losses, the calibration minimizer, and
% measured calibration timings. These are not full-period forecasts or a
% decision about the final rolling-search budget. A documents the exact
% historical prefix, predictor names, folds, API environment and check.
%
% The fitting and prediction helpers are byte-identical to the user-tested
% Step 6 core. Every independent fit resets seed 451. The caller's RNG and
% lastwarn state are restored on normal return or error.
%
% Official MATLAB R2025a API references:
% https://www.mathworks.com/help/releases/R2025a/stats/fitrensemble.html
% https://www.mathworks.com/help/releases/R2025a/stats/templatetree.html
% https://www.mathworks.com/help/releases/R2025a/stats/classreg.learning.regr.compactregressionensemble.predict.html

    if nargin<3, opts=struct(); end
    if ~isstruct(opts) || ~isscalar(opts) || ~isfield(opts,'IsObserved')
        error('VAREN:BoostCalProvenance', ...
            'Supply scalar opts with IsObserved for the truncated historical prefix.');
    end
    cfg=fixed_options(opts);
    if ~isdatetime(timeUTC) || isempty(timeUTC.TimeZone)
        error('VAREN:BoostCalTimeZone','timeUTC must be a timezone-aware datetime.');
    end
    if ~isnumeric(loadMW) || ~isreal(loadMW)
        error('VAREN:BoostCalData','Load values must be real numeric values.');
    end
    timeUTC=timeUTC(:);
    timeUTC.TimeZone='UTC';
    loadMW=double(loadMW(:));
    observed=opts.IsObserved(:);
    if numel(timeUTC)~=9481 || numel(loadMW)~=9481 || numel(observed)~=9481
        error('VAREN:BoostCalPrefix', ...
            'Supply exactly raw rows 1:9481. Evaluation targets must be excluded.');
    end
    if ~islogical(observed) || ~all(observed) || any(isnat(timeUTC)) || ...
            any(~isfinite(loadMW)) || any(loadMW<=0)
        error('VAREN:BoostCalData', ...
            'Require finite positive observed loads, logical true flags, and no NaT.');
    end
    firstExpected=datetime(2019,1,1,8,0,0,'TimeZone','UTC');
    lastExpected=datetime(2020,1,31,8,0,0,'TimeZone','UTC');
    firstForecast=datetime(2020,1,31,9,0,0,'TimeZone','UTC');
    if timeUTC(1)~=firstExpected || timeUTC(end)~=lastExpected || ...
            any(diff(timeUTC)~=hours(1)) || timeUTC(end)>=firstForecast
        error('VAREN:BoostCalPrefix', ...
            'The historical prefix must run hourly from 2019-01-01T08Z to 2020-01-31T08Z.');
    end
    if exist('fitrensemble','file')~=2 || exist('templateTree','file')~=2
        error('VAREN:BoostCalToolbox', ...
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
    ready=all(isfinite(X),2) & isfinite(z) & year(timeLocal)==2020;
    rawRows=find(ready);
    if ~isequal(rawRows,(8762:9481)')
        error('VAREN:BoostCalRows', ...
            'Expected exactly the original 720 calibration rows, raw indices 8762:9481.');
    end
    Xcal=X(rawRows,:);
    zcal=z(rawRows);
    ycal=loadMW(rawRows);
    previousY=loadMW(rawRows-1);
    folds=chronological_folds(numel(rawRows),cfg);
    nStructures=size(cfg.StructureGrid,1);
    nCycles=numel(cfg.CycleGrid);
    scores=NaN(nStructures,nCycles);
    foldMSE=NaN(nStructures,cfg.K,nCycles);
    foldFitSeconds=NaN(nStructures,cfg.K);
    foldPredictSeconds=NaN(nStructures,cfg.K);
    foldTotalSeconds=NaN(nStructures,cfg.K);
    structureSeconds=NaN(nStructures,1);
    prefixCheck=struct();
    calibrationClock=tic;
    for j=1:nStructures
        structureClock=tic;
        [score,foldLoss,check,fitTime,predictTime,totalTime]= ...
            calibration_cv(Xcal,zcal,ycal,previousY,cfg.StructureGrid(j,:),cfg,j==1);
        scores(j,:)=score;
        foldMSE(j,:,:)=reshape(foldLoss,[1,cfg.K,nCycles]);
        foldFitSeconds(j,:)=fitTime(:)';
        foldPredictSeconds(j,:)=predictTime(:)';
        foldTotalSeconds(j,:)=totalTime(:)';
        extraCheckSeconds=0;
        if j==1
            prefixCheck=check;
            extraCheckSeconds=check.Seconds;
        end
        structureSeconds(j)=toc(structureClock)-extraCheckSeconds;
        if cfg.Verbose
            fprintf('LSBoost calibration-only: structure %d/%d; splits=%d, leaf=%d complete.\n', ...
                j,nStructures,cfg.StructureGrid(j,1),cfg.StructureGrid(j,2));
        end
    end
    calibrationSeconds=toc(calibrationClock)-prefixCheck.Seconds;
    bestByStructure=min(scores,[],2);
    selectedStructure=find(bestByStructure==min(bestByStructure),1,'first');
    selectedCycle=find(scores(selectedStructure,:)==min(scores(selectedStructure,:)),1,'first');
    selected=cfg.StructureGrid(selectedStructure,:);

    C=struct();
    C.Protocol='LSBoost bounded pre-forecast calibration audit v1.0';
    C.ForecastOrigins=0;
    C.Scores=scores;
    C.CalibrationTable=make_calibration_table(cfg,scores,selectedStructure,selectedCycle);
    C.FoldMSE=foldMSE;
    C.SelectedStructureIndex=selectedStructure;
    C.SelectedMaxNumSplits=selected(1);
    C.SelectedMinLeafSize=selected(2);
    C.LearnRate=cfg.LearnRate;
    C.SelectedNumLearningCycles=cfg.CycleGrid(selectedCycle);
    C.SelectedCycleAtLowerBoundary=selectedCycle==1;
    C.SelectedCycleAtUpperBoundary=selectedCycle==nCycles;
    C.CalibrationSeconds=calibrationSeconds;
    C.FoldFitSeconds=foldFitSeconds;
    C.FoldPredictSeconds=foldPredictSeconds;
    C.FoldTotalSeconds=foldTotalSeconds;
    C.StructureSeconds=structureSeconds;
    C.FittedMaximumCycles=max(cfg.CycleGrid);

    A=struct();
    A.Config=cfg;
    A.InputTimeUTC=timeUTC;
    A.InputLoadMW=loadMW;
    A.InputIsObserved=observed;
    A.PredictorNames=predictorNames;
    A.CalibrationRawIndex=rawRows;
    A.CalibrationStartUTC=timeUTC(rawRows(1));
    A.CalibrationEndUTC=timeUTC(rawRows(end));
    A.FirstForecastUTC=firstForecast;
    A.Folds=folds;
    A.PrefixCheck=prefixCheck;
    A.MATLABVersion=version;
    A.ToolboxVersions=ver('stats');
    A.ComputerArchitecture=computer;
    A.Scope=['Calibration-only search using exactly 720 model-ready rows before ' ...
        'the first evaluation target. No evaluation observations, forecasts, ' ...
        'errors or rolling evaluation windows are supplied to this function. ' ...
        'All four original structures and LearnRate=.1 are retained. ' ...
        'The sole bounded extension is CycleGrid=[25 50 100 200 400 800 1600]. ' ...
        'Timing measures calibration on these folds, not a full-period runtime.'];
    C.TotalSeconds=toc(coreClock);
end

function cfg=fixed_options(opts)
    allowed={'IsObserved','Verbose'};
    names=fieldnames(opts);
    for j=1:numel(names)
        if ~ismember(names{j},allowed)
            error('VAREN:BoostCalOption', ...
                'Unknown or non-editable calibration option: %s.',names{j});
        end
    end
    cfg=struct('InputRows',9481,'CalibrationRows',720, ...
        'CalibrationFirstRawIndex',8762,'CalibrationLastRawIndex',9481, ...
        'K',5,'InitialTrainFraction',.5, ...
        'StructureGrid',[3 20;3 5;15 20;15 5], ...
        'LearnRate',.1,'CycleGrid',[25 50 100 200 400 800 1600], ...
        'RandomSeed',451,'Verbose',true);
    if isfield(opts,'Verbose'), cfg.Verbose=opts.Verbose; end
    if ~(isscalar(cfg.Verbose) && ...
            (islogical(cfg.Verbose) || (isnumeric(cfg.Verbose) && ismember(cfg.Verbose,[0 1]))))
        error('VAREN:BoostCalVerbose','Verbose must be a logical scalar or 0/1.');
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

function [scores,foldMSE,check,fitSeconds,predictSeconds,totalSeconds]= ...
        calibration_cv(X,z,y,previousY,structure,cfg,doPrefixCheck)
    folds=chronological_folds(size(X,1),cfg);
    nCycles=numel(cfg.CycleGrid);
    squaredSum=zeros(1,nCycles);
    count=0;
    foldMSE=NaN(cfg.K,nCycles);
    fitSeconds=NaN(cfg.K,1);
    predictSeconds=NaN(cfg.K,1);
    totalSeconds=NaN(cfg.K,1);
    check=struct();
    for k=1:cfg.K
        foldClock=tic;
        training=folds.TrainStart(k):folds.TrainEnd(k);
        validation=folds.ValidationStart(k):folds.ValidationEnd(k);
        if max(training)>=min(validation)
            error('VAREN:BoostCalCausality','Chronological training must precede validation.');
        end
        fitClock=tic;
        mdl=fit_ensemble(X(training,:),z(training),structure,max(cfg.CycleGrid),cfg);
        fitSeconds(k)=toc(fitClock);
        predictionClock=tic;
        zhat=NaN(numel(validation),nCycles);
        for j=1:nCycles
            zhat(:,j)=predict_checked(mdl,X(validation,:),cfg.CycleGrid(j));
        end
        predictSeconds(k)=toc(predictionClock);
        extraCheckSeconds=0;
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
            extraCheckSeconds=check.Seconds;
            if ~check.Pass
                error('VAREN:BoostCalPrefixPrediction', ...
                    'Staged prediction differs from separately fitted prefix: %.17g > %.17g.', ...
                    discrepancy,tolerance);
            end
        end
        predicted=previousY(validation).*exp(zhat/100);
        if any(~isfinite(predicted(:))) || any(predicted(:)<=0)
            error('VAREN:BoostCalCVForecast','Non-finite or nonpositive validation prediction.');
        end
        losses=(y(validation)-predicted).^2;
        foldMSE(k,:)=mean(losses,1);
        squaredSum=squaredSum+sum(losses,1);
        count=count+numel(validation);
        totalSeconds(k)=toc(foldClock)-extraCheckSeconds;
        if cfg.Verbose
            fprintf('  Calibration fold %d/%d; train=%d, validate=%d, path=%d trees.\n', ...
                k,cfg.K,numel(training),numel(validation),max(cfg.CycleGrid));
        end
    end
    scores=squaredSum/count;
    if any(~isfinite(scores)) || any(~isfinite(foldMSE(:)))
        error('VAREN:BoostCalCVLoss','Non-finite chronological calibration loss.');
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

