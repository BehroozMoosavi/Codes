function [runDir, zipPath] = run_varen_step5_full(dataFile, pilotLocation)
%RUN_VAREN_STEP5_FULL Validate resume, then run all 34343 forecast origins.
%
% NEW RUN:
%   [step5Folder,step5Zip] = run_varen_step5_full(step3DataFile,step4Folder)
%   [step5Folder,step5Zip] = run_varen_step5_full() % select two MAT files
%
% RESUME AN EXISTING RUN (even after restarting MATLAB):
%   [step5Folder,step5Zip] = run_varen_step5_full('/full/path/to/step5_full_...')
%
% Put this file and varen_checkpointed_core.m in the MATLAB current folder.
% The original tested varen_corrected_core.m is not modified or overwritten.
% The full tuning grids, chronological folds and forecast rules are retained.
% A 12-origin pause / 24-origin resume check must reproduce the supplied
% Step 4 predictions and parameters before the full run begins.
%
% Revision 1.0. Source-reviewed; MATLAB execution is verified by the built-in
% equivalence gate on the user's machine. No execution elsewhere is claimed.
% Full-run checkpoints are saved every 250 completed origins and at the end.
% Use only one MATLAB session for a given output folder. Keep the source
% files unchanged while a run is in progress or will be resumed.

    corePath = which('varen_checkpointed_core');
    if isempty(corePath) || exist('varen_checkpointed_core','file')~=2
        error('VAREN:Step5Core','Put varen_checkpointed_core.m beside this wrapper.');
    end
    wrapperPath = [mfilename('fullpath') '.m'];
    coreSource = fileread(corePath);
    resumeRun = nargin==1 && isfolder(dataFile) && ...
        isfile(fullfile(dataFile,'step5_inputs.mat'));

    if resumeRun
        runDir = char(dataFile);
        B = load(fullfile(runDir,'step5_inputs.mat'), ...
            'D','reference','coreSourceAtStart','inputPath','pilotMatPath');
        must(isequal(B.coreSourceAtStart,coreSource), ...
            'The checkpointed core changed. Restore the source used for this run.');
        D = B.D; reference = B.reference;
        inputPath = B.inputPath; pilotMatPath = B.pilotMatPath;
    else
        if nargin==1 && isfolder(dataFile)
            error('VAREN:Step5Resume','This folder has no step5_inputs.mat to resume.');
        end
        if nargin<1 || isempty(dataFile)
            dataFile = chooseMat('Select the Step 3 caiso_oasis_inputs_2019_2023.mat');
        end
        if nargin<2 || isempty(pilotLocation)
            pilotLocation = chooseMat('Select varen_corrected_results.mat from the Step 4 folder');
        end
        inputPath = char(dataFile);
        if isfolder(pilotLocation)
            pilotMatPath = fullfile(pilotLocation,'varen_corrected_results.mat');
        else
            pilotMatPath = char(pilotLocation);
        end
        must(isfile(inputPath) && isfile(pilotMatPath),'Input or Step 4 reference MAT file is missing.');
        D = load(inputPath,'timeUTC','intervalEndUTC','loadMW','isObserved', ...
            'opts','protocol','sourceMetadata');
        reference = load(pilotMatPath,'R','A');
        validateInputs(D,reference);
        baseDir = fileparts(wrapperPath);
        [~,uniqueName] = fileparts(tempname(baseDir));
        runDir = fullfile(baseDir,['step5_full_' uniqueName]);
        [ok,msg] = mkdir(runDir);
        if ~ok, error('VAREN:Step5Folder','Cannot create output folder: %s',msg); end
        coreSourceAtStart = coreSource;
        save(fullfile(runDir,'step5_inputs.mat'),'D','reference', ...
            'coreSourceAtStart','inputPath','pilotMatPath','-v7.3');
        copyfile(wrapperPath,fullfile(runDir,'run_varen_step5_full.m'));
        copyfile(corePath,fullfile(runDir,'varen_checkpointed_core.m'));
    end

    zipPath = fullfile(runDir,'step5_full_results.zip');
    oldDiaryFile = get(0,'DiaryFile'); oldDiaryState = get(0,'Diary');
    diaryCleanup = onCleanup(@()restoreDiary(oldDiaryFile,oldDiaryState)); %#ok<NASGU>
    diary(fullfile(runDir,'step5_console.txt'));
    boundaryState = warning('query','VAREN:LambdaBoundary');
    warningCleanup = onCleanup(@()warning(boundaryState)); %#ok<NASGU>
    warning('on','VAREN:LambdaBoundary');

    try
        fprintf('\nVA-REN REVISION: STEP 5 FULL EVALUATION\nRevision 1.0; MATLAB %s\n',version);
        fprintf('Existing run resumed: %d\nRun folder: %s\n',resumeRun,runDir);
        fprintf('Input source: %s\nStep 4 reference: %s\n',inputPath,pilotMatPath);
        validateInputs(D,reference);
        opts = struct('IsObserved',D.isObserved,'VolCalibrationYear',2019, ...
            'FirstEvaluationYear',2020,'LastEvaluationYear',2023,'Verbose',true);

        gatePath = fullfile(runDir,'step5_equivalence_check.mat');
        if isfile(gatePath)
            gateData = load(gatePath,'G'); G = gateData.G;
            must(strcmp(G.Status,'PASS') && isequal(G.CoreSource,coreSource) && ...
                strcmp(G.MATLABVersion,version),'The saved equivalence gate does not match this run.');
            fprintf('Checkpoint/resume equivalence: prior PASS retained.\n');
        else
            fprintf('\nChecking resume using the same 24 Step 4 forecast origins.\n');
            checkDir = fullfile(runDir,'resume_check');
            opts.OutputDir = checkDir; opts.MaxForecasts = 24;
            phasePath = fullfile(runDir,'step5_resume_phase1.mat');
            if ~isfile(phasePath)
                control = struct('CheckpointEvery',12,'StopAfter',12);
                [~,~,phase1] = varen_checkpointed_core(D.timeUTC,D.loadMW,opts,control);
                must(phase1.CompletedOrigins==12 && strcmp(phase1.Status,'paused'), ...
                    'Resume test did not pause after 12 complete origins.');
                savedPhase = load(phase1.CheckpointFile,'checkpoint');
                phase1Prefix = struct('Predicted',savedPhase.checkpoint.Predicted, ...
                    'Seconds',savedPhase.checkpoint.Seconds, ...
                    'AlphaCalibrationSeconds',savedPhase.checkpoint.AlphaCalibrationSeconds);
                save(phasePath,'phase1','phase1Prefix','-v7');
            else
                phaseData = load(phasePath,'phase1','phase1Prefix');
                must(phaseData.phase1.CompletedOrigins==12,'Invalid saved resume-test phase.');
                phase1Prefix = phaseData.phase1Prefix;
            end
            control = struct('CheckpointEvery',12,'StopAfter',Inf);
            [Rp,Ap,phase2] = varen_checkpointed_core(D.timeUTC,D.loadMW,opts,control);
            must(phase2.Resumed && phase2.StartedFromCompletedOrigins>=12 && ...
                phase2.CompletedOrigins==24 && strcmp(phase2.Status,'complete'), ...
                'The 24-origin continuation did not use the saved checkpoint.');
            comparisons = comparePilot(reference.R,reference.A,Rp,Ap);
            compareNumbers(phase1Prefix.Predicted,Rp.PredictionMW(1:12,:),0,'Retained forecast prefix');
            compareNumbers(phase1Prefix.Seconds,Rp.RuntimeSeconds(1:12,:),0,'Retained prefix runtimes');
            compareNumbers(phase1Prefix.AlphaCalibrationSeconds,Rp.StaticAlphaCalibrationSeconds, ...
                0,'Retained one-time calibration');
            comparisons = [comparisons; table( ...
                {'RetainedPrefixPredictions';'RetainedPrefixRuntimes';'RetainedCalibrationTime'}, ...
                zeros(3,1),zeros(3,1),true(3,1),'VariableNames',comparisons.Properties.VariableNames)];
            G = struct('Status','PASS','MATLABVersion',version,'CoreSource',coreSource, ...
                'ReferenceForecastOrigins',24,'PauseAfterOrigins',12, ...
                'PredictionToleranceMW',1e-7,'Comparisons',comparisons,'ResumeState',phase2);
            writetable(comparisons,fullfile(runDir,'step5_equivalence_checks.csv'));
            save(gatePath,'G','-v7');
            fprintf('CHECKPOINT_EQUIVALENCE: PASS (12 paused, 24 completed, Step 4 reproduced).\n');
        end

        fprintf('\nStarting/continuing all 34343 forecast origins with the fixed settings.\n');
        fprintf('Progress is printed by the core; checkpoints save every 250 complete origins.\n');
        fprintf('The pilot suggested about 5.3 hours of serial fitting; actual time may differ.\n');
        fprintf('To resume, run run_varen_step5_full with the run folder printed above.\n\n');
        fullDir = fullfile(runDir,'full');
        opts.OutputDir = fullDir; opts.MaxForecasts = Inf;
        control = struct('CheckpointEvery',250,'StopAfter',Inf);
        callClock = tic;
        [R,A,S] = varen_checkpointed_core(D.timeUTC,D.loadMW,opts,control);
        fullCallWallSeconds = toc(callClock);
        validateFull(D,reference,R,A,S);
        save(fullfile(runDir,'step5_execution_state.mat'),'S','fullCallWallSeconds','-v7');

        endpointReport = table(R.ModelNames(3:7)', ...
            sum(R.LambdaAtGridBoundary(:,3:7),1)', ...
            mean(R.LambdaAtGridBoundary(:,3:7),1)', ...
            sum(R.Lambda(:,3:7)==min(A.Config.LambdaGrid),1)', ...
            sum(R.Lambda(:,3:7)==max(A.Config.LambdaGrid),1)', ...
            'VariableNames',{'Model','EndpointSelections','EndpointFraction', ...
            'LowerEndpointSelections','UpperEndpointSelections'});
        writetable(endpointReport,fullfile(runDir,'step5_lambda_endpoints.csv'));
        originSeconds = sum(R.RuntimeSeconds,2);
        timing = table(R.ModelNames(:),sum(R.RuntimeSeconds,1)', ...
            mean(R.RuntimeSeconds,1)',median(R.RuntimeSeconds,1)', ...
            min(R.RuntimeSeconds,[],1)',max(R.RuntimeSeconds,[],1)', ...
            'VariableNames',{'Model','FullForecastLoopSeconds','MeanSecondsPerOrigin', ...
            'MedianSecondsPerOrigin','MinSecondsPerOrigin','MaxSecondsPerOrigin'});
        writetable(timing,fullfile(runDir,'step5_timing_by_model.csv'));
        timestampCheck = table(R.RawIndex,utcText(R.TimeUTC), ...
            utcText(D.intervalEndUTC(R.RawIndex)),R.TargetMW, ...
            'VariableNames',{'RawIndex','IntervalStartUTC','IntervalEndUTC','TargetMW'});
        writetable(timestampCheck,fullfile(runDir,'step5_forecast_intervals.csv'));

        summary = sprintf([ ...
            'VA-REN REVISION: STEP 5 FULL EVALUATION\nRevision 1.0; MATLAB %s\n\n' ...
            'STEP5_STATUS: PASS\nCHECKPOINT_EQUIVALENCE: PASS\n' ...
            'FULL_PERIOD_FORECASTS_RUN: YES\nFORECAST_ORIGINS: 34343\nCONFIGURED_OUTPUTS: 7\n\n' ...
            'First forecast interval start: %s\nLast forecast interval start: %s\n' ...
            'All origins, targets, naive forecasts and model parameter domains passed validation.\n' ...
            'The first 24 full-run forecasts also reproduce the original Step 4 pilot.\n' ...
            'Selected static alpha: %.6g (selected once using the initial calibration only).\n' ...
            'LASSO and calibrated StaticEN coincide; FixedEN05 remains the fixed-alpha=.5 control.\n\n' ...
            'One-time static-alpha calibration: %.3f seconds\n' ...
            'Total recorded per-model forecast-loop fitting time: %.3f seconds (%.3f hours)\n' ...
            'Mean / median origin fitting time: %.3f / %.3f seconds\n' ...
            'Core cumulative compute time: %.3f seconds\n' ...
            'Wall time of the latest full-core invocation: %.3f seconds\n' ...
            'Latest invocation started from %d saved origins; computed %d additional origins.\n' ...
            'Latest invocation checkpoint write time: %.3f seconds\n' ...
            'Latest invocation checkpoint read time: %.3f seconds\n' ...
            'Run pauses are excluded from core compute/fit timings.\n' ...
            'Compute totals cover retained work; a process crash can leave abandoned work after the last save unrecorded.\n' ...
            'Initial resume-equivalence check and final ZIP packaging are outside full-run fitting times.\n\n' ...
            'Lambda endpoint counts/fractions are retained for review; no grid changes were made.\n' ...
            'Target: hourly system demand, not solar-subtracted net load.\n' ...
            'Historical published actuals may contain revisions; this is not an as-of vintage archive.\n' ...
            'Completed-hour load is assumed available at the next forecast origin.\n' ...
            'Annual lag: 8760 elapsed hours. Calendar predictors: America/Los_Angeles.\n' ...
            'No DM/GW significance tests or nonlinear benchmark are added by this run.\n' ...
            'A successful run is not a finding of accuracy superiority.\n' ...
            'Return step5_full_results.zip for the next results/revision review.\n'], ...
            version,char(utcText(R.TimeUTC(1))),char(utcText(R.TimeUTC(end))),R.StaticAlpha, ...
            R.StaticAlphaCalibrationSeconds,sum(originSeconds),sum(originSeconds)/3600, ...
            mean(originSeconds),median(originSeconds),R.TotalSeconds,fullCallWallSeconds, ...
            S.StartedFromCompletedOrigins,S.OriginsComputedThisCall, ...
            S.CurrentCallCheckpointSeconds,S.CurrentCallResumeReadSeconds);
        writeText(fullfile(runDir,'step5_summary.txt'),summary);
        fprintf('\n%s\n',summary);
        disp(R.Metrics); disp(endpointReport);
        diary('off');
        packageResults(runDir,zipPath,false);
        fprintf('Upload ZIP: %s\n',zipPath);
    catch ME
        report = getReport(ME,'extended','hyperlinks','off');
        fprintf('\nSTEP5_STATUS: FAIL\n%s\nRun folder retained: %s\n',report,runDir);
        try
            writeText(fullfile(runDir,'step5_summary.txt'),sprintf([ ...
                'STEP5_STATUS: FAIL\nCompleted-origin checkpoints are retained locally.\n' ...
                'Do not change the grid or suppress a fitting warning to force completion.\n\n%s\n'],report));
            diary('off');
            failureZip = fullfile(runDir,'step5_failure.zip');
            packageResults(runDir,failureZip,true);
            fprintf('Upload failure ZIP: %s\n',failureZip);
        catch packagingError
            diary('off');
            fprintf('Could not package the failure: %s\nReturn step5_console.txt and the error.\n', ...
                packagingError.message);
        end
        rethrow(ME);
    end
end

function validateInputs(D,reference)
    must(all(isfield(D,{'timeUTC','intervalEndUTC','loadMW','isObserved','opts', ...
        'protocol','sourceMetadata'})) && all(isfield(reference,{'R','A'})), ...
        'Select the Step 3 input MAT and the Step 4 core-result MAT.');
    must(isdatetime(D.timeUTC) && strcmp(D.timeUTC.TimeZone,'UTC') && ...
        iscolumn(D.timeUTC) && numel(D.timeUTC)==43824 && ...
        numel(D.loadMW)==43824 && numel(D.intervalEndUTC)==43824 && ...
        iscolumn(D.loadMW) && iscolumn(D.intervalEndUTC) && ...
        all(diff(D.timeUTC)==hours(1)) && all(D.intervalEndUTC-D.timeUTC==hours(1)), ...
        'Unexpected prepared physical-hour series.');
    must(islogical(D.isObserved) && all(D.isObserved) && ...
        all(isfinite(D.loadMW)) && all(D.loadMW>0) && ...
        D.timeUTC(1)==datetime(2019,1,1,8,0,0,'TimeZone','UTC') && ...
        D.intervalEndUTC(end)==datetime(2024,1,1,8,0,0,'TimeZone','UTC'), ...
        'Invalid input values, provenance flags or series coverage.');
    must(strcmp(D.sourceMetadata.canonical_csv_sha256, ...
        '249120565acab4af01f3dab853da075562b809edaf758ff0f7f7033774b33881') && ...
        strcmp(D.protocol.TimestampProvenance,'provider_supplied_UTC_interval_start_and_end') && ...
        D.protocol.ExpectedForecasts==34343 && D.protocol.VolatilityCalibrationWindows==8711, ...
        'Input provenance or forecast protocol differs from the reviewed data.');
    must(isequal(reference.A.InputTimeUTC,D.timeUTC) && ...
        isequal(reference.A.InputLoadMW,D.loadMW) && ...
        isequal(reference.A.InputIsObserved,D.isObserved), ...
        'The Step 4 reference used different input data.');
    must(isequal(reference.R.RawIndex,(9482:9505)') && ...
        isequal(size(reference.R.PredictionMW),[24 7]) && ...
        isequal(reference.R.ModelNames,modelNames()) && reference.R.StaticAlpha==1 && ...
        reference.A.Config.MaxForecasts==24,'Unexpected Step 4 reference run.');
    expected = struct('W',720,'Wv',48,'K',5,'CalibrationRows',720, ...
        'EvaluationWarmupRows',720,'InitialTrainFraction',.5,'AlphaFloor',.001, ...
        'FixedAlpha',.5,'AlphaGrid',[.05 .1 .2 .3 .4 .5 .6 .7 .8 .9 1], ...
        'LambdaGrid',logspace(-5,2,60),'VolCalibrationYear',2019, ...
        'FirstEvaluationYear',2020,'LastEvaluationYear',2023,'RelTol',1e-6,'MaxIter',100000);
    fields = fieldnames(expected);
    for j=1:numel(fields)
        f=fields{j};
        must(isequal(reference.A.Config.(f),expected.(f)),['Unexpected Step 4 setting: ' f]);
    end
end

function checks = comparePilot(R0,A0,R1,A1)
    must(isequal(R0.RawIndex,R1.RawIndex) && isequal(R0.TimeUTC,R1.TimeUTC) && ...
        isequal(R0.TimeLocal,R1.TimeLocal) && isequal(R0.TargetMW,R1.TargetMW) && ...
        isequal(R0.ModelNames,R1.ModelNames),'Resume test changed forecast alignment.');
    fields = fieldnames(A0.Config);
    for j=1:numel(fields)
        f=fields{j};
        if strcmp(f,'OutputDir'), continue; end
        must(isequal(A0.Config.(f),A1.Config.(f)),['Resume test changed setting: ' f]);
    end
    names = {'PredictionMW','ErrorMW','Alpha','Lambda','EffectiveL1', ...
        'EffectiveL2','ActiveSetSize','StaticAlpha'};
    tolerances = [1e-7 1e-7 1e-12 1e-12 1e-12 1e-12 0 0];
    label = cell(0,1); maxDifference = zeros(0,1); tolerance = zeros(0,1);
    for j=1:numel(names)
        label{end+1,1}=['R.' names{j}]; %#ok<AGROW>
        maxDifference(end+1,1)=compareNumbers(R0.(names{j}),R1.(names{j}),tolerances(j),label{end}); %#ok<AGROW>
        tolerance(end+1,1)=tolerances(j); %#ok<AGROW>
    end
    names = {'ForecastVolatility','ForecastNormalizedVolatility','ForecastUnclippedVolatility', ...
        'ForecastAlpha','VolatilityCalibrationMin','VolatilityCalibrationMax','VolatilityCalibrationCount'};
    for j=1:numel(names)
        label{end+1,1}=['A.' names{j}]; %#ok<AGROW>
        maxDifference(end+1,1)=compareNumbers(A0.(names{j}),A1.(names{j}),1e-10,label{end}); %#ok<AGROW>
        tolerance(end+1,1)=1e-10; %#ok<AGROW>
    end
    compareNumbers(R0.AlphaCalibrationTable.Alpha,R1.AlphaCalibrationTable.Alpha,0,'Calibration alpha candidates');
    delta=compareNumbers(R0.AlphaCalibrationTable.MinCVLevelMSE, ...
        R1.AlphaCalibrationTable.MinCVLevelMSE,1e-7,'Initial calibration CV losses');
    label{end+1,1}='StaticAlphaCalibrationCVLoss';
    maxDifference(end+1,1)=delta; tolerance(end+1,1)=1e-7;
    must(isequal(R0.LambdaAtGridBoundary,R1.LambdaAtGridBoundary), ...
        'Resume test changed lambda endpoint flags.');
    for f={'VolatilityCalibrationStartUTC','VolatilityCalibrationEndUTC', ...
            'StaticCalibrationStartUTC','StaticCalibrationEndUTC'}
        must(isequal(A0.(f{1}),A1.(f{1})),['Resume test changed ' f{1}]);
    end
    checks=table(label,maxDifference,tolerance,true(numel(label),1), ...
        'VariableNames',{'Quantity','MaxAbsoluteDifference','AbsoluteTolerance','Passed'});
end

function validateFull(D,reference,R,A,S)
    N=34343; rows=(9482:43824)';
    must(strcmp(S.Status,'complete') && S.CompletedOrigins==N && S.TotalOrigins==N, ...
        'Full core did not report all 34343 completed origins.');
    must(isequal(R.ModelNames,modelNames()) && isequal(R.RawIndex,rows) && ...
        isequal(size(R.PredictionMW),[N 7]) && isequal(R.TimeUTC,D.timeUTC(rows)) && ...
        isequal(R.TargetMW,D.loadMW(rows)),'Full-run origins, models or targets are misaligned.');
    must(all(isfinite(R.PredictionMW(:))) && all(R.PredictionMW(:)>0) && ...
        isequal(R.ErrorMW,R.TargetMW-R.PredictionMW) && ...
        isequal(R.PredictionMW(:,1),D.loadMW(rows-1)) && ...
        isequal(R.PredictionMW(:,2),D.loadMW(rows-168)), ...
        'Invalid full-run predictions, errors or naive-reference alignment.');
    compareNumbers(R.PredictionMW(1:24,:),reference.R.PredictionMW,1e-7,'Full-run first 24 forecasts');
    compareNumbers(R.Alpha(1:24,:),reference.R.Alpha,1e-12,'Full-run first 24 alphas');
    compareNumbers(R.Lambda(1:24,:),reference.R.Lambda,1e-12,'Full-run first 24 lambdas');
    must(isequal(A.InputTimeUTC,D.timeUTC) && isequal(A.InputLoadMW,D.loadMW) && ...
        isequal(A.InputIsObserved,D.isObserved),'Full-run inputs changed.');
    fields=fieldnames(reference.A.Config);
    for j=1:numel(fields)
        f=fields{j};
        if any(strcmp(f,{'OutputDir','MaxForecasts'})), continue; end
        must(isequal(A.Config.(f),reference.A.Config.(f)),['Full-run setting changed: ' f]);
    end
    must(isinf(A.Config.MaxForecasts) && R.StaticAlpha==reference.R.StaticAlpha && ...
        A.VolatilityCalibrationCount==8711 && numel(A.PredictorNames)==70, ...
        'Unexpected full-run calibration or predictor configuration.');
    compareNumbers(R.PredictionMW(:,4),R.PredictionMW(:,5),1e-7,'LASSO/StaticEN identity at alpha 1');
    for f={'Alpha','Lambda','ActiveSetSize','RuntimeSeconds','LambdaAtGridBoundary'}
        must(isequal(size(R.(f{1})),[N 7]),['Wrong output dimensions: ' f{1}]);
    end
    a=R.Alpha(:,3:7); lam=R.Lambda(:,3:7); active=R.ActiveSetSize(:,3:7);
    must(all(isfinite(a(:))) && all(a(:)>=0 & a(:)<=1) && ...
        all(isfinite(lam(:))) && all(ismember(lam(:),A.Config.LambdaGrid)) && ...
        all(isfinite(active(:))) && all(active(:)>=0 & active(:)<=70) && ...
        all(active(:)==floor(active(:))),'Invalid full-run fitted parameters.');
    must(all(R.Alpha(:,3)==0) && all(R.Alpha(:,4)==1) && all(R.Alpha(:,5)==1) && ...
        all(R.Alpha(:,6)==.5),'Fixed model alphas changed.');
    compareNumbers(R.Alpha(:,7),A.ForecastAlpha,0,'Saved VA alpha identity');
    compareNumbers(A.ForecastAlpha,max(.001,1-min(1,max(0,A.ForecastUnclippedVolatility))), ...
        1e-12,'VA alpha mapping');
    must(isequal(R.LambdaAtGridBoundary(:,3:7), ...
        lam==min(A.Config.LambdaGrid) | lam==max(A.Config.LambdaGrid)), ...
        'Endpoint flags disagree with selected penalties.');
    must(all(isfinite(R.RuntimeSeconds(:))) && all(R.RuntimeSeconds(:)>=0) && ...
        height(R.Metrics)==7 && all(R.Metrics.N==N) && ...
        all(isfinite(R.Metrics.RMSE_MW)),'Invalid full-run timing or accuracy summary.');
end

function delta = compareNumbers(a,b,tol,label)
    must(isequal(size(a),size(b)) && isequal(isnan(a),isnan(b)), ...
        [label ': dimensions or NaN pattern differ.']);
    keep=~isnan(a);
    must(all(isfinite(a(keep))) && all(isfinite(b(keep))),[label ': non-finite values.']);
    if ~any(keep(:)), delta=0; else, delta=max(abs(a(keep)-b(keep))); end
    must(delta<=tol,sprintf('%s: difference %.12g exceeds %.12g.',label,delta,tol));
end

function names = modelNames()
    names={'Persistence','SeasonalNaive','Ridge','LASSO','StaticEN','FixedEN05','VAREN'};
end

function file = chooseMat(titleText)
    [name,folder]=uigetfile('*.mat',titleText);
    if isequal(name,0), error('VAREN:Step5Cancelled','No MAT file was selected.'); end
    file=fullfile(folder,name);
end

function must(condition,message)
    if ~isscalar(condition) || ~condition
        error('VAREN:Step5Validation','%s',message);
    end
end

function s = utcText(t)
    t.TimeZone='UTC'; t.Format='yyyy-MM-dd HH:mm:ss';
    s=replace(string(t),' ','T')+"Z";
end

function writeText(path,value)
    [fid,msg]=fopen(path,'w','n','UTF-8');
    if fid<0, error('VAREN:Step5Write','Cannot write report: %s',msg); end
    closeFile=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s',value);
end

function packageResults(folder,zipPath,isFailure)
    listing=dir(folder); names={listing(~[listing.isdir]).name};
    names=names(~endsWith(names,'.zip','IgnoreCase',true));
    for part={'resume_check','full'}
        sub=fullfile(folder,part{1});
        if ~isfolder(sub), continue; end
        files=dir(sub); leaf={files(~[files.isdir]).name};
        resultNames={'varen_corrected_results.mat','aligned_predictions.csv', ...
            'parameters_and_runtime.csv','accuracy_common_sample.csv', ...
            'accuracy_by_clock_window.csv','volatility_diagnostics.csv','static_alpha_calibration.csv'};
        for j=1:numel(leaf)
            n=leaf{j};
            isCheckpoint=strcmp(n,'varen_checkpoint.mat') || strcmp(n,'varen_checkpoint_previous.mat');
            if ~any(strcmp(n,resultNames)) && ~(isCheckpoint && isFailure), continue; end
            names{end+1}=fullfile(part{1},n); %#ok<AGROW>
        end
    end
    zip(zipPath,names,folder);
end

function restoreDiary(oldFile,oldState)
    diary('off');
    if strcmp(oldState,'on'), diary(oldFile); end
end
