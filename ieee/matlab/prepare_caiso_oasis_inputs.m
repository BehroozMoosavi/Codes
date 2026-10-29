function [dataFile, reportDir] = prepare_caiso_oasis_inputs(dataDir)
%PREPARE_CAISO_OASIS_INPUTS Step 3: validate provider-stamped physical hours.
% [dataFile,reportDir] = prepare_caiso_oasis_inputs(dataDir)
% dataDir defaults to pwd and must contain 451_step3_oasis_source.zip.
% No model is fitted. Original EMS spreadsheets are not used or altered.
%
% The supplied archive contains a canonical CSV from CAISO OASIS SLD_FCST,
% ACTUAL, CA ISO-TAC, SYS_FCST_ACT_MW, plus its acquisition manifest/metadata.
% Both interval endpoints are provider-supplied UTC, not reconstructed from
% local hour labels. SourceArchive/SourceMember/SourceRow identify each value.
% A separate source audit compares the canonical values with raw replies.
% This helper checks the delivered series and manifest; it does not refetch
% OASIS or perform that full comparison with every raw response member.
%
% This is a NEW consistent system-load source series, not a reproduction of
% the submitted EMS experiment or a recovery of averaged 2019 EMS hours.
% Published historical actuals can include revisions; these are not as-of
% data vintages. The subsequent backtest assumes each completed hour's load
% is available at the next forecast origin. The target is system demand,
% not solar-subtracted net load. No interpolation or deduplication is done.
%
% Protocol: local 2019 volatility calibration; local 2020-2023 evaluation
% era; elapsed 8760-hour annual lag; 720 initial model-ready calibration rows.
% Expected: 43824 input rows, 35063 model-ready rows, 34343 forecast origins.
%
% Prepared for MATLAB R2025a; source-reviewed, not executed in MATLAB in the
% authoring environment. Return the validation ZIP before fitting forecasts.
% Revision 1.0, 2026-10-27.

    if nargin<1 || isempty(dataDir), dataDir=pwd; end
    dataDir=char(string(dataDir));
    inputZip=fullfile(dataDir,'451_step3_oasis_source.zip');
    if ~isfile(inputZip)
        error('VAREN:MissingArchive','Place 451_step3_oasis_source.zip in %s.',dataDir);
    end
    fprintf('\nSTEP3: validate CAISO OASIS provider UTC intervals. No forecasts.\n');
    started=tic;
    [~,suffix]=fileparts(tempname);
    reportDir=fullfile(pwd,['step3_oasis_' suffix]);
    if isfolder(reportDir) || isfile(reportDir)
        error('VAREN:OutputExists','Output already exists: %s',reportDir);
    end
    [ok,msg]=mkdir(reportDir);
    if ~ok, error('VAREN:OutputFolder','%s',msg); end
    sourceDir=fullfile(reportDir,'source');
    [ok,msg]=mkdir(sourceDir);
    if ~ok, error('VAREN:OutputFolder','%s',msg); end
    unzip(inputZip,sourceDir);
    required=["caiso_oasis_actual_2019_2023.csv";"source_manifest.csv"; ...
        "source_metadata.json";"README_SOURCE.txt"];
    if any(~isfile(fullfile(string(sourceDir),required)))
        error('VAREN:ArchiveLayout','The supplied archive is missing required root files.');
    end

    dataNames={'IntervalStartUTC','IntervalEndUTC','LoadMW','Report', ...
        'MarketRunID','TACArea','XMLDataItem','SourceArchive','SourceMember','SourceRow'};
    sourceRows=readTypedCSV(fullfile(sourceDir,required(1)),dataNames,{'LoadMW','SourceRow'});
    manifestNames={'SourceArchive','RequestStartUTC','RequestEndUTC','URL', ...
        'RetrievedUTC','ArchiveSHA256','ArchiveBytes','AllAreasSourceRows', ...
        'SelectedCAISORows','SelectedRowsOutsideRequest'};
    sourceManifest=readTypedCSV(fullfile(sourceDir,required(2)),manifestNames, ...
        {'ArchiveBytes','AllAreasSourceRows','SelectedCAISORows','SelectedRowsOutsideRequest'});
    sourceMetadata=jsondecode(fileread(fullfile(sourceDir,required(3))));
    N=height(sourceRows);
    if N~=43824, error('VAREN:Count','Expected 43824 rows; read %d.',N); end
    timeUTC=parseUTC(sourceRows.IntervalStartUTC,'interval starts');
    intervalEndUTC=parseUTC(sourceRows.IntervalEndUTC,'interval ends');
    first=datetime(2019,1,1,8,0,0,'TimeZone','UTC');
    stop=datetime(2024,1,1,8,0,0,'TimeZone','UTC');
    if timeUTC(1)~=first || intervalEndUTC(end)~=stop || ...
            any(diff(timeUTC)~=hours(1)) || any(intervalEndUTC-timeUTC~=hours(1))
        error('VAREN:UTCGrid','Require exactly the complete, unique, ordered physical-hour grid.');
    end
    loadMW=sourceRows.LoadMW;
    fieldsMatch=sourceRows.Report=="SLD_FCST" & sourceRows.MarketRunID=="ACTUAL" & ...
        sourceRows.TACArea=="CA ISO-TAC" & sourceRows.XMLDataItem=="SYS_FCST_ACT_MW";
    lineagePresent=nonblank(sourceRows.SourceArchive) & nonblank(sourceRows.SourceMember) & ...
        isfinite(sourceRows.SourceRow) & sourceRows.SourceRow>=2 & ...
        sourceRows.SourceRow==fix(sourceRows.SourceRow);
    finitePositive=isfinite(loadMW) & loadMW>0;
    if ~all(fieldsMatch & lineagePresent & finitePositive)
        error('VAREN:Records','Invalid load, report identity, or original-row provenance.');
    end

    % Match every retained row to its independently inventoried request block.
    if height(sourceManifest)~=66 || numel(unique(sourceManifest.SourceArchive))~=66
        error('VAREN:Manifest','Expected 66 uniquely identified source archives.');
    end
    archiveNames="oasis_actual_"+compose('%03d',(1:66)')+".zip";
    if ~isequal(sourceManifest.SourceArchive,archiveNames)
        error('VAREN:ManifestOrder','Unexpected source archive inventory or order.');
    end
    requestStart=parseUTC(sourceManifest.RequestStartUTC,'request starts');
    requestEnd=parseUTC(sourceManifest.RequestEndUTC,'request ends');
    if requestStart(1)~=first || requestEnd(end)~=stop || ...
            any(requestStart(2:end)~=requestEnd(1:end-1)) || any(requestEnd<=requestStart)
        error('VAREN:RequestCoverage','Manifest request intervals do not partition the study period.');
    end
    numericManifest=[sourceManifest.ArchiveBytes sourceManifest.AllAreasSourceRows ...
        sourceManifest.SelectedCAISORows sourceManifest.SelectedRowsOutsideRequest];
    if any(~isfinite(numericManifest(:))) || any(numericManifest(:)<0) || ...
            any(numericManifest(:)~=fix(numericManifest(:))) || ...
            any(sourceManifest.ArchiveBytes==0) || ...
            any(sourceManifest.AllAreasSourceRows<sourceManifest.SelectedCAISORows) || ...
            ~all(nonblank(sourceManifest.URL) & nonblank(sourceManifest.RetrievedUTC)) || ...
            ~all(matchesPattern(sourceManifest.ArchiveSHA256,'^[0-9a-fA-F]{64}$'))
        error('VAREN:ManifestFields','Invalid source manifest provenance fields.');
    end
    [listed,archiveIndex]=ismember(sourceRows.SourceArchive,sourceManifest.SourceArchive);
    if ~all(listed), error('VAREN:Lineage','An input row names an unlisted archive.'); end
    actualCounts=accumarray(archiveIndex,1,[66 1]);
    countsMatch=actualCounts==sourceManifest.SelectedCAISORows & ...
        actualCounts==hours(requestEnd-requestStart);
    withinRequest=timeUTC>=requestStart(archiveIndex) & intervalEndUTC<=requestEnd(archiveIndex);
    if ~all(countsMatch) || ~all(withinRequest)
        error('VAREN:ArchiveCounts','Per-archive row counts or interval bounds disagree.');
    end
    lineageKey=sourceRows.SourceArchive+"|"+sourceRows.SourceMember+"|"+string(sourceRows.SourceRow);
    if numel(unique(lineageKey))~=N
        error('VAREN:LineageDuplicate','The same original source record is reused.');
    end
    verifyMetadata(sourceMetadata,N,height(sourceManifest));
    archiveChecks=table(sourceManifest.SourceArchive,actualCounts, ...
        sourceManifest.SelectedCAISORows,countsMatch, ...
        sourceManifest.SelectedRowsOutsideRequest,'VariableNames', ...
        {'SourceArchive','RetainedRows','ManifestSelectedRows','CountsMatch','ReportedRowsOutsideRequest'});
    % Reported-actual qualification follows the validated source records.
    % Availability at a historical forecast origin is a separate assumption.
    isObserved=finitePositive & fieldsMatch & lineagePresent & listed & withinRequest;
    sourceRows.IsObserved=isObserved;
    if ~all(isObserved), error('VAREN:Observed','Some source rows are not qualified.'); end

    timeLocal=timeUTC; timeLocal.TimeZone='America/Los_Angeles';
    localYear=year(timeLocal);
    localDate=datetime(localYear,month(timeLocal),day(timeLocal));
    [days,~,dayGroup]=unique(localDate);
    dailyCount=accumarray(dayGroup,1);
    midnight=datetime(year(days),month(days),day(days),0,0,0,'TimeZone','America/Los_Angeles');
    nextMidnight=midnight+caldays(1);
    midnight.TimeZone='UTC'; nextMidnight.TimeZone='UTC';
    expectedDaily=hours(nextMidnight-midnight);
    if numel(days)~=1826 || any(dailyCount~=expectedDaily) || ...
            any(~ismember(dailyCount,[23 24 25]))
        error('VAREN:DailyCoverage','Pacific calendar-day coverage check failed.');
    end
    dailyChecks=table(days,dailyCount,expectedDaily,dailyCount==expectedDaily, ...
        'VariableNames',{'LocalDate','ObservedPhysicalHours','ExpectedPhysicalHours','CountMatches'});
    dailyChecks.LocalDate.Format='yyyy-MM-dd';
    years=(2019:2023)'; expectedYearly=[8760;8784;8760;8760;8760];
    yearlyCount=zeros(5,1); days23=zeros(5,1); days24=zeros(5,1); days25=zeros(5,1);
    for j=1:5
        yearlyCount(j)=sum(localYear==years(j)); inYear=year(days)==years(j);
        days23(j)=sum(inYear & dailyCount==23);
        days24(j)=sum(inYear & dailyCount==24);
        days25(j)=sum(inYear & dailyCount==25);
    end
    if any(yearlyCount~=expectedYearly) || any(days23~=1) || any(days25~=1)
        error('VAREN:YearCoverage','Expected annual physical-hour or DST-day counts failed.');
    end
    yearlyChecks=table(years,yearlyCount,expectedYearly,days23,days24,days25, ...
        'VariableNames',{'LocalYear','ObservedPhysicalHours','ExpectedPhysicalHours','Days23Hours','Days24Hours','Days25Hours'});

    % Mirror the core's input-readiness rules; this estimates no model.
    logY=log(loadMW); z=[NaN;100*diff(logY)]; zPast=[NaN;z(1:end-1)];
    V=movvar(zPast,[47 0],0,'omitnan');
    V(movsum(double(isfinite(zPast)),[47 0])~=48)=NaN;
    cal=localYear==2019 & isfinite(V);
    if sum(cal)~=8711 || ~(max(V(cal))>min(V(cal)))
        error('VAREN:VolCalibration','Expected 8711 complete, nonconstant calibration windows.');
    end
    ready=isfinite(z) & isfinite(V) & localYear>=2020 & localYear<=2023;
    for L=[1:24 168 8760]
        ready=ready & [false(L,1);isfinite(z(1:end-L))];
    end
    for L=[1 24 168 8760]
        ready=ready & [false(L,1);isfinite(logY(1:end-L))];
    end
    readyRows=find(ready); expectedForecasts=numel(readyRows)-720;
    if numel(readyRows)~=35063 || any(diff(readyRows)~=1) || expectedForecasts~=34343 || ...
            timeUTC(readyRows(1))~=datetime(2020,1,1,9,0,0,'TimeZone','UTC') || ...
            timeUTC(readyRows(721))~=datetime(2020,1,31,9,0,0,'TimeZone','UTC') || ...
            timeUTC(find(cal,1,'last'))>=timeUTC(readyRows(1))
        error('VAREN:Readiness','Default core readiness, calibration order, or forecast count differs.');
    end
    opts=struct('IsObserved',isObserved,'VolCalibrationYear',2019, ...
        'FirstEvaluationYear',2020,'LastEvaluationYear',2023);
    protocol=struct('Revision','1.0','MATLABVersion',version, ...
        'Source','CAISO OASIS SLD_FCST ACTUAL CA ISO-TAC SYS_FCST_ACT_MW', ...
        'TimestampProvenance','provider_supplied_UTC_interval_start_and_end', ...
        'Target','hourly system demand in MW, not solar-subtracted net load', ...
        'SourceChange','New OASIS series replaces EMS values; not a reproduction or recovery of averaged EMS hours.', ...
        'VintageLimitation','Published historical actuals may contain revisions; not an as-of vintage archive.', ...
        'AvailabilityAssumption','Completed-hour load is assumed available at the next forecast origin.', ...
        'IsObservedMeaning','Finite positive published ACTUAL value with valid report identity and original-row/request provenance.', ...
        'CalendarTimeZone','America/Los_Angeles','AnnualLagHours',8760, ...
        'VolatilityCalibrationYear',2019,'VolatilityWindowHours',48, ...
        'VolatilityCalibrationWindows',sum(cal),'EvaluationYears',2020:2023, ...
        'WarmupRows',720,'InputRows',N,'ModelReadyRows',numel(readyRows), ...
        'ExpectedForecasts',expectedForecasts,'ForecastsRun',false);
    protocol.FirstModelReadyUTC=timeUTC(readyRows(1));
    protocol.FirstExpectedForecastUTC=timeUTC(readyRows(721));
    protocol.LastExpectedForecastUTC=timeUTC(readyRows(end));
    dataFile=fullfile(reportDir,'caiso_oasis_inputs_2019_2023.mat');
    save(dataFile,'timeUTC','intervalEndUTC','loadMW','isObserved','sourceRows', ...
        'sourceManifest','sourceMetadata','opts','protocol','-v7');
    writetable(yearlyChecks,fullfile(reportDir,'yearly_checks.csv'),'Encoding','UTF-8');
    writetable(dailyChecks,fullfile(reportDir,'daily_checks.csv'),'Encoding','UTF-8');
    writetable(archiveChecks,fullfile(reportDir,'archive_checks.csv'),'Encoding','UTF-8');
    summary=sprintf([ ...
        'CAISO OASIS INPUT VALIDATION: STEP 3\nRevision 1.0; MATLAB %s\n\n' ...
        'STEP3_STATUS: PASS\nFORECASTS_RUN: NO\n' ...
        'TIMESTAMP_PROVENANCE: PROVIDER_SUPPLIED_UTC_INTERVAL_START_AND_END\n\n' ...
        '43824 unique, contiguous, finite positive hourly system-load records.\n' ...
        '66 source archives; every retained source row has valid lineage.\n' ...
        'First interval start: %s\nLast interval end: %s\n' ...
        'Every Pacific calendar day has its expected 23, 24, or 25 physical hours.\n' ...
        '2019 causal volatility calibration: 8711 complete 48-hour windows.\n' ...
        '2020-2023 model-ready rows: 35063; initial calibration rows: 720.\n' ...
        'Expected default-core forecasts: 34343.\n' ...
        'First model-ready start: %s (2020-01-01 01:00 Pacific).\n' ...
        'First forecast interval start: %s (2020-01-31 01:00 Pacific).\n\n' ...
        '%s\n%s\n%s\n' ...
        'Target: system demand, not solar-subtracted net load.\n' ...
        'The annual predictor uses 8760 elapsed hours, not a calendar-year match.\n' ...
        'No forecasts, interpolation, deduplication, or EMS repairs were performed.\n' ...
        'Return step3_validation.zip for review before the forecasting step.\n'], ...
        version,char(utcText(timeUTC(1))),char(utcText(intervalEndUTC(end))), ...
        char(utcText(protocol.FirstModelReadyUTC)),char(utcText(protocol.FirstExpectedForecastUTC)), ...
        protocol.SourceChange,protocol.VintageLimitation,protocol.AvailabilityAssumption);
    writeText(fullfile(reportDir,'step3_summary.txt'),summary);
    for j=2:4, copyfile(fullfile(sourceDir,required(j)),fullfile(reportDir,required(j))); end
    zipNames={'step3_summary.txt','yearly_checks.csv','daily_checks.csv','archive_checks.csv', ...
        'source_manifest.csv','source_metadata.json','README_SOURCE.txt','caiso_oasis_inputs_2019_2023.mat'};
    codePath=[mfilename('fullpath') '.m'];
    if isfile(codePath)
        copyfile(codePath,fullfile(reportDir,'prepare_caiso_oasis_inputs.m'));
        zipNames{end+1}='prepare_caiso_oasis_inputs.m';
    end
    zipPath=fullfile(reportDir,'step3_validation.zip');
    zip(zipPath,zipNames,reportDir);
    fprintf('\n%s\nMAT file: %s\nUpload ZIP: %s\nElapsed: %.2f seconds\n', ...
        summary,dataFile,zipPath,toc(started));
end

function T=readTypedCSV(path,names,numericNames)
    io=detectImportOptions(path,'Delimiter',',','TextType','string','VariableNamingRule','preserve');
    if ~isequal(io.VariableNames,names)
        error('VAREN:CSVSchema','Unexpected CSV columns in %s.',path);
    end
    io.DataLines=[2 Inf];
    io=setvartype(io,numericNames,'double');
    textNames=setdiff(names,numericNames,'stable');
    io=setvartype(io,textNames,'string');
    T=readtable(path,io);
end

function t=parseUTC(s,label)
    if ~all(matchesPattern(s,'^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$'))
        error('VAREN:UTCText','Invalid explicit UTC text in %s.',label);
    end
    t=datetime(s,'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z''','TimeZone','UTC','Locale','en_US');
    if any(isnat(t)) || ~isequal(utcText(t),s)
        error('VAREN:UTCParse','UTC parse/roundtrip failed in %s.',label);
    end
end

function q=matchesPattern(s,pattern)
    q=false(size(s)); present=nonblank(s);
    q(present)=~cellfun('isempty',regexp(cellstr(s(present)),pattern,'once'));
end

function q=nonblank(s)
    q=~ismissing(s) & strlength(strtrim(s))>0;
end

function s=utcText(t)
    t.TimeZone='UTC'; t.Format='yyyy-MM-dd HH:mm:ss';
    s=replace(string(t),' ','T')+"Z";
end

function verifyMetadata(m,n,blocks)
    required={'status','rows','query_blocks','interpolations','duplicate_rows_removed','missing_hours','invalid_loads'};
    if ~isstruct(m) || ~all(isfield(m,required)) || ...
            ~strcmp(m.status,'COMPLETE_WITH_NO_GAPS_OR_DUPLICATES') || ...
            ~isequal(m.rows,n) || ~isequal(m.query_blocks,blocks) || ...
            any([m.interpolations m.duplicate_rows_removed m.missing_hours m.invalid_loads]~=0)
        error('VAREN:Metadata','Acquisition metadata does not attest the required complete unfilled series.');
    end
end

function writeText(path,value)
    [fid,msg]=fopen(path,'w','n','UTF-8');
    if fid<0, error('VAREN:Write','Cannot write summary: %s',msg); end
    cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s',value);
end
