function [summary, reportDir] = audit_caiso_source_files(dataDir)
%AUDIT_CAISO_SOURCE_FILES  Step 2: inspect the original CAISO Excel files.
%
%   [summary, reportDir] = audit_caiso_source_files;
%       Opens a folder picker for the five ORIGINAL Excel files.
%
%   [summary, reportDir] = audit_caiso_source_files('/path/to/data');
%       Uses the specified folder instead.
%
%   Designed for MATLAB R2025a. No additional toolbox is needed.
%   This function reads the workbooks and writes a NEW report directory under
%   the current MATLAB folder. It does not edit any source workbook, discard
%   duplicate records, fill gaps, assign a time zone, or fit forecasting models.
%
%   IMPORTANT: NominalClockLabel = reported Date + hours(HE_or_HR - 1)
%   is used ONLY to compare labels with the original importer's convention.
%   It is NOT an assertion about hour beginning/ending or true UTC instants.
%   An absent nominal label is NOT necessarily a missing physical observation.
%   Likewise, a repeated label is NOT necessarily an erroneous duplicate.
%
%   All imported rows below the detected header remain in source row order in
%   all_source_records.csv, including invalid and blank-core rows if present.
%   RawDate/RawHour/RawLoad preserve the values read by MATLAB as text; numeric
%   raw values use 17 significant digits. Original Excel formatting is not
%   represented. DateParseMethod documents interpretation of each date cell.
%
%   Output: step2_data_audit.zip contains the summary, source-label diagnostics,
%   provider notes, and this code. If ZIP creation is unavailable, individual
%   reports remain available in reportDir.
%
%   Revision 1.1, 2026-10-27. This audit does not certify the forecasting data.

    if nargin < 1 || isempty(dataDir)
        dataDir = uigetdir(pwd, ...
            'Select the folder containing the five ORIGINAL CAISO XLSX files');
        if isequal(dataDir, 0)
            error('VAREN:AuditCancelled', 'No input folder was selected.');
        end
    end
    dataDir = char(string(dataDir));
    if ~isfolder(dataDir)
        error('VAREN:MissingDataFolder', 'Folder does not exist: %s', dataDir);
    end

    files = ["historicalemshourlyload-2019.xlsx"; ...
             "historicalemshourlyload-2020.xlsx"; ...
             "historicalemshourlyload-2021.xlsx"; ...
             "historicalemshourlyload-2022.xlsx"; ...
             "historicalemshourlyloadfor2023.xlsx"];
    years = (2019:2023)';
    paths = fullfile(string(dataDir), files);
    missingFiles = files(~isfile(paths));
    if ~isempty(missingFiles)
        error('VAREN:MissingSourceFiles', ...
            'Selected folder is missing these original files:\n%s', ...
            char(strjoin(missingFiles, newline)));
    end

    % tempname supplies a unique suffix; previous reports are never reused.
    [~, uniqueSuffix] = fileparts(tempname);
    reportDir = fullfile(pwd, ['step2_data_audit_' uniqueSuffix]);
    if isfolder(reportDir) || isfile(reportDir)
        error('VAREN:AuditOutputExists', 'Choose a new output directory.');
    end
    [made, message] = mkdir(reportDir);
    if ~made
        error('VAREN:AuditOutputFolder', '%s', message);
    end

    fprintf('\n============================================\n');
    fprintf('VA-REN REVISION: STEP 2 SOURCE-DATA AUDIT\n');
    fprintf('============================================\n');
    fprintf('MATLAB: %s\n', version);
    fprintf('Source folder: %s\n', dataDir);
    fprintf('Reading original labels without applying corrections.\n\n');
    auditTimer = tic;

    recordsByFile = cell(numel(files), 1);
    summaryByFile = cell(numel(files), 1);
    dailyByFile = cell(numel(files), 1);
    absentByFile = cell(numel(files), 1);
    orderByFile = cell(numel(files), 1);
    notesByFile = cell(numel(files), 1);
    inventoryByFile = cell(numel(files), 1);

    for k = 1:numel(files)
        fprintf('Reading %s ...\n', char(files(k)));
        sheets = sheetnames(paths(k));
        if isempty(sheets)
            error('VAREN:NoWorksheet', 'No worksheet found in %s.', files(k));
        end
        % The original workbooks each contain their data on the first sheet.
        % Do not select a sheet by year: the 2020 sheet is named "2019 ...".
        selectedSheet = sheets(1);
        % A1 anchors the import at worksheet row 1. MATLAB does not allow
        % NumHeaderLines with a starting-cell Range; detectHeader handles it.
        C = readcell(paths(k), 'Sheet', selectedSheet, ...
            'Range', 'A1', 'UseExcel', false);
        [headerRow, dateCol, hourCol, loadCol] = detectHeader(C, files(k));
        if headerRow == size(C, 1)
            error('VAREN:EmptySourceData', 'No data below header in %s.', files(k));
        end
        [T, S, D, A, O] = auditOneSheet(C, headerRow, dateCol, ...
            hourCol, loadCol, files(k), selectedSheet, years(k));
        N = collectNotes(C, headerRow, [dateCol hourCol loadCol], ...
            files(k), selectedSheet, years(k));
        S.NoteCells = height(N);
        info = dir(paths(k));
        inventoryByFile{k} = table(files(k), years(k), selectedSheet, ...
            strjoin(sheets, ' | '), info.bytes, headerRow, size(C,1), ...
            size(C,2), rawText(C{headerRow,dateCol}), ...
            rawText(C{headerRow,hourCol}), rawText(C{headerRow,loadCol}), ...
            'VariableNames', {'File','WorkbookYear','SelectedSheet', ...
            'AllSheets','FileBytes','HeaderRow','ImportedWorksheetRows', ...
            'ImportedWorksheetColumns','DateHeader','HourHeader','LoadHeader'});
        recordsByFile{k} = T;
        summaryByFile{k} = S;
        dailyByFile{k} = D;
        absentByFile{k} = A;
        orderByFile{k} = O;
        notesByFile{k} = N;
        fprintf(['  %d records; %d repeated-label group(s); ' ...
            '%d absent nominal hour(s); %d source-order step flag(s).\n'], ...
            S.SourceRecords, S.DuplicateLabelGroups, ...
            S.AbsentNominalHours, S.SourceOrderStepFlags);
    end

    summary = vertcat(summaryByFile{:});
    records = vertcat(recordsByFile{:});
    daily = vertcat(dailyByFile{:});
    absent = vertcat(absentByFile{:});
    orderFlags = vertcat(orderByFile{:});
    notes = vertcat(notesByFile{:});
    inventory = vertcat(inventoryByFile{:});
    duplicateRows = records(records.LabelMultiplicity > 1, :);
    clockReview = records(strlength(records.ReviewReason) > 0, :);
    dailyExceptions = daily(daily.SourceRecords ~= 24 | ...
        daily.UniqueValidHourLabels ~= 24, :);

    reportNames = { ...
        'source_file_summary.csv', 'all_source_records.csv', ...
        'duplicate_source_labels.csv', 'absent_nominal_labels.csv', ...
        'daily_record_counts.csv', 'daily_label_exceptions.csv', ...
        'source_order_anomalies.csv', 'clock_review_rows.csv', ...
        'source_notes.csv', 'source_inventory.csv'};
    reportTables = {summary, records, duplicateRows, absent, daily, ...
        dailyExceptions, orderFlags, clockReview, notes, inventory};
    for j = 1:numel(reportNames)
        writetable(reportTables{j}, fullfile(reportDir, reportNames{j}), ...
            'Encoding', 'UTF-8');
    end
    % Keep the exact audit code with the output when available on disk.
    auditCode = [mfilename('fullpath') '.m'];
    if isfile(auditCode)
        copyfile(auditCode, fullfile(reportDir, 'audit_caiso_source_files.m'));
        reportNames{end+1} = 'audit_caiso_source_files.m';
    end

    summaryText = makeSummaryText(summary, dailyExceptions, absent, notes);
    textPath = fullfile(reportDir, 'step2_summary.txt');
    writeText(textPath, summaryText);
    reportNames{end+1} = 'step2_summary.txt';
    fprintf('\n%s\n', summaryText);

    zipPath = fullfile(reportDir, 'step2_data_audit.zip');
    try
        zip(zipPath, reportNames, reportDir);
        fprintf('Upload this ZIP: %s\n', zipPath);
    catch ME
        fprintf(['ZIP creation was unavailable: %s\n' ...
            'The individual audit reports have been saved successfully.\n' ...
            'Send step2_summary.txt and source_file_summary.csv instead.\n'], ...
            ME.message);
    end
    fprintf('Report folder: %s\n', reportDir);
    fprintf('Elapsed time: %.2f seconds\n', toc(auditTimer));
    fprintf('STEP2_STATUS: AUDIT_COMPLETED\n');
    fprintf('DATA_CLOCK_STATUS: NOT_YET_VERIFIED\n');
    fprintf('============================================\n');
end

function [headerRow, dateCol, hourCol, loadCol] = detectHeader(C, fileName)
    headerRow = [];
    dateCol = [];
    hourCol = [];
    loadCol = [];
    for r = 1:min(30, size(C,1))
        labels = strings(1, size(C,2));
        for c = 1:size(C,2)
            labels(c) = strtrim(rawText(C{r,c}));
        end
        d = find(strcmpi(labels, 'Date'));
        h = find(strcmpi(labels, 'HE') | strcmpi(labels, 'HR'));
        l = find(strcmpi(labels, 'CAISO Total') | strcmpi(labels, 'CAISO'));
        if numel(d) == 1 && numel(h) == 1 && numel(l) == 1
            headerRow = r;
            dateCol = d;
            hourCol = h;
            loadCol = l;
            break
        end
    end
    if isempty(headerRow)
        error('VAREN:SourceHeader', ...
            'Could not identify unique Date, HE/HR, and CAISO columns in %s.', ...
            fileName);
    end
end

function [T, S, D, A, O] = auditOneSheet(C, headerRow, dc, hc, lc, file, sheet, y)
    n = size(C,1) - headerRow;
    T = table;
    T.File = repmat(file, n, 1);
    T.WorkbookYear = repmat(y, n, 1);
    T.Sheet = repmat(sheet, n, 1);
    T.SourceRow = (headerRow+1:size(C,1))';
    T.RawDate = strings(n,1);
    T.RawHour = strings(n,1);
    T.RawLoad = strings(n,1);
    T.DateCellClass = strings(n,1);
    T.HourCellClass = strings(n,1);
    T.LoadCellClass = strings(n,1);
    T.DateParseMethod = strings(n,1);
    T.ParsedDate = NaT(n,1);
    T.HourLabel = nan(n,1);
    T.LoadMW = nan(n,1);
    T.IsDataCandidate = false(n,1);
    for i = 1:n
        r = T.SourceRow(i);
        T.RawDate(i) = rawText(C{r,dc});
        T.RawHour(i) = rawText(C{r,hc});
        T.RawLoad(i) = rawText(C{r,lc});
        T.DateCellClass(i) = string(class(C{r,dc}));
        T.HourCellClass(i) = string(class(C{r,hc}));
        T.LoadCellClass(i) = string(class(C{r,lc}));
        [T.ParsedDate(i), T.DateParseMethod(i)] = parseDate(C{r,dc});
        T.HourLabel(i) = parseNumber(C{r,hc});
        T.LoadMW(i) = parseNumber(C{r,lc});
        T.IsDataCandidate(i) = ~isBlank(C{r,dc}) || ...
            ~isBlank(C{r,hc}) || ~isBlank(C{r,lc});
    end
    T.ParsedDate.Format = 'yyyy-MM-dd HH:mm:ss';
    dayLabel = dateshift(T.ParsedDate, 'start', 'day');
    T.DateValid = ~isnat(T.ParsedDate);
    T.DateAtMidnight = T.DateValid & T.ParsedDate == dayLabel;
    T.DateInWorkbookYear = T.DateValid & year(T.ParsedDate) == y;
    T.HourValid = isfinite(T.HourLabel) & T.HourLabel == fix(T.HourLabel) & ...
        T.HourLabel >= 1 & T.HourLabel <= 24;
    T.LoadValid = isfinite(T.LoadMW) & T.LoadMW > 0;
    validLabel = T.DateAtMidnight & T.HourValid;
    T.NominalClockLabel = NaT(n,1);
    T.NominalClockLabel(validLabel) = T.ParsedDate(validLabel) + ...
        hours(T.HourLabel(validLabel)-1);
    T.NominalClockLabel.Format = 'yyyy-MM-dd HH:mm:ss';
    T.LabelMultiplicity = zeros(n,1);
    multiplicities = zeros(0,1);
    if any(validLabel)
        [~,~,g] = unique(T.NominalClockLabel(validLabel));
        multiplicities = accumarray(g, 1);
        T.LabelMultiplicity(validLabel) = multiplicities(g);
    end
    % Diagnostics only: grouping labels above does not remove any source row.
    T.NominalStepFromPriorRowHours = nan(n,1);
    if n > 1
        T.NominalStepFromPriorRowHours(2:end) = ...
            hours(diff(T.NominalClockLabel));
    end
    stepFlag = isfinite(T.NominalStepFromPriorRowHours) & ...
        abs(T.NominalStepFromPriorRowHours-1) > 1e-9;
    T.HE24AfterHE23DateChanged = false(n,1);
    if n > 1
        T.HE24AfterHE23DateChanged(2:end) = ...
            T.HourLabel(1:end-1) == 23 & T.HourLabel(2:end) == 24 & ...
            T.DateValid(1:end-1) & T.DateValid(2:end) & ...
            T.ParsedDate(1:end-1) ~= T.ParsedDate(2:end);
    end

    expected = (datetime(y,1,1):hours(1): ...
        datetime(y+1,1,1)-hours(1))';
    absent = setdiff(expected, T.NominalClockLabel(validLabel));
    absent.Format = 'yyyy-MM-dd HH:mm:ss';
    A = table(repmat(file,numel(absent),1), repmat(y,numel(absent),1), ...
        dateshift(absent,'start','day'), hour(absent)+1, absent, ...
        'VariableNames', {'File','WorkbookYear','ReportedDate', ...
        'AbsentHourLabel','NominalClockLabel'});
    A.ReportedDate.Format = 'yyyy-MM-dd';

    expectedDays = (datetime(y,1,1):caldays(1):datetime(y,12,31))';
    dayList = unique([expectedDays; dayLabel(T.DateValid)]);
    dayList.Format = 'yyyy-MM-dd';
    m = numel(dayList);
    D = table(repmat(file,m,1), repmat(y,m,1), dayList, ...
        zeros(m,1), zeros(m,1), strings(m,1), strings(m,1), ...
        'VariableNames', {'File','WorkbookYear','ReportedDate', ...
        'SourceRecords','UniqueValidHourLabels','AbsentHourLabels', ...
        'RepeatedHourLabels'});
    for j = 1:m
        onDay = T.DateValid & dayLabel == dayList(j) & T.IsDataCandidate;
        h = T.HourLabel(onDay & T.HourValid);
        D.SourceRecords(j) = sum(onDay);
        D.UniqueValidHourLabels(j) = numel(unique(h));
        D.AbsentHourLabels(j) = strjoin(string(setdiff((1:24)',h))', ' ');
        if ~isempty(h)
            hc = accumarray(h,1,[24 1]);
            D.RepeatedHourLabels(j) = strjoin(string(find(hc>1))', ' ');
        end
    end

    q = find(stepFlag);
    O = table(T.File(q), T.WorkbookYear(q), T.SourceRow(q-1), ...
        T.ParsedDate(q-1), T.HourLabel(q-1), T.LoadMW(q-1), ...
        T.SourceRow(q), T.ParsedDate(q), T.HourLabel(q), T.LoadMW(q), ...
        T.NominalStepFromPriorRowHours(q), T.HE24AfterHE23DateChanged(q), ...
        'VariableNames', {'File','WorkbookYear','PreviousSourceRow', ...
        'PreviousDate','PreviousHourLabel','PreviousLoadMW','SourceRow', ...
        'ReportedDate','HourLabel','LoadMW','NominalStepHours', ...
        'HE24AfterHE23DateChanged'});

    T.ReviewReason = strings(n,1);
    T.ReviewReason = addReason(T.ReviewReason, ~T.IsDataCandidate, 'blank_core_row');
    bad = T.IsDataCandidate & (~T.DateValid | ~T.DateAtMidnight | ...
        ~T.DateInWorkbookYear | ~T.HourValid | ~T.LoadValid);
    T.ReviewReason = addReason(T.ReviewReason, bad, 'invalid_or_unexpected_core_value');
    T.ReviewReason = addReason(T.ReviewReason, T.LabelMultiplicity > 1, 'repeated_label');
    T.ReviewReason = addReason(T.ReviewReason, stepFlag, 'nominal_step_not_one_hour');
    T.ReviewReason = addReason(T.ReviewReason, T.HE24AfterHE23DateChanged, ...
        'date_changes_between_H23_and_H24');
    previous = false(n,1);
    previous(q-1) = true;
    T.ReviewReason = addReason(T.ReviewReason, previous, 'precedes_flagged_step');
    exceptionalDays = D.ReportedDate(D.SourceRecords~=24 | D.UniqueValidHourLabels~=24);
    T.ReviewReason = addReason(T.ReviewReason, ismember(dayLabel,exceptionalDays), ...
        'day_has_nonstandard_label_count');
    % Include the transition dates and one neighboring day on either side.
    % These dates identify records for inspection, not a time-zone conversion.
    march1 = datetime(y,3,1);
    november1 = datetime(y,11,1);
    springDay = march1 + caldays(mod(1-weekday(march1),7)+7);
    fallDay = november1 + caldays(mod(1-weekday(november1),7));
    dstContext = abs(days(dayLabel-springDay)) <= 1 | ...
        abs(days(dayLabel-fallDay)) <= 1;
    T.ReviewReason = addReason(T.ReviewReason, dstContext, 'DST_calendar_context');

    S = table(y, file, sheet, n, sum(T.IsDataCandidate), ...
        sum(~T.IsDataCandidate), sum(T.IsDataCandidate & ~T.DateValid), ...
        sum(T.DateValid & ~T.DateAtMidnight), ...
        sum(T.DateValid & ~T.DateInWorkbookYear), ...
        sum(T.IsDataCandidate & ~T.HourValid), ...
        sum(T.IsDataCandidate & ~T.LoadValid), ...
        sum(multiplicities>1), sum(max(multiplicities-1,0)), ...
        numel(absent), sum(D.SourceRecords~=24), sum(stepFlag), ...
        sum(T.HE24AfterHE23DateChanged), ...
        'VariableNames', {'Year','File','Sheet','RowsBelowHeader', ...
        'SourceRecords','BlankCoreRows','InvalidDateRows','NonMidnightDates', ...
        'WrongYearRows','InvalidHourRows','InvalidLoadRows', ...
        'DuplicateLabelGroups','ExtraDuplicateRows','AbsentNominalHours', ...
        'DaysWithNon24Rows','SourceOrderStepFlags','HE24DateChangeRows'});
end

function N = collectNotes(C, headerRow, coreCols, file, sheet, y)
    nr = zeros(0,1);
    nc = zeros(0,1);
    nt = strings(0,1);
    % Preserve textual notes outside the Date/HE-or-HR/CAISO columns,
    % including notes beside the actual transition-day records.
    for c = 1:size(C,2)
        for r = 1:size(C,1)
            if r == headerRow || (r > headerRow && ismember(c,coreCols))
                continue
            end
            v = C{r,c};
            if (ischar(v) || (isstring(v) && isscalar(v))) && ~isBlank(v)
                nr(end+1,1) = r; %#ok<AGROW>
                nc(end+1,1) = c; %#ok<AGROW>
                nt(end+1,1) = string(v); %#ok<AGROW>
            end
        end
    end
    N = table(repmat(file,numel(nr),1), repmat(y,numel(nr),1), ...
        repmat(sheet,numel(nr),1), nr, nc, nt, ...
        'VariableNames', {'File','WorkbookYear','Sheet','SourceRow', ...
        'SourceColumnNumber','NoteText'});
end

function reason = addReason(reason, mask, value)
    empty = mask & strlength(reason)==0;
    used = mask & strlength(reason)>0;
    reason(empty) = string(value);
    reason(used) = reason(used) + "; " + string(value);
end

function tf = isBlank(v)
    if isempty(v)
        tf = true;
    elseif ischar(v)
        tf = isempty(strtrim(v));
    elseif isstring(v) && isscalar(v)
        tf = ismissing(v) || strlength(strtrim(v))==0;
    else
        % Missing numbers/dates still count as present core cells if another
        % field identifies a record. This avoids silently discarding bad rows.
        try
            tf = isscalar(v) && ismissing(v);
        catch
            tf = false;
        end
    end
end

function s = rawText(v)
    if isempty(v)
        s = "";
    elseif isnumeric(v) && isscalar(v)
        s = string(sprintf('%.17g',v));
    elseif isdatetime(v) && isscalar(v)
        if isnat(v)
            s = "NaT";
        else
            v.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
            s = string(v);
        end
    elseif ischar(v) || (isstring(v) && isscalar(v))
        s = string(v);
        if ismissing(s), s = ""; end
    elseif isBlank(v)
        s = "";
    else
        s = "<unsupported:" + string(class(v)) + ">";
    end
end

function x = parseNumber(v)
    x = NaN;
    if isnumeric(v) && isscalar(v) && isreal(v)
        x = double(v);
    elseif ischar(v) || (isstring(v) && isscalar(v))
        % Do not strip suffixes, punctuation, or DST markers from hour labels.
        x = str2double(strtrim(string(v)));
        if ~isscalar(x) || ~isreal(x), x = NaN; end
    end
end

function [d, method] = parseDate(v)
    d = NaT;
    method = "unparsed";
    if isBlank(v)
        method = "missing";
        return
    end
    if isdatetime(v) && isscalar(v)
        if isempty(v.TimeZone)
            d = v;
            method = "native_datetime";
        else
            method = "zoned_datetime_requires_review";
        end
        return
    end
    if isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v)
        % The original five XLSX workbooks use the Excel 1900 date system.
        % Unexpected years are reported, never silently adjusted.
        try
            d = datetime(v, 'ConvertFrom', 'excel');
            method = "excel_1900_serial";
        catch
            method = "invalid_excel_serial";
        end
        return
    end
    if ischar(v) || (isstring(v) && isscalar(v))
        s = strtrim(string(v));
        % Explicit formats only: slash dates in these US files mean M/d/yyyy.
        % Do not use locale-dependent automatic guessing or swap day/month.
        formats = ["yyyy-MM-dd HH:mm:ss"; "yyyy-MM-dd'T'HH:mm:ss"; ...
            "yyyy-MM-dd"; "M/d/yyyy HH:mm:ss"; "M/d/yyyy"; ...
            "d-MMM-yyyy HH:mm:ss"; "d-MMM-yyyy"];
        for j = 1:numel(formats)
            try
                candidate = datetime(s, 'InputFormat', formats(j), ...
                    'Locale', 'en_US');
                if ~isnat(candidate)
                    d = candidate;
                    method = "text:" + formats(j);
                    return
                end
            catch
                % Retain the raw value and flag the row if no explicit format fits.
            end
        end
    end
end

function report = makeSummaryText(S, D, A, N)
    lines = [ ...
        "VA-REN REVISION: STEP 2 SOURCE-DATA AUDIT"; ...
        "Audit revision: 1.1 (2026-10-27)"; ...
        "MATLAB: " + string(version); ...
        "Created: " + string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')); ...
        ""; ...
        "STEP2_STATUS: AUDIT_COMPLETED"; ...
        "DATA_CLOCK_STATUS: NOT_YET_VERIFIED"; ...
        ""; ...
        "This reports the supplied source labels; it does not certify physical hours."; ...
        "NominalClockLabel = reported Date + hours(HE_or_HR - 1), with no time zone."; ...
        "Absent/repeated labels can arise from DST or source date conventions."; ...
        "Source workbooks were not changed. No rows were repaired, deleted, or filled."; ...
        "No forecasting model or accuracy test was run."; ...
        ""; ...
        "PER-FILE COUNTS"; ...
        "Year | Records | Invalid date/hour/load | Repeat groups | Extra repeats | Absent labels | Non24 days | Order flags | H24 date changes"];
    for k = 1:height(S)
        lines(end+1,1) = string(sprintf( ...
            '%d | %d | %d/%d/%d | %d | %d | %d | %d | %d | %d', ...
            S.Year(k), S.SourceRecords(k), S.InvalidDateRows(k), ...
            S.InvalidHourRows(k), S.InvalidLoadRows(k), ...
            S.DuplicateLabelGroups(k), S.ExtraDuplicateRows(k), ...
            S.AbsentNominalHours(k), S.DaysWithNon24Rows(k), ...
            S.SourceOrderStepFlags(k), S.HE24DateChangeRows(k))); %#ok<AGROW>
    end
    lines(end+1,1) = "Total source records: " + string(sum(S.SourceRecords));
    lines(end+1,1) = "Blank-core rows retained: " + string(sum(S.BlankCoreRows));
    lines(end+1,1) = "Dates with non-midnight values: " + string(sum(S.NonMidnightDates));
    lines(end+1,1) = "Dates outside workbook year: " + string(sum(S.WrongYearRows));
    lines = [lines; ""; "REPORTED DATES WITH INCOMPLETE OR REPEATED HOUR LABELS"];
    if isempty(D), lines(end+1,1) = "None detected."; end
    for k = 1:height(D)
        lines(end+1,1) = string(sprintf( ...
            '%s: %s; records=%d; absent labels=[%s]; repeated labels=[%s]', ...
            D.File(k), string(D.ReportedDate(k)), D.SourceRecords(k), ...
            D.AbsentHourLabels(k), D.RepeatedHourLabels(k))); %#ok<AGROW>
    end
    lines = [lines; ""; "ABSENT NOMINAL LABELS (NOT CONFIRMED PHYSICAL GAPS)"];
    if isempty(A), lines(end+1,1) = "None detected."; end
    for k = 1:height(A)
        lines(end+1,1) = A.File(k) + ": " + string(A.ReportedDate(k)) + ...
            " hour label " + string(A.AbsentHourLabel(k)); %#ok<AGROW>
    end
    lines = [lines; ""; "SOURCE NOTES (COPIED AS READ)"];
    if isempty(N), lines(end+1,1) = "No external note cells detected."; end
    for k = 1:height(N)
        lines(end+1,1) = N.File(k) + ", worksheet row " + ...
            string(N.SourceRow(k)) + ", column " + ...
            string(N.SourceColumnNumber(k)) + ":"; %#ok<AGROW>
        lines(end+1,1) = N.NoteText(k); %#ok<AGROW>
    end
    lines = [lines; ""; "NEXT STEP"; ...
        "Return step2_data_audit.zip for review before constructing forecasting inputs."; ...
        "If ZIP creation is unavailable, return this text file and source_file_summary.csv."; ...
        "all_source_records.csv preserves imported source-row order and flags."; ...
        "clock_review_rows.csv includes affected rows and DST-date context."; ...
        "A completed audit is not approval to interpolate or relabel observations."];
    report = char(strjoin(lines, newline));
end

function writeText(path, content)
    [fid, message] = fopen(path, 'w', 'n', 'UTF-8');
    if fid < 0
        error('VAREN:AuditWriteFailed', 'Cannot write report: %s', message);
    end
    fileCleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '%s\n', content);
end
