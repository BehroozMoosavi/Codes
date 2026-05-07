%% ============================================================
%  data_clean.m
%  FILE 01 of 30
%
%  Build the CAISO hourly demand dataset for VA-REN.
%
%  TARGET VARIABLE
%  ---------------
%      z_t = 100 * [log(y_t) - log(y_{t-1})]
%
%  PREDICTOR SET  (70 features)
%  ----------------------------
%      26 lagged log-changes : z_{t-1}...z_{t-24}, z_{t-168}, z_{t-8760}
%       4 log-level anchors  : log(y_{t-1}), log(y_{t-24}),
%                              log(y_{t-168}), log(y_{t-8760})
%      23 hour-of-day dummies: hour 0 is reference
%       6 weekday dummies    : Sunday is reference
%      11 month dummies      : January is reference
%
%  VOLATILITY INDEX  (for VA-REN adaptive alpha)
%  -----------------------------------------------
%      Wv = 48  backward rolling variance of z_t
%      Normalized to [0,1] using 2019 calibration
%      alpha_t = 1 - NormalizedVolatility_t
%
%  STANDARDIZATION
%  ---------------
%      Fixed calibration: first 720 usable obs of 2020.
%      Frozen after calibration — no rolling standardization.
%
%  INPUT FILES
%  -----------
%      historicalemshourlyload-2019.xlsx
%      historicalemshourlyload-2020.xlsx
%      historicalemshourlyload-2021.xlsx
%      historicalemshourlyload-2022.xlsx
%      historicalemshourlyloadfor2023.xlsx
%
%  OUTPUT
%  ------
%      caiso_final_logchange_dataset_2019_2023.mat
%      caiso_hourly_full_2019_2023.csv
%      caiso_hourly_model_ready_logchange_2019_2023.csv
%      caiso_hourly_evaluation_logchange_2020_2023.csv
%% ============================================================

clear; clc;

files = {
    "historicalemshourlyload-2019.xlsx",   2019;
    "historicalemshourlyload-2020.xlsx",   2020;
    "historicalemshourlyload-2021.xlsx",   2021;
    "historicalemshourlyload-2022.xlsx",   2022;
    "historicalemshourlyloadfor2023.xlsx", 2023
};

allData = table();

%% ------------------------------------------------------------
%  Read each yearly Excel file
%% ------------------------------------------------------------

for i = 1:size(files,1)

    fileName = files{i,1};
    fileYear = files{i,2};

    fprintf("Reading %s ...\n", fileName);

    if ~isfile(fileName)
        error("File not found: %s", fileName);
    end

    C = readcell(fileName);

    headerRow = [];
    for r = 1:min(20, size(C,1))
        rowText  = strtrim(string(C(r,:)));
        hasDate  = any(strcmpi(rowText,"Date"));
        hasHour  = any(strcmpi(rowText,"HE")) || any(strcmpi(rowText,"HR"));
        hasCAISO = any(strcmpi(rowText,"CAISO Total")) || any(strcmpi(rowText,"CAISO"));
        if hasDate && hasHour && hasCAISO
            headerRow = r; break;
        end
    end

    if isempty(headerRow)
        error("Could not detect header row in %s.", fileName);
    end

    header  = strtrim(string(C(headerRow,:)));
    dateIdx = find(strcmpi(header,"Date"),1);
    hourIdx = find(strcmpi(header,"HE"),1);
    if isempty(hourIdx), hourIdx = find(strcmpi(header,"HR"),1); end
    loadIdx = find(strcmpi(header,"CAISO Total"),1);
    if isempty(loadIdx), loadIdx = find(strcmpi(header,"CAISO"),1); end

    if isempty(dateIdx), error("No Date column in %s.",       fileName); end
    if isempty(hourIdx), error("No HE/HR column in %s.",      fileName); end
    if isempty(loadIdx), error("No CAISO load column in %s.", fileName); end

    dateRaw = C(headerRow+1:end, dateIdx);
    hourRaw = C(headerRow+1:end, hourIdx);
    loadRaw = C(headerRow+1:end, loadIdx);

    dateVal = convertDateColumn(dateRaw);
    hourVal = convertNumericColumn(hourRaw);
    loadVal = convertNumericColumn(loadRaw);

    valid   = ~isnat(dateVal) & ~isnan(hourVal) & ~isnan(loadVal) & loadVal > 0;
    dateVal = dateVal(valid);
    hourVal = hourVal(valid);
    loadVal = loadVal(valid);

    timestamp = dateVal + hours(hourVal - 1);

    TY               = table();
    TY.Timestamp     = timestamp(:);
    TY.Date          = dateshift(timestamp(:),"start","day");
    TY.Year          = year(timestamp(:));
    TY.Month         = month(timestamp(:));
    TY.Day           = day(timestamp(:));
    TY.Hour          = hour(timestamp(:));
    TY.Weekday       = weekday(timestamp(:));
    TY.HE            = hourVal(:);
    TY.CAISO_Load_MW = loadVal(:);

    TY      = sortrows(TY,"Timestamp");
    allData = [allData; TY];

    fprintf("  Loaded %d rows for %d\n", height(TY), fileYear);

end

%% Sort and deduplicate
allData = sortrows(allData,"Timestamp");
[~,ia]  = unique(allData.Timestamp,"stable");
if numel(ia) < height(allData)
    warning("Duplicate timestamps found. Keeping first occurrence.");
    allData = allData(ia,:);
end

fprintf("\nMerged rows : %d\n", height(allData));
fprintf("Start       : %s\n",  string(allData.Timestamp(1)));
fprintf("End         : %s\n",  string(allData.Timestamp(end)));

%% Fill missing hourly timestamps
expectedTime = (allData.Timestamp(1):hours(1):allData.Timestamp(end))';
missingTimes = setdiff(expectedTime, allData.Timestamp);

fprintf("Expected rows : %d\n", numel(expectedTime));
fprintf("Missing rows  : %d\n", numel(missingTimes));

if ~isempty(missingTimes)
    warning("Missing timestamps detected. Interpolating.");
    fullTable           = table();
    fullTable.Timestamp = expectedTime;
    allData = outerjoin(fullTable, allData, "Keys","Timestamp", ...
                        "MergeKeys",true,"Type","left");
    allData = sortrows(allData,"Timestamp");
    allData.CAISO_Load_MW = fillmissing(allData.CAISO_Load_MW,"linear");
    allData.Date    = dateshift(allData.Timestamp,"start","day");
    allData.Year    = year(allData.Timestamp);
    allData.Month   = month(allData.Timestamp);
    allData.Day     = day(allData.Timestamp);
    allData.Hour    = hour(allData.Timestamp);
    allData.Weekday = weekday(allData.Timestamp);
    allData.HE      = allData.Hour + 1;
end

%% ------------------------------------------------------------
%  Transformed response
%% ------------------------------------------------------------

loadSeries             = allData.CAISO_Load_MW(:);
allData.Log_Load       = log(loadSeries);
allData.Delta_Load_MW  = [NaN; diff(loadSeries)];
allData.Delta_Log_Load = [NaN; 100*diff(allData.Log_Load)];

%% ------------------------------------------------------------
%  Lagged log-change predictors
%% ------------------------------------------------------------

zSeries = allData.Delta_Log_Load(:);

for L = 1:24
    allData.(sprintf("LogDeltaLag_%d",L)) = localLag(zSeries,L);
end
allData.LogDeltaLag_168  = localLag(zSeries,168);
allData.LogDeltaLag_8760 = localLag(zSeries,8760);

%% ------------------------------------------------------------
%  Log-level anchors + raw level lags
%% ------------------------------------------------------------

logLoadSeries = allData.Log_Load(:);

allData.LogLoadLag_1    = localLag(logLoadSeries,1);
allData.LogLoadLag_24   = localLag(logLoadSeries,24);
allData.LogLoadLag_168  = localLag(logLoadSeries,168);
allData.LogLoadLag_8760 = localLag(logLoadSeries,8760);

allData.LoadLag_1    = localLag(loadSeries,1);
allData.LoadLag_24   = localLag(loadSeries,24);
allData.LoadLag_168  = localLag(loadSeries,168);
allData.LoadLag_8760 = localLag(loadSeries,8760);

%% ------------------------------------------------------------
%  Calendar dummies
%% ------------------------------------------------------------

for h = 1:23
    allData.(sprintf("HourDummy_%02d",h)) = double(allData.Hour == h);
end
for d = 2:7
    allData.(sprintf("WeekdayDummy_%d",d)) = double(allData.Weekday == d);
end
for m = 2:12
    allData.(sprintf("MonthDummy_%02d",m)) = double(allData.Month == m);
end

%% ------------------------------------------------------------
%  Volatility index  (VA-REN adaptive alpha)
%% ------------------------------------------------------------

Wv     = 48;
rawVol = movvar(allData.Delta_Log_Load(:),[Wv-1,0],0,"omitnan");
allData.RawVolatility = rawVol(:);

cal2019 = rawVol(allData.Year == 2019);
cal2019 = cal2019(~isnan(cal2019));

if isempty(cal2019), error("No 2019 volatility values found."); end

Vmin = min(cal2019);
Vmax = max(cal2019);

if Vmax == Vmin, error("Vmax equals Vmin."); end

normVol = (rawVol - Vmin)./(Vmax - Vmin);
normVol = max(0, min(1, normVol));

allData.NormalizedVolatility = normVol(:);
allData.Alpha_VAREN          = 1 - allData.NormalizedVolatility;

fprintf("\nVolatility:\n");
fprintf("  Wv          = %d hrs\n", Wv);
fprintf("  Vmin (2019) = %.8f\n",   Vmin);
fprintf("  Vmax (2019) = %.8f\n",   Vmax);

%% ------------------------------------------------------------
%  Predictor and response names
%% ------------------------------------------------------------

responseName   = "Delta_Log_Load";
predictorNames = strings(0);

for L = 1:24
    predictorNames(end+1) = sprintf("LogDeltaLag_%d",L);
end
predictorNames(end+1) = "LogDeltaLag_168";
predictorNames(end+1) = "LogDeltaLag_8760";
predictorNames(end+1) = "LogLoadLag_1";
predictorNames(end+1) = "LogLoadLag_24";
predictorNames(end+1) = "LogLoadLag_168";
predictorNames(end+1) = "LogLoadLag_8760";
for h = 1:23
    predictorNames(end+1) = sprintf("HourDummy_%02d",h);
end
for d = 2:7
    predictorNames(end+1) = sprintf("WeekdayDummy_%d",d);
end
for m = 2:12
    predictorNames(end+1) = sprintf("MonthDummy_%02d",m);
end

fprintf("\nNumber of predictors: %d\n", numel(predictorNames));

%% ------------------------------------------------------------
%  Model-ready dataset
%% ------------------------------------------------------------

usableIdx = true(height(allData),1);
usableIdx = usableIdx & ~isnan(allData.(responseName));
for j = 1:numel(predictorNames)
    usableIdx = usableIdx & ~isnan(allData.(predictorNames(j)));
end
usableIdx = usableIdx & ...
    ~isnan(allData.RawVolatility)        & ...
    ~isnan(allData.NormalizedVolatility) & ...
    ~isnan(allData.Alpha_VAREN)          & ...
    ~isnan(allData.LoadLag_1)            & ...
    ~isnan(allData.LoadLag_168);

modelData = allData(usableIdx,:);

fprintf("\nModel-ready rows : %d\n", height(modelData));
fprintf("First row        : %s\n",  string(modelData.Timestamp(1)));
fprintf("Last  row        : %s\n",  string(modelData.Timestamp(end)));

%% ------------------------------------------------------------
%  Standardization  (first 720 obs of 2020, frozen)
%% ------------------------------------------------------------

W_standardize    = 720;
evalCandIdx      = find(modelData.Year >= 2020 & modelData.Year <= 2023);

if numel(evalCandIdx) < W_standardize
    error("Not enough 2020-2023 obs for standardization.");
end

stdIdx = false(height(modelData),1);
stdIdx(evalCandIdx(1:W_standardize)) = true;

fprintf("\nStandardization sample:\n");
fprintf("  Start : %s\n", string(modelData.Timestamp(evalCandIdx(1))));
fprintf("  End   : %s\n", string(modelData.Timestamp(evalCandIdx(W_standardize))));
fprintf("  Rows  : %d\n", W_standardize);

X_raw   = modelData{:,predictorNames};
X_calib = X_raw(stdIdx,:);
X_mean  = mean(X_calib,1,"omitnan");
X_std   = std(X_calib,0,1,"omitnan");
X_std(X_std == 0 | isnan(X_std)) = 1;
X_standardized = (X_raw - X_mean)./X_std;

standardizedNames = strings(size(predictorNames));
for j = 1:numel(predictorNames)
    standardizedNames(j)             = "Z_" + predictorNames(j);
    modelData.(standardizedNames(j)) = X_standardized(:,j);
end

z_raw   = modelData.(responseName);
z_calib = z_raw(stdIdx);
z_mean  = mean(z_calib,"omitnan");
z_std   = std(z_calib,0,"omitnan");

if z_std == 0 || isnan(z_std)
    error("Response std is zero or NaN.");
end

fprintf("\nz_mean = %.8f\n", z_mean);
fprintf("z_std  = %.8f\n",  z_std);

modelData.Response_Standardized = (z_raw - z_mean)./z_std;

%% ------------------------------------------------------------
%  Evaluation sample 2020-2023
%% ------------------------------------------------------------

evaluationData    = modelData(modelData.Year >= 2020 & modelData.Year <= 2023,:);
X                 = evaluationData{:,standardizedNames};
y_eval_z_std      = evaluationData.Response_Standardized(:);
y_eval_z_raw      = evaluationData.Delta_Log_Load(:);
y_eval_level      = evaluationData.CAISO_Load_MW(:);
y_lag1_eval_level = evaluationData.LoadLag_1(:);
timestamps_eval   = evaluationData.Timestamp;

fprintf("\nEval rows  : %d\n", height(evaluationData));
fprintf("Eval start : %s\n",  string(evaluationData.Timestamp(1)));
fprintf("Eval end   : %s\n",  string(evaluationData.Timestamp(end)));

%% Year summaries
fprintf("\nFull dataset by year:\n");
for yr = unique(allData.Year)'
    idx = allData.Year == yr;
    fprintf("  %d: N=%d  meanLoad=%.0f  minLoad=%.0f  maxLoad=%.0f\n", ...
        yr, sum(idx), ...
        mean(allData.CAISO_Load_MW(idx),"omitnan"), ...
        min(allData.CAISO_Load_MW(idx),[],"omitnan"), ...
        max(allData.CAISO_Load_MW(idx),[],"omitnan"));
end

fprintf("\nEvaluation by year:\n");
for yr = unique(evaluationData.Year)'
    idx = evaluationData.Year == yr;
    fprintf("  %d: N=%d  meanLoad=%.0f  meanZ=%.4f  stdZ=%.4f\n", ...
        yr, sum(idx), ...
        mean(evaluationData.CAISO_Load_MW(idx),"omitnan"), ...
        mean(evaluationData.Delta_Log_Load(idx),"omitnan"), ...
        std(evaluationData.Delta_Log_Load(idx),"omitnan"));
end

%% Yearly volatility summary
fprintf("\nVolatility by year:\n");
for yr = unique(allData.Year)'
    idx     = allData.Year == yr;
    vol_yr  = allData.NormalizedVolatility(idx);
    fprintf("  %d: mean=%.4f  std=%.4f  max=%.4f\n", ...
        yr, mean(vol_yr,"omitnan"), std(vol_yr,"omitnan"), ...
        max(vol_yr,[],"omitnan"));
end

%% Evaluation sanity check
fprintf("\n=== EVALUATION BOUNDARIES ===\n");
fprintf("  First forecast origin (t=W+1) : obs %d\n", W_standardize+1);
fprintf("  Timestamp                     : %s\n", ...
        string(evaluationData.Timestamp(W_standardize+1)));
fprintf("  Calibration ends at           : %s\n", ...
        string(evaluationData.Timestamp(W_standardize)));

%% ------------------------------------------------------------
%  Save
%% ------------------------------------------------------------

writetable(allData,        "caiso_hourly_full_2019_2023.csv");
writetable(modelData,      "caiso_hourly_model_ready_logchange_2019_2023.csv");
writetable(evaluationData, "caiso_hourly_evaluation_logchange_2020_2023.csv");

save("caiso_final_logchange_dataset_2019_2023.mat", ...
    "allData","modelData","evaluationData", ...
    "X","y_eval_z_std","y_eval_z_raw","y_eval_level", ...
    "y_lag1_eval_level","timestamps_eval", ...
    "predictorNames","standardizedNames","responseName", ...
    "X_mean","X_std","z_mean","z_std","Vmin","Vmax","Wv");

fprintf("\nSaved: caiso_final_logchange_dataset_2019_2023.mat\n");

%% ============================================================
%  Local helpers
%% ============================================================

function x = convertNumericColumn(rawCol)
    n = numel(rawCol); x = NaN(n,1);
    for k = 1:n
        v = rawCol{k};
        if isnumeric(v), x(k) = double(v);
        elseif islogical(v), x(k) = double(v);
        elseif ischar(v)||isstring(v)
            s = strtrim(string(v));
            if s==""||strcmpi(s,"NaN")||strcmpi(s,"NA"), x(k)=NaN;
            else, x(k)=str2double(s); end
        else, x(k)=NaN;
        end
    end
end

function dt = convertDateColumn(rawCol)
    n  = numel(rawCol); dt = NaT(n,1);
    formats = ["MM/dd/yyyy";"M/d/yyyy";"yyyy-MM-dd"; ...
               "dd-MMM-yyyy";"MM/dd/yy";"M/d/yy"];
    for k = 1:n
        v = rawCol{k};
        if isdatetime(v), dt(k)=v;
        elseif isnumeric(v)
            try, dt(k)=datetime(v,"ConvertFrom","excel"); catch, dt(k)=NaT; end
        elseif ischar(v)||isstring(v)
            s = strtrim(string(v));
            if s==""||strcmpi(s,"NaN")||strcmpi(s,"NA"), dt(k)=NaT;
            else
                parsed = NaT;
                for f=1:numel(formats)
                    try, parsed=datetime(s,"InputFormat",formats(f));
                    catch, end
                    if ~isnat(parsed), break; end
                end
                if isnat(parsed), try, parsed=datetime(s); catch, end; end
                dt(k)=parsed;
            end
        else, dt(k)=NaT;
        end
    end
    dt = dateshift(dt,"start","day");
end

function ylag = localLag(y,L)
    y=y(:); n=numel(y);
    if L>=n, ylag=NaN(n,1);
    else, ylag=[NaN(L,1); y(1:end-L)]; end
end
