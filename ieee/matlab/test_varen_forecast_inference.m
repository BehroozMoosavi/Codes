function Report = test_varen_forecast_inference()
%TEST_VAREN_FORECAST_INFERENCE Deterministic synthetic checks; no real results.
% Independent reference values were calculated in Python/NumPy/SciPy using
% direct, unscaled Bartlett autocovariances and joint loss-covariance delta
% formulas. No MATLAB execution was available when this test was authored.
% The wrapper must execute it successfully before processing real forecasts.
    tick=tic;before=rng;
    checks=struct('Check',{},'Pass',{},'MaxAbsDifference',{});
    n=1008;t=(1:n)';
    y=25000+20*(mod(t,24)-12)+3*(mod(t,168)-84);
    ea=3*(mod(t,17)-8)+2*(mod(t,7)-3)+.1*mod(t,31);
    eb=ea+.05*(mod(t,13)-6)+.001;
    ec=.8*ea+.2*(mod(t,11)-5);
    ed=1.3*ea+mod(t,19)-9-1;
    pred=y-[ea eb eb ec ed];
    names={'VAREN','StaticEN','LASSO','FixedEN05','LSBoost24h'};
    hr=mod(t-1,24);
    H=[ones(n,1),.25+.5*mod(t,97)/96,double(hr>=10&hr<=14),double(hr>=16&hr<=21)];
    opts=struct('Instruments',H,'InstrumentNames',{{'Constant','Volatility','Midday','Evening'}});
    [T,D]=varen_forecast_inference(y,pred,names,opts);
    record('No RNG changes',isequal(before,rng),0);
    record('No fitting or random-number use',D.ForecastsFitted==0&&~D.RandomNumbersUsed,0);
    record('Expected row geometry',height(T.DM)==24&&height(T.Effects)==12&&height(T.GW)==24,0);
    record('All regular fixture tests defined',all(strcmp(T.DM.Status,'OK'))&&all(strcmp(T.GW.Status,'OK')),0);
    record('Planned Holm families preserved',all(T.DM.HolmFamilySize==3)&&all(T.GW.HolmFamilySize==3),0);
    textcheck('DM text columns are string columns',T.DM,{'ModelA','ModelB','Loss','Status','AnalysisRole'});
    textcheck('Effect text columns are string columns',T.Effects,{'ModelA','ModelB','RMSE_Status', ...
        'RelativeRMSE_Status','MAE_Status','AnalysisRole'});
    textcheck('GW text columns are string columns',T.GW,{'ModelA','ModelB','Loss','Status','AnalysisRole'});

    % Columns: lag, MSE gain/LRV/SE/Z/p/Holm, GW Wald/p, RMSE gain/SE/CI,
    % relative RMSE gain/SE/CI, MAE gain/SE/CI. Each row is independent Python output.
    expected=[ ...
    168 .04448195238136373 3.5160231068453109 .059060292596234978 .75316173398367825 .45135273098538431 .45135273098538431 2.0003808196556445 .73568883448225519 .0014533141428696439 .0019287505562105617 -.0023269674824646543 .005233595768203942 .0094960836937496573 .012596444276033006 -.015192493420540749 .034184660808040064 .0012500000000092394 .0024466511789316241 -.0035453481934292069 .0060453481934476857; ...
    24 .04448195238136373 76.761484249625184 .27595700049221422 .16119160703306287 .87194249129637669 .87194249129637669 .70960485459954881 .95013918583854373 .0014533141428696439 .0090159348955896421 -.016217593539443947 .019124221825183235 .0094960836937496573 .058907602152047309 -.10596069493987725 .12495286232737657 .0012500000000092394 .010101194852880871 -.018547978112458636 .021047978112477116; ...
    48 .04448195238136373 12.5197039288709 .11144658451149002 .39913248644042287 .68979558625367532 .68979558625367532 .80167826339102899 .93822292823662112 .0014533141428696439 .0036411112334404565 -.0056831327383778637 .0085897610241171506 .0094960836937496573 .023789885764779074 -.037131235601539449 .056123402989038763 .0012500000000092394 .0041807191890680002 -.0069440590400395427 .0094440590400580206; ...
    720 .04448195238136373 2.428575587333679 .049084632808144707 .90622970645881562 .36481428261259607 .36481428261259607 8.0791683318873186 .088720828818886308 .0014533141428696439 .0016024291794541476 -.0016873893366365564 .0045940176223758442 .0094960836937496573 .01046172610225953 -.011008522682801621 .030000690070300935 .0012500000000092394 .0017208906426683821 -.0021228836809527771 .0046228836809712563];
    dmFields={'MeanLossGain_BminusA','LongRunVariance','SE_MeanLossGain','DM_Z','P_TwoSided','P_Holm'};
    effectFields={'RMSE_Gain_MW','RMSE_SE_MW','RMSE_CI95_Lower_MW','RMSE_CI95_Upper_MW', ...
        'RMSE_Gain_pct','RMSE_SE_pct','RMSE_CI95_Lower_pct','RMSE_CI95_Upper_pct', ...
        'MAE_Gain_MW','MAE_SE_MW','MAE_CI95_Lower_MW','MAE_CI95_Upper_MW'};
    for j=1:size(expected,1)
        L=expected(j,1);
        dr=T.DM(T.DM.HACLag==L&strcmp(T.DM.ModelB,'StaticEN')&strcmp(T.DM.Loss,'MSE'),:);
        gr=T.GW(T.GW.HACLag==L&strcmp(T.GW.ModelB,'StaticEN')&strcmp(T.GW.Loss,'MSE'),:);
        er=T.Effects(T.Effects.HACLag==L&strcmp(T.Effects.ModelB,'StaticEN'),:);
        actual=[dr{1,dmFields},gr.WaldStatistic,gr.P_Joint,er{1,effectFields}];
        closecheck(sprintf('Independent numerical fixture lag%d',L),actual,expected(j,2:end),2e-8,2e-8);
    end
    closecheck('Benefit and harm signs',T.DM.MeanLossGain_BminusA(T.DM.HACLag==168&strcmp(T.DM.Loss,'MSE')), ...
        [.04448195238136373;-83.9031285714286;191.31848253968235],2e-8,2e-8);

    % Reversing a pair reverses DM/effect direction without changing two-sided p.
    reverse=opts;reverse.HACLags=168;reverse.PrimaryHACLag=168;
    reverse.PairNames=fliplr(D.Config.PairNames);
    [R,~]=varen_forecast_inference(y,pred,names,reverse);
    main=T.DM(T.DM.HACLag==168,:);mainEffects=T.Effects(T.Effects.HACLag==168,:);
    closecheck('Pair reversal DM sign',R.DM.DM_Z,-main.DM_Z,2e-8,2e-8);
    closecheck('Pair reversal p invariant',R.DM.P_TwoSided,main.P_TwoSided,2e-8,2e-8);
    closecheck('Pair reversal RMSE gain',R.Effects.RMSE_Gain_MW,-mainEffects.RMSE_Gain_MW,2e-8,2e-8);

    % Unit changes must leave all p-values and percent effects invariant.
    small=opts;small.HACLags=168;small.PrimaryHACLag=168;
    [S,~]=varen_forecast_inference(2*y,2*pred,names,small);
    closecheck('Unit scale p invariant',S.DM.P_TwoSided,main.P_TwoSided,2e-8,2e-8);
    closecheck('Unit scale RMSE CI',S.Effects.RMSE_CI95_Lower_MW,2*mainEffects.RMSE_CI95_Lower_MW,2e-8,2e-8);
    closecheck('Unit scale percent CI invariant',S.Effects.RMSE_CI95_Upper_pct,mainEffects.RMSE_CI95_Upper_pct,2e-8,2e-8);
    Q=[1 3 0 0;0 7 0 0;0 0 2 0;0 0 0 .5];small.Instruments=H*Q;
    [Qout,~]=varen_forecast_inference(y,pred,names,small);
    closecheck('GW nonsingular instrument transformation',Qout.GW.WaldStatistic,T.GW.WaldStatistic(T.GW.HACLag==168),1e-6,2e-8);
    one=small;one.Instruments=ones(n,1);one.InstrumentNames={'Constant'};
    [One,~]=varen_forecast_inference(y,pred,names,one);
    closecheck('Constant-only GW equals squared DM',One.GW.WaldStatistic,One.DM.DM_Z.^2,1e-6,2e-8);
    closecheck('Constant-only GW p equals two-sided DM p',One.GW.P_Joint,One.DM.P_TwoSided,2e-8,2e-8);

    % Exact losses give undefined tests. Their displayed zero differences are descriptive.
    dup=struct('HACLags',0,'PrimaryHACLag',0,'PairNames',{{'StaticEN','LASSO'}});
    [U,~]=varen_forecast_inference(y,pred,names,dup);
    record('Duplicate losses not assigned p=1',all(isnan(U.DM.P_TwoSided))&&all(isnan(U.DM.P_Holm))&& ...
        all(strcmp(U.DM.Status,'IDENTICAL_LOSSES'))&&all(U.DM.IdenticalPredictions),0);
    textcheck('Empty GW keeps string column types',U.GW,{'ModelA','ModelB','Loss','Status','AnalysisRole'});
    textcheck('Single-row effects keep string column types',U.Effects,{'ModelA','ModelB','RMSE_Status', ...
        'RelativeRMSE_Status','MAE_Status','AnalysisRole'});
    constant=struct('HACLags',0,'PrimaryHACLag',0,'PairNames',{{'A','B'}});
    [C,~]=varen_forecast_inference(100*ones(40,1),[99*ones(40,1),98*ones(40,1)],{'A','B'},constant);
    record('Constant nonzero differential undefined',all(isnan(C.DM.P_TwoSided))&&all(strcmp(C.DM.Status,'ZERO_LONG_RUN_VARIANCE')),0);
    [Z,~]=varen_forecast_inference(100*ones(40,1),[100*ones(40,1),100-mod((1:40)',3)],{'A','B'},constant);
    record('Zero MSE blocks nonsmooth delta intervals',strcmp(Z.Effects.RMSE_Status,'NONSMOOTH_ZERO_MSE')&& ...
        isnan(Z.Effects.RMSE_CI95_Lower_MW)&&isnan(Z.Effects.RMSE_CI95_Lower_pct),0);
    singular=opts;singular.HACLags=0;singular.PrimaryHACLag=0;
    singular.Instruments=[ones(n,1),ones(n,1)];singular.InstrumentNames={'Constant','DuplicateConstant'};
    [U,~]=varen_forecast_inference(y,pred,names,singular);
    record('Singular GW is reported, not inverted',all(isnan(U.GW.P_Joint))&&all(U.GW.HolmFamilySize==3),0);

    % Construct three known two-sided normal p-values at lag0; check full Holm order.
    nh=70;u=repmat((-3:3)',10,1);se=sqrt(mean(u.^2)/nh);
    targetP=[.04 .001 .02];z=sqrt(2)*erfcinv(targetP);
    yh=1000*ones(nh,1);ph=yh-[10*ones(nh,1),sqrt(100+u+se*z)];
    hnames={'VAREN','StaticEN','FixedEN05','LSBoost24h'};
    hopts=struct('HACLags',0,'PrimaryHACLag',0);
    [Hout,~]=varen_forecast_inference(yh,ph,hnames,hopts);
    rows=strcmp(Hout.DM.Loss,'MSE');
    closecheck('Known-p fixture',Hout.DM.P_TwoSided(rows),targetP',2e-8,2e-8);
    closecheck('Holm ordered cumulative adjustment',Hout.DM.P_Holm(rows),[.04;.003;.04],2e-8,2e-8);
    ph(:,3)=ph(:,1);[Hout,~]=varen_forecast_inference(yh,ph,hnames,hopts);
    rows=strcmp(Hout.DM.Loss,'MSE');
    pp=Hout.DM.P_Holm(rows);
    record('Undefined comparison retains planned family',all(Hout.DM.HolmFamilySize==3)&&isnan(pp(2)),0);
    closecheck('Undefined comparison conservatively counts',pp([1 3]),[.08;.06],2e-8,2e-8);

    badNames=names;badNames{3}=names{2};
    expecterror('Duplicate column labels rejected',@()varen_forecast_inference(y,pred,badNames,opts));
    bad=opts;bad.HACLags=n;
    expecterror('Excess HAC lag rejected',@()varen_forecast_inference(y,pred,names,bad));
    bad=opts;bad.Instruments(1,1)=0;
    expecterror('Missing constant instrument rejected',@()varen_forecast_inference(y,pred,names,bad));
    badPred=pred;badPred(1,1)=NaN;
    expecterror('Nonfinite prediction rejected',@()varen_forecast_inference(y,badPred,names,opts));
    record('Caller RNG preserved through all checks',isequal(before,rng),0);
    Report=struct('Status','PASS','Checks',struct2table(checks),'CheckCount',numel(checks), ...
        'FixtureN',n,'Seconds',toc(tick),'CorePath',which('varen_forecast_inference'), ...
        'Provenance','Deterministic synthetic numerical fixture and invariants; no empirical forecasts fitted.');
    fprintf('PASS: %d saved-forecast inference checks; %.3f seconds.\n',Report.CheckCount,Report.Seconds);

    function record(label,pass,difference)
        checks(end+1)=struct('Check',label,'Pass',logical(pass),'MaxAbsDifference',difference); %#ok<AGROW>
        if ~pass,error('VAREN:InferenceSelfTest','Failed: %s.',label);end
    end
    function closecheck(label,actual,wanted,absolute,relative)
        delta=abs(actual-wanted);bound=absolute+relative*abs(wanted);
        pass=isequal(size(actual),size(wanted))&&all(isfinite(actual),'all')&&all(delta<=bound,'all');
        record(label,pass,max(delta,[],'all'));
    end
    function textcheck(label,V,fields)
        pass=all(cellfun(@(f)isstring(V.(f))&&isequal(size(V.(f)),[height(V) 1]),fields));
        record(label,pass,0);
    end
    function expecterror(label,fn)
        caught=false;
        try,fn();catch ME,caught=startsWith(ME.identifier,'VAREN:Inference');end
        record(label,caught,0);
    end
end
