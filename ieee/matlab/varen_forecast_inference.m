function [T,D] = varen_forecast_inference(y,pred,names,opts)
%VAREN_FORECAST_INFERENCE Paired, serial-dependence-aware saved-forecast analysis.
% [T,D] = varen_forecast_inference(y,pred,names,opts)
% T.DM: asymptotic two-sided loss-differential tests and pointwise 95% CIs.
% T.Effects: paired RMSE/relative-RMSE delta-method CIs and MAE-gain CIs.
% T.GW: optional secondary joint conditional-moment Wald tests.
%
% Rows MUST be the common, contiguous, physical-hour forecast sequence.
% The caller verifies timestamps; this numerical core receives no dates.
% Gains are comparator (ModelB) minus focal model (ModelA): positive favors A.
% This is revision-stage exploratory inference, not a preregistered analysis.
% No model fitting, selection, resampling, random draws, or file I/O occurs.
%
% Options (all optional):
% HACLags          [168 24 48 720], hours; primary lag168 is listed first.
% PrimaryHACLag    168; other bandwidths are sensitivity checks.
% PairNames        {'VAREN','StaticEN';'VAREN','FixedEN05';'VAREN','LSBoost24h'}.
% Instruments      [] or N-by-K origin-available instruments, first column1.
% InstrumentNames  corresponding K distinct names (required with instruments).
%
% Bartlett HAC: Omega=Gamma0+sum_{k=1}^L (1-k/(L+1))*(Gammak+Gammak').
% The centered sample autocovariances use divisor N at every lag. The mean
% has estimated covariance Omega/N. No small-sample DM adjustment is used.
% RMSE effects use the joint paired squared-loss influence function.
% Relative RMSE gain is100*(1-RMSE_A/RMSE_B), in percent.
% CIs are pointwise; Holm adjusts only p-values within the requested pairs,
% separately by test type, loss and bandwidth, preserving the planned family.
% Undefined tests count as p=1 internally but retain NaN reported p-values.
% No rejection is selected by comparing bandwidths. Inference requires
% suitable weak dependence/stability; HAC is not a remedy for arbitrary
% nonstationarity. A failure to reject does not establish equivalence.
%
% GW-style test: E[instruments_t * (loss_Bt-loss_At)]=0. This tests only the
% supplied moments, not all conditional predictive ability, and a rejection
% does not determine which model should be preferred. Instruments must be
% available at the forecast origin; the caller establishes their provenance.
% Positive column scaling conditions the Wald calculation without changing
% its null. Rank-deficient/ill-conditioned covariance gives undefined tests.
%
% References: Diebold & Mariano(1995), JBES13:253-263;
% Newey & West(1987), Econometrica55:703-708;
% Giacomini & White(2006), Econometrica74:1545-1578;
% Holm(1979), Scandinavian Journal of Statistics6:65-70.
% Source-reviewed; MATLAB execution must be checked by the supplied tests.

    if nargin<4, opts=struct(); end
    [y,pred,names,cfg]=validate_inputs(y,pred,names,opts);
    n=numel(y); pairN=size(cfg.PairNames,1); lagN=numel(cfg.HACLags);
    err=y-pred; sq=err.^2; ab=abs(err);
    if any(~isfinite(sq),'all')
        error('VAREN:InferenceOverflow','Squared forecast errors are not finite.');
    end
    z95=sqrt(2)*erfcinv(.05);
    dm=repmat(dm_record(),pairN*2*lagN,1);
    effects=repmat(effect_record(),pairN*lagN,1);
    gw=repmat(gw_record(),pairN*2*lagN*double(~isempty(cfg.Instruments)),1);
    di=0; ei=0; gi=0;
    for p=1:pairN
        a=find(strcmp(names,cfg.PairNames{p,1}));
        b=find(strcmp(names,cfg.PairNames{p,2}));
        mA=mean(sq(:,a)); mB=mean(sq(:,b));
        rA=sqrt(mA); rB=sqrt(mB);
        dM=sq(:,b)-sq(:,a); dA=ab(:,b)-ab(:,a);
        if any(~isfinite([dM;dA]))
            error('VAREN:InferenceOverflow','Paired loss differences are not finite.');
        end
        % Column1=MSE differential;2=MAE differential;3=RMSE influence;
        % column4=relative RMSE gain influence. Smoothness fails at zero MSE.
        smooth=(rA>0 && rB>0);
        influence=zeros(n,2);
        if smooth
            cA=sq(:,a)-mA; cB=sq(:,b)-mB;
            influence(:,1)=cB/(2*rB)-cA/(2*rA);
            influence(:,2)=50*(rA/rB)*(cB/mB-cA/mA);
        end
        [omega,omegaStatus]=scalar_hac([dM dA influence],cfg.HACLags);
        identicalPrediction=isequal(pred(:,a),pred(:,b));
        for j=1:lagN
            L=cfg.HACLags(j);
            for lossIndex=1:2
                if lossIndex==1
                    lossName='MSE'; delta=dM; meanA=mA; meanB=mB;
                else
                    lossName='MAE'; delta=dA; meanA=mean(ab(:,a)); meanB=mean(ab(:,b));
                end
                [se,lo,hi,z,pval,status]=scalar_test(delta,omega(j,lossIndex), ...
                    omegaStatus{j,lossIndex},z95);
                di=di+1;
                dm(di)=struct('PairIndex',p,'ModelA',names{a},'ModelB',names{b}, ...
                    'Loss',lossName,'HACLag',L,'N',n,'MeanLoss_A',meanA, ...
                    'MeanLoss_B',meanB,'MeanLossGain_BminusA',mean(delta), ...
                    'LongRunVariance',omega(j,lossIndex),'SE_MeanLossGain',se, ...
                    'CI95_Lower',lo,'CI95_Upper',hi,'DM_Z',z,'P_TwoSided',pval, ...
                    'P_Holm',NaN,'HolmFamilySize',0,'Status',status, ...
                    'IdenticalPredictions',identicalPrediction, ...
                    'AnalysisRole',analysis_role(lossName,L,cfg.PrimaryHACLag,false));
            end
            ei=ei+1; e=effect_record();
            e.PairIndex=p; e.ModelA=names{a}; e.ModelB=names{b};e.HACLag=L;e.N=n;
            e.RMSE_A_MW=rA;e.RMSE_B_MW=rB;e.RMSE_Gain_MW=rB-rA;
            if rB>0, e.RMSE_Gain_pct=100*(1-rA/rB);end
            e.MAE_Gain_MW=mean(dA);
            [e.MAE_SE_MW,e.MAE_CI95_Lower_MW,e.MAE_CI95_Upper_MW,~,~,e.MAE_Status]= ...
                scalar_test(dA,omega(j,2),omegaStatus{j,2},z95);
            if smooth
                [e.RMSE_SE_MW,e.RMSE_CI95_Lower_MW,e.RMSE_CI95_Upper_MW,e.RMSE_Status]= ...
                    effect_ci(e.RMSE_Gain_MW,omega(j,3),omegaStatus{j,3},n,z95,all(dM==0));
                [e.RMSE_SE_pct,e.RMSE_CI95_Lower_pct,e.RMSE_CI95_Upper_pct,e.RelativeRMSE_Status]= ...
                    effect_ci(e.RMSE_Gain_pct,omega(j,4),omegaStatus{j,4},n,z95,all(dM==0));
            end
            e.AnalysisRole=analysis_role('MSE',L,cfg.PrimaryHACLag,false);
            effects(ei)=e;
        end
        if ~isempty(cfg.Instruments)
            for lossIndex=1:2
                if lossIndex==1, delta=dM;lossName='MSE';else,delta=dA;lossName='MAE';end
                moments=cfg.Instruments.*delta;
                [wald,pval,rc,status]=joint_hac(moments,cfg.HACLags,all(delta==0));
                for j=1:lagN
                    gi=gi+1;
                    gw(gi)=struct('PairIndex',p,'ModelA',names{a},'ModelB',names{b}, ...
                        'Loss',lossName,'HACLag',cfg.HACLags(j),'N',n, ...
                        'InstrumentCount',size(moments,2),'WaldStatistic',wald(j), ...
                        'DF',size(moments,2),'P_Joint',pval(j),'P_Holm',NaN, ...
                        'HolmFamilySize',0,'ScaledCovarianceRcond',rc(j),'Status',status{j}, ...
                        'AnalysisRole',analysis_role(lossName,cfg.HACLags(j),cfg.PrimaryHACLag,true));
                end
            end
        end
    end
    T=struct();T.DM=struct2table(dm);T.Effects=struct2table(effects);T.GW=struct2table(gw);
    % Force text to one string per row. struct2table can otherwise collect
    % equal-length character rows into a character matrix, making indexing
    % depend on the lengths of model names and labels.
    T.DM=normalize_text(T.DM,dm,{'ModelA','ModelB','Loss','Status','AnalysisRole'});
    T.Effects=normalize_text(T.Effects,effects,{'ModelA','ModelB','RMSE_Status', ...
        'RelativeRMSE_Status','MAE_Status','AnalysisRole'});
    T.GW=normalize_text(T.GW,gw,{'ModelA','ModelB','Loss','Status','AnalysisRole'});
    T.DM=adjust_families(T.DM,'P_TwoSided');
    if ~isempty(T.GW),T.GW=adjust_families(T.GW,'P_Joint');end
    D=struct('Protocol','Saved-forecast HAC inference v1.0','N',n, ...
        'ModelNames',{names},'Config',rmfield(cfg,'Instruments'), ...
        'MeanSquaredError',mean(sq,1),'MeanAbsoluteError',mean(ab,1), ...
        'ConfidenceLevel',.95,'NormalCriticalValue',z95, ...
        'ForecastsFitted',0,'RandomNumbersUsed',false, ...
        'InferenceScope','Revision-stage exploratory; not preregistered confirmatory inference.', ...
        'EffectDirection','Positive gain means ModelB loss/error exceeds ModelA; positive favors ModelA.', ...
        'CIInterpretation','Pointwise asymptotic 95% CIs; no simultaneous-CI coverage claim.', ...
        'Multiplicity','Holm preserves all requested pairs separately for each test type, loss and bandwidth; undefined p-values count as1 internally and remain NaN in the report.', ...
        'PrimaryReport','MSE paired comparisons at PrimaryHACLag; other losses/bandwidths are secondary/sensitivity.', ...
        'GWInterpretation','Secondary GW-style Wald test of E[h_t*(loss_B-loss_A)]=0 for supplied origin-available instruments only; no direction follows from rejection.', ...
        'Assumptions','Suitable weak dependence, moments, and stability of the relevant loss/moment process; HAC does not resolve arbitrary structural change.', ...
        'Degeneracy','Identical empirical losses have no loss-differential test. Zero-width displayed differences for exact duplicates are descriptive, not proof of equivalence.', ...
        'GWConditioning','Moment columns scaled by positive RMS; covariance Cholesky must succeed and rcond exceed 1e-12. No pseudoinverse.', ...
        'HACConvention','Demeaned moments, common divisor N, Bartlett weights 1-k/(L+1), variance of mean Omega/N; no finite-sample DM correction.');
end

function T=normalize_text(T,records,fields)
    for j=1:numel(fields)
        f=fields{j};T.(f)=reshape(string({records.(f)}),[],1);
    end
end

function [y,pred,names,cfg]=validate_inputs(y,pred,names,opts)
    if ~isnumeric(y)||~isreal(y)||~isvector(y)||numel(y)<3||any(~isfinite(y))||any(y<=0)
        error('VAREN:InferenceTarget','Target must be a finite positive real vector.');end
    y=double(y(:));n=numel(y);
    if ~isnumeric(pred)||~isreal(pred)||~ismatrix(pred)||size(pred,1)~=n|| ...
            size(pred,2)<2||any(~isfinite(pred),'all')
        error('VAREN:InferencePredictions','Predictions must be a finite real N-by-M matrix with M>=2.');end
    pred=double(pred);names=text_cells(names,'model names');names=names(:)';
    if numel(names)~=size(pred,2)||numel(unique(names))~=numel(names)
        error('VAREN:InferenceNames','Supply one unique model name per prediction column.');end
    cfg=struct('HACLags',[168 24 48 720],'PrimaryHACLag',168, ...
        'PairNames',{{'VAREN','StaticEN';'VAREN','FixedEN05';'VAREN','LSBoost24h'}}, ...
        'Instruments',[],'InstrumentNames',{{}});
    if ~isstruct(opts)||~isscalar(opts),error('VAREN:InferenceOptions','opts must be a scalar structure.');end
    fields=fieldnames(opts);
    for j=1:numel(fields)
        f=fields{j};if ~isfield(cfg,f),error('VAREN:InferenceOptions','Unknown option: %s.',f);end
        cfg.(f)=opts.(f);
    end
    L=cfg.HACLags;
    if ~isnumeric(L)||~isreal(L)||~isvector(L)||isempty(L)||any(~isfinite(L))|| ...
            any(L<0)||any(L~=floor(L))||numel(unique(L))~=numel(L)||max(L)>n-2
        error('VAREN:InferenceBandwidth','HACLags must be distinct integer lags from0 toN-2.');end
    cfg.HACLags=double(L(:)');
    if ~isnumeric(cfg.PrimaryHACLag)||~isscalar(cfg.PrimaryHACLag)|| ...
            ~isfinite(cfg.PrimaryHACLag)||~ismember(cfg.PrimaryHACLag,cfg.HACLags)
        error('VAREN:InferenceBandwidth','PrimaryHACLag must be in HACLags.');end
    cfg.PairNames=text_cells(cfg.PairNames,'pair names');
    if size(cfg.PairNames,2)~=2||isempty(cfg.PairNames)
        error('VAREN:InferencePairs','PairNames must be a nonempty two-column cell/string array.');end
    for p=1:size(cfg.PairNames,1)
        if any(~ismember(cfg.PairNames(p,:),names))||strcmp(cfg.PairNames{p,1},cfg.PairNames{p,2})
            error('VAREN:InferencePairs','Each pair must name two distinct existing model columns.');end
        for q=1:p-1
            if isequal(cfg.PairNames(p,:),cfg.PairNames(q,:))||isequal(cfg.PairNames(p,:),fliplr(cfg.PairNames(q,:)))
                error('VAREN:InferencePairs','Repeated or reversed duplicate comparison.');end
        end
    end
    H=cfg.Instruments;
    if isempty(H)
        if ~isempty(cfg.InstrumentNames),error('VAREN:InferenceInstruments','Instrument names supplied without instruments.');end
    else
        if ~isnumeric(H)||~isreal(H)||~ismatrix(H)||size(H,1)~=n|| ...
                size(H,2)>=n||any(~isfinite(H),'all')||any(H(:,1)~=1)
            error('VAREN:InferenceInstruments','Instruments must be finite N-by-K with K<N and first column exactly1.');end
        cfg.Instruments=double(H);cfg.InstrumentNames=text_cells(cfg.InstrumentNames,'instrument names');
        cfg.InstrumentNames=cfg.InstrumentNames(:)';
        if numel(cfg.InstrumentNames)~=size(H,2)||numel(unique(cfg.InstrumentNames))~=size(H,2)
            error('VAREN:InferenceInstruments','Supply a distinct name for each instrument.');end
    end
end

function c=text_cells(c,label)
    if isstring(c),if any(ismissing(c),'all'),error('VAREN:InferenceText','Missing %s.',label);end;c=cellstr(c);end
    if ~iscell(c)||any(~cellfun(@(s)ischar(s)&&isrow(s)&&~isempty(strtrim(s)),c),'all')
        error('VAREN:InferenceText','%s must be nonempty character/string entries.',label);end
end

function [omega,status]=scalar_hac(Z,lags)
    if any(~isfinite(Z),'all'),error('VAREN:InferenceOverflow','Nonfinite influence series.');end
    n=size(Z,1);X=Z-mean(Z,1);scale=max(abs(X),[],1);scale(scale==0)=1;X=X./scale;
    gamma0=mean(X.^2,1);omega=ones(numel(lags),1)*gamma0;bound=omega;
    for k=1:max(lags)
        gamma=sum(X(k+1:end,:).*X(1:end-k,:),1)/n;
        term=2*max(0,1-k./(lags(:)+1))*gamma;
        omega=omega+term;bound=bound+abs(term);
    end
    tol=100*eps*bound;status=repmat({'OK'},size(omega));
    bad=omega< -tol;
    if any(bad,'all'),error('VAREN:InferenceHAC','Bartlett HAC variance is materially negative.');end
    zero=omega<=tol;omega(zero)=0;status(zero)={'ZERO_LONG_RUN_VARIANCE'};
    omega=omega.*(scale.^2);
    if any(~isfinite(omega),'all'),error('VAREN:InferenceOverflow','Nonfinite HAC variance.');end
end

function [se,lo,hi,z,p,status]=scalar_test(d,omega,hacStatus,z95)
    se=NaN;lo=NaN;hi=NaN;z=NaN;p=NaN;status=hacStatus;
    if all(d==0),se=0;lo=0;hi=0;status='IDENTICAL_LOSSES';return;end
    if ~strcmp(status,'OK'),return;end
    se=sqrt(omega/numel(d));gain=mean(d);lo=gain-z95*se;hi=gain+z95*se;
    z=gain/se;p=erfc(abs(z)/sqrt(2));
end

function [se,lo,hi,status]=effect_ci(gain,omega,hacStatus,n,z95,identical)
    se=NaN;lo=NaN;hi=NaN;status=hacStatus;
    if identical,se=0;lo=gain;hi=gain;status='IDENTICAL_LOSSES';return;end
    if ~strcmp(status,'OK'),return;end
    se=sqrt(omega/n);lo=gain-z95*se;hi=gain+z95*se;
end

function [wald,pval,rc,status]=joint_hac(G,lags,identical)
    n=size(G,1);q=size(G,2);J=numel(lags);
    wald=NaN(J,1);pval=NaN(J,1);rc=NaN(J,1);status=repmat({'UNDEFINED'},J,1);
    if identical,status(:)={'IDENTICAL_LOSSES'};return;end
    if any(~isfinite(G),'all'),error('VAREN:InferenceOverflow','Nonfinite conditional moments.');end
    mx=max(abs(G),[],1);mx(mx==0)=1;G=G./mx;
    rmsScale=sqrt(mean(G.^2,1));
    if any(rmsScale==0),status(:)={'SINGULAR_MOMENT_COVARIANCE'};return;end
    G=G./rmsScale;mu=mean(G,1)';X=G-mu';gamma0=(X'*X)/n;
    Omega=repmat(gamma0,1,1,J);
    for k=1:max(lags)
        gamma=(X(k+1:end,:)'*X(1:end-k,:))/n;
        for j=1:J
            if k<=lags(j),Omega(:,:,j)=Omega(:,:,j)+(1-k/(lags(j)+1))*(gamma+gamma');end
        end
    end
    for j=1:J
        V=(Omega(:,:,j)+Omega(:,:,j)')/2;
        if any(~isfinite(V),'all'),error('VAREN:InferenceOverflow','Nonfinite moment covariance.');end
        rc(j)=rcond(V);[R,flag]=chol(V);
        if flag~=0||~isfinite(rc(j))||rc(j)<=1e-12
            status{j}='SINGULAR_OR_ILL_CONDITIONED_COVARIANCE';continue;
        end
        u=R'\mu;wald(j)=n*(u'*u);
        if ~isfinite(wald(j)),status{j}='NONFINITE_WALD';wald(j)=NaN;continue;end
        pval(j)=gammainc(wald(j)/2,q/2,'upper');status{j}='OK';
    end
end

function T=adjust_families(T,pfield)
    lagList=unique(T.HACLag);lossList={'MSE','MAE'};
    for j=1:numel(lagList)
        for k=1:numel(lossList)
            idx=find(T.HACLag==lagList(j)&strcmp(T.Loss,lossList{k}));
            m=numel(idx);T.HolmFamilySize(idx)=m;
            p=T.(pfield)(idx);valid=isfinite(p);p(~valid)=1;
            [p,order]=sort(p);
            adjusted=min(1,cummax((m-(1:m)'+1).*p));
            mapped=NaN(m,1);mapped(order)=adjusted;
            T.P_Holm(idx(valid))=mapped(valid);
        end
    end
end

function role=analysis_role(loss,L,primary,isGW)
    if isGW,role='Secondary conditional-moment exploratory';
    elseif L~=primary,role='Bandwidth sensitivity exploratory';
    elseif strcmp(loss,'MSE'),role='Main exploratory';
    else,role='Secondary loss exploratory';end
    if isGW&&L~=primary,role='Conditional-moment bandwidth sensitivity exploratory';end
end

function r=dm_record()
    r=struct('PairIndex',0,'ModelA','','ModelB','','Loss','','HACLag',0,'N',0, ...
        'MeanLoss_A',NaN,'MeanLoss_B',NaN,'MeanLossGain_BminusA',NaN, ...
        'LongRunVariance',NaN,'SE_MeanLossGain',NaN,'CI95_Lower',NaN,'CI95_Upper',NaN, ...
        'DM_Z',NaN,'P_TwoSided',NaN,'P_Holm',NaN,'HolmFamilySize',0,'Status','', ...
        'IdenticalPredictions',false,'AnalysisRole','');
end

function r=effect_record()
    r=struct('PairIndex',0,'ModelA','','ModelB','','HACLag',0,'N',0, ...
        'RMSE_A_MW',NaN,'RMSE_B_MW',NaN,'RMSE_Gain_MW',NaN,'RMSE_SE_MW',NaN, ...
        'RMSE_CI95_Lower_MW',NaN,'RMSE_CI95_Upper_MW',NaN,'RMSE_Gain_pct',NaN, ...
        'RMSE_SE_pct',NaN,'RMSE_CI95_Lower_pct',NaN,'RMSE_CI95_Upper_pct',NaN, ...
        'MAE_Gain_MW',NaN,'MAE_SE_MW',NaN,'MAE_CI95_Lower_MW',NaN,'MAE_CI95_Upper_MW',NaN, ...
        'RMSE_Status','NONSMOOTH_ZERO_MSE','RelativeRMSE_Status','NONSMOOTH_ZERO_MSE', ...
        'MAE_Status','','AnalysisRole','');
end

function r=gw_record()
    r=struct('PairIndex',0,'ModelA','','ModelB','','Loss','','HACLag',0,'N',0, ...
        'InstrumentCount',0,'WaldStatistic',NaN,'DF',0,'P_Joint',NaN,'P_Holm',NaN, ...
        'HolmFamilySize',0,'ScaledCovarianceRcond',NaN,'Status','','AnalysisRole','');
end
