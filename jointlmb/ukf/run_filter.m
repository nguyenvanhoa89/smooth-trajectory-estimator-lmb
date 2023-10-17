function [est,est_re] = run_filter(model,meas)


%=== Setup

%output variables
est.X= cell(model.K,1);
est.N= zeros(model.K,1);
est.L= cell(model.K,1);

%filter parameters
filter.T_max= 100;                  %maximum number of tracks
filter.track_threshold= 1e-3;       %threshold to prune tracks
filter.H_upd= 10000;                  %requested number of updated components/hypotheses (for GLMB update)

filter.L_max= 10;                   %limit on number of Gaussians in each track 
filter.elim_threshold= 1e-5;        %pruning threshold for Gaussians in each track 
filter.merge_threshold= 4;          %merging threshold for Gaussians in each track

filter.ukf_alpha= 1;                %scale parameter for UKF - choose alpha=1 ensuring lambda=beta and offset of first cov weight is beta for numerical stability
filter.ukf_beta= 2;                 %scale parameter for UKF
filter.ukf_kappa= 2;                %scale parameter for UKF (alpha=1 preferred for stability, giving lambda=1, offset of beta for first cov weight)

filter.P_G= 0.9999999;                           %gate size in percentage
filter.gamma= chi2inv(filter.P_G,model.z_dim);   %inv chi^2 dn gamma value
filter.gate_flag= 1;                             %gating on or off 1/0

filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

est.filter= filter;
est.strategy = 'LMB';
est_re =  init_est_re(meas);


%=== Filtering

%initial prior
tt_lmb_update= cell(0,1);      %track table for LMB (cell array of structs for individual tracks)
filtering_time = zeros(model.K,1);
est_time = zeros(model.K,1);
est_re_time = zeros(model.K,1);

%recursive filtering
for k=1:model.K
    s_filtering_time = tic;
    %joint predict and update, results in GLMB, convert to LMB
    glmb_update= jointlmbpredictupdate(tt_lmb_update,model,filter,meas,k);              T_predict= length(tt_lmb_update)+model.T_birth;
    tt_lmb_update= glmb2lmb(glmb_update,'include_ah',true);                                               T_posterior= length(tt_lmb_update);
       
    %pruning, truncation and track cleanup
    tt_lmb_update= clean_lmb(tt_lmb_update,filter);                                     T_clean= length(tt_lmb_update);
    filtering_time(k) = toc(s_filtering_time);
    %state estimation
    s_est_time = tic;
    [est.X{k},est.N(k),est.L{k}]= extract_estimates(tt_lmb_update,model);
    est_time(k) = toc(s_est_time);
    s_est_re_time = tic;
    est_re=extract_lmb_estimates_smooth(model,tt_lmb_update,est_re,meas);
    est_re_time(k) = toc(s_est_re_time);

    %display diagnostics
    display_diaginfo(tt_lmb_update,k,est,filter,T_predict,T_posterior,T_clean);
    
end
est.est_time = est_time;
est.filtering_time = filtering_time;
est_re.est_time = est_re_time;
est_re.filtering_time = filtering_time;
end



function glmb_nextupdate= jointlmbpredictupdate(tt_lmb_update,model,filter,meas,k)
%---generate birth tracks
tt_birth= cell(length(model.r_birth),1);                                           %initialize cell array
for tabbidx=1:length(model.r_birth)
    tt_birth{tabbidx}.r= model.r_birth(tabbidx);                                   %birth prob for birth track
    tt_birth{tabbidx}.m= model.m_birth{tabbidx};                                   %means of Gaussians for birth track
    tt_birth{tabbidx}.P= model.P_birth{tabbidx};                                   %covs of Gaussians for birth track
    tt_birth{tabbidx}.w= model.w_birth{tabbidx}(:);                                %weights of Gaussians for birth track
    tt_birth{tabbidx}.l= [k;tabbidx];                                              %track label
    tt_birth{tabbidx}.ah= [];
end

%---generate surviving tracks
tt_survive= cell(length(tt_lmb_update),1);                                                                              %initialize cell array
for tabsidx=1:length(tt_lmb_update)
    tt_survive{tabsidx}.r= model.P_S*tt_lmb_update{tabsidx}.r;                                                          %predicted existence probability for surviving track
    [mtemp_predict,Ptemp_predict]= model.filter.predict_func(model,tt_lmb_update{tabsidx}.m,tt_lmb_update{tabsidx}.P);      %kalman prediction for GM
    tt_survive{tabsidx}.m= mtemp_predict;                                                                               %means of Gaussians for surviving track
    tt_survive{tabsidx}.P= Ptemp_predict;                                                                               %covs of Gaussians for predicted track
    tt_survive{tabsidx}.w= tt_lmb_update{tabsidx}.w;                                                                    %weights of Gaussians for predicted track
    tt_survive{tabsidx}.l= tt_lmb_update{tabsidx}.l;                                                                    %track label
    tt_survive{tabsidx}.ah= tt_lmb_update{tabsidx}.ah;
end

%create predicted tracks - concatenation of birth and survival
tt_predict= cat(1,tt_birth,tt_survive);                                                                                %copy track table back to GLMB struct

%gating by tracks
if filter.gate_flag
    for tabidx=1:length(tt_predict)
        tt_predict{tabidx}.gatemeas= model.filter.gate_func(meas.Z{k},filter.gamma,model,tt_predict{tabidx}.m,tt_predict{tabidx}.P,filter.ukf_alpha,filter.ukf_kappa,filter.ukf_beta);
    end
else
    for tabidx=1:length(tt_predict)
        tt_predict{tabidx}.gatemeas= 1:size(meas.Z{k},2);
    end
end

%precalculation loop for average survival/death probabilities
avps= zeros(length(tt_predict),1);
for tabidx=1:length(tt_predict)
    avps(tabidx)= tt_predict{tabidx}.r;
end
avqs= 1-avps;

%precalculation loop for average detection/missed probabilities
avpd= zeros(length(tt_predict),1);
for tabidx=1:length(tt_predict)
    avpd(tabidx)= model.P_D;
end
avqd= 1-avpd;

%create updated tracks (single target Bayes update)
m= size(meas.Z{k},2);                                   %number of measurements
tt_update= cell((1+m)*length(tt_predict),1);            %initialize cell array
%missed detection tracks (legacy tracks)
for tabidx= 1:length(tt_predict)
    tt_update{tabidx}= tt_predict{tabidx};              %same track table
    tt_update{tabidx}.ah = [tt_update{tabidx}.ah; 0];   %track association history (updated for missed detection)
end
%measurement updated tracks (all pairs)
allcostm= zeros(length(tt_predict),m);
for tabidx= 1:length(tt_predict)
    for emm= tt_predict{tabidx}.gatemeas
            stoidx= length(tt_predict)*emm + tabidx; %index of predicted track i updated with measurement j is (number_predicted_tracks*j + i)
            [qz_temp,m_temp,P_temp] = model.filter.update_func(meas.Z{k}(:,emm),model,tt_predict{tabidx}.m,tt_predict{tabidx}.P);   %kalman update for this track and this measurement
            w_temp= qz_temp.*tt_predict{tabidx}.w+eps;                                                                                      %unnormalized updated weights
            tt_update{stoidx}.m= m_temp;                                                                                                    %means of Gaussians for updated track
            tt_update{stoidx}.P= P_temp;                                                                                                    %covs of Gaussians for updated track
            tt_update{stoidx}.w= w_temp/sum(w_temp);                                                                                        %weights of Gaussians for updated track
            tt_update{stoidx}.l = tt_predict{tabidx}.l;                                                                                     %track label
            tt_update{stoidx}.ah = [tt_predict{tabidx}.ah; emm];                                                                            %assign measurement association
            allcostm(tabidx,emm)= sum(w_temp);                                                                                              %predictive likelihood
    end
end
glmb_nextupdate.tt= tt_update;                                                                                                              %copy track table back to GLMB struct
%joint cost matrix
jointcostm= [diag(avqs) ...
             diag(avps.*avqd) ...
             repmat(avps.*avpd,[1 m]).*allcostm/(model.lambda_c*model.pdf_c)];
%gated measurement index matrix
gatemeasidxs= zeros(length(tt_predict),m);
for tabidx= 1:length(tt_predict)
    gatemeasidxs(tabidx,1:length(tt_predict{tabidx}.gatemeas))= tt_predict{tabidx}.gatemeas;
end
gatemeasindc= gatemeasidxs>0;
         

%component updates

    %calculate best updated hypotheses/components
    cpreds= length(tt_predict);
    nbirths= model.T_birth;
    nexists= length(tt_lmb_update);
    ntracks= nbirths + nexists;
    tindices= [(1:nbirths) nbirths+(1:nexists)];                                                                                          %indices of all births and existing tracks  for current component
    lselmask= false(length(tt_predict),m); lselmask(tindices,:)= gatemeasindc(tindices,:);                                              %logical selection mask to index gating matrices
    mindices= unique_faster(gatemeasidxs(lselmask));                                                                                    %union indices of gated measurements for corresponding tracks
    costm= jointcostm(tindices,[tindices cpreds+tindices 2*cpreds+mindices]);                                                           %cost matrix - [no_birth/is_death | born/survived+missed | born/survived+detected]
    neglogcostm= -log(costm);                                                                                                           %negative log cost
    [uasses,nlcost]= gibbswrap_jointpredupdt_custom(neglogcostm,round(filter.H_upd));                                                   %murty's algo/gibbs sampling to calculate m-best assignment hypotheses/components
    uasses(uasses<=ntracks)= -inf;                                                                                                      %set not born/track deaths to -inf assignment
    uasses(uasses>ntracks & uasses<= 2*ntracks)= 0;                                                                                     %set survived+missed to 0 assignment
    uasses(uasses>2*ntracks)= uasses(uasses>2*ntracks)-2*ntracks;                                                                       %set survived+detected to assignment of measurement index from 1:|Z|    
    uasses(uasses>0)= mindices(uasses(uasses>0));                                                                                       %restore original indices of gated measurements
    
    %generate corrresponding jointly predicted/updated hypotheses/components
    for hidx=1:length(nlcost)
        update_hypcmp_tmp= uasses(hidx,:)'; 
        update_hypcmp_idx= cpreds.*update_hypcmp_tmp+[(1:nbirths)'; nbirths+(1:nexists)'];
        glmb_nextupdate.w(hidx)= -model.lambda_c+m*log(model.lambda_c*model.pdf_c)-nlcost(hidx);                                             %hypothesis/component weight
        glmb_nextupdate.I{hidx}= update_hypcmp_idx(update_hypcmp_idx>0);                                                                                              %hypothesis/component tracks (via indices to track table)
        glmb_nextupdate.n(hidx)= sum(update_hypcmp_idx>0);                                                                                                            %hypothesis/component cardinality
    end

glmb_nextupdate.w= exp(glmb_nextupdate.w-logsumexp(glmb_nextupdate.w));                                                                                                                 %normalize weights

%extract cardinality distribution
for card=0:max(glmb_nextupdate.n)
    glmb_nextupdate.cdn(card+1)= sum(glmb_nextupdate.w(glmb_nextupdate.n==card));                                                                                                       %extract probability of n targets
end

%remove duplicate entries and clean track table
glmb_nextupdate= clean_update(clean_predict(glmb_nextupdate));
end



function glmb_temp= clean_predict(glmb_raw)
%hash label sets, find unique ones, merge all duplicates
for hidx= 1:length(glmb_raw.w)
    glmb_raw.hash{hidx}= sprintf('%i*',sort(glmb_raw.I{hidx}(:)'));
end

[cu,~,ic]= unique(glmb_raw.hash);

glmb_temp.tt= glmb_raw.tt;
glmb_temp.w= zeros(length(cu),1);
glmb_temp.I= cell(length(cu),1);
glmb_temp.n= zeros(length(cu),1);
for hidx= 1:length(ic)
        glmb_temp.w(ic(hidx))= glmb_temp.w(ic(hidx))+glmb_raw.w(hidx);
        glmb_temp.I{ic(hidx)}= glmb_raw.I{hidx};
        glmb_temp.n(ic(hidx))= glmb_raw.n(hidx);
end
glmb_temp.cdn= glmb_raw.cdn;
end



function glmb_clean= clean_update(glmb_temp)
%flag used tracks
usedindicator= zeros(length(glmb_temp.tt),1);
for hidx= 1:length(glmb_temp.w)
    usedindicator(glmb_temp.I{hidx})= usedindicator(glmb_temp.I{hidx})+1;
end
trackcount= sum(usedindicator>0);

%remove unused tracks and reindex existing hypotheses/components
newindices= zeros(length(glmb_temp.tt),1); newindices(usedindicator>0)= 1:trackcount;
glmb_clean.tt= glmb_temp.tt(usedindicator>0);
glmb_clean.w= glmb_temp.w;
for hidx= 1:length(glmb_temp.w)
    glmb_clean.I{hidx}= newindices(glmb_temp.I{hidx});
end
glmb_clean.n= glmb_temp.n;
glmb_clean.cdn= glmb_temp.cdn;
end



function tt_lmb= glmb2lmb(glmb, varargin)
    %  Convert from GLMB density to LMB density
    
    % --- Input Parser
    p = inputParser;
    addParameter(p,'include_ah',false);                                                         %track threshold
    parse(p, varargin{:});
    include_ah = p.Results.include_ah;

    %find unique labels (with different possibly different association histories)
    lmat= zeros(2,length(glmb.tt),1);
    for tabidx= 1:length(glmb.tt)
        lmat(:,tabidx)= glmb.tt{tabidx}.l;
    end
    lmat= lmat';
    
    [cu,~,ic]= unique(lmat,'rows'); cu= cu';
    
    %initialize LMB struct
    tt_lmb= cell(size(cu,2),1);
    for tabidx=1:length(tt_lmb)
        tt_lmb{tabidx}.r= 0;
        tt_lmb{tabidx}.m= [];
        tt_lmb{tabidx}.P= [];
        tt_lmb{tabidx}.w= [];
        tt_lmb{tabidx}.l= cu(:,tabidx);
        if include_ah
            tt_lmb{tabidx}.ah= []; 
            tt_lmb{tabidx}.trkidx= []; 
        end
    end
    
    %extract individual tracks
    for hidx=1:length(glmb.w)
        for t= 1:glmb.n(hidx)
            trkidx= glmb.I{hidx}(t);
            newidx= ic(trkidx);
            tt_lmb{newidx}.m= cat(2,tt_lmb{newidx}.m,glmb.tt{trkidx}.m);
            tt_lmb{newidx}.P= cat(3,tt_lmb{newidx}.P,glmb.tt{trkidx}.P);
            tt_lmb{newidx}.w= cat(1,tt_lmb{newidx}.w,glmb.w(hidx)*glmb.tt{trkidx}.w);
            if include_ah
                tt_lmb{newidx}.ah= cat(1,tt_lmb{newidx}.ah,glmb.tt{trkidx}.ah(end)); 
                tt_lmb{newidx}.trkidx= cat(1,tt_lmb{newidx}.trkidx,trkidx*ones(length(glmb.tt{trkidx}.w),1)); 
            end
        end
    end
    
    %extract the maximum association history for each track
    if include_ah
        for tabidx=1:length(tt_lmb)
            % 1st method
            [~,sel_ah_idx] = max(tt_lmb{tabidx}.w);
            try
                tt_lmb{tabidx}.ah = glmb.tt{tt_lmb{tabidx}.trkidx(sel_ah_idx)}.ah;
            catch err
               disp(err.message);
            end
            tt_lmb{tabidx}.trkidx= []; 

        end
    end
    
    %extract existence probabilities and normalize track weights
    for tabidx=1:length(tt_lmb)
        tt_lmb{tabidx}.r= sum(tt_lmb{tabidx}.w);
        tt_lmb{tabidx}.w= tt_lmb{tabidx}.w/tt_lmb{tabidx}.r;
        tt_lmb{tabidx}.r = limit_range(tt_lmb{tabidx}.r);
    end
    
end




function tt_lmb_out= clean_lmb(tt_lmb_in,filter)
%prune tracks with low existence probabilities
rvect= get_rvals(tt_lmb_in);
idxkeep= find(rvect > filter.track_threshold);
rvect= rvect(idxkeep);
tt_lmb_out= tt_lmb_in(idxkeep);

%enforce cap on maximum number of tracks
if length(tt_lmb_out) > filter.T_max
    [~,idxkeep]= sort(rvect,'descend');
    tt_lmb_out= tt_lmb_out(idxkeep);   
end

%cleanup tracks
for tabidx=1:length(tt_lmb_out)
    [tt_lmb_out{tabidx}.w,tt_lmb_out{tabidx}.m,tt_lmb_out{tabidx}.P]= gaus_prune(tt_lmb_out{tabidx}.w,tt_lmb_out{tabidx}.m,tt_lmb_out{tabidx}.P,filter.elim_threshold);
    [tt_lmb_out{tabidx}.w,tt_lmb_out{tabidx}.m,tt_lmb_out{tabidx}.P]= gaus_merge(tt_lmb_out{tabidx}.w,tt_lmb_out{tabidx}.m,tt_lmb_out{tabidx}.P,filter.merge_threshold);
    [tt_lmb_out{tabidx}.w,tt_lmb_out{tabidx}.m,tt_lmb_out{tabidx}.P]= gaus_cap(tt_lmb_out{tabidx}.w,tt_lmb_out{tabidx}.m,tt_lmb_out{tabidx}.P,filter.L_max);
end
end



function rvect= get_rvals(tt_lmb)                           %function to extract vector of existence probabilities from LMB track table
rvect= zeros(length(tt_lmb),1);
for tabidx=1:length(tt_lmb)
   rvect(tabidx)= tt_lmb{tabidx}.r; 
end
end



function [X,N,L]=extract_estimates(tt_lmb,model)
%extract estimates via MAP cardinality and corresponding tracks
rvect= get_rvals(tt_lmb); rvect= min(rvect,0.999); rvect= max(rvect,0.001);
cdn= prod(1-rvect)*esf(rvect./(1-rvect));
[~,mode] = max(cdn);
N = min(length(rvect),mode-1);
X= zeros(model.x_dim,N);
L= zeros(2,N);

[~,idxcmp]= sort(rvect,'descend');
for n=1:N
    [~,idxtrk]= max(tt_lmb{idxcmp(n)}.w);
    X(:,n)= tt_lmb{idxcmp(n)}.m(:,idxtrk);
    L(:,n)= tt_lmb{idxcmp(n)}.l;
end
end



function display_diaginfo(tt_lmb,k,est,filter,T_predict,T_posterior,T_clean)
rvect= get_rvals(tt_lmb); rvect= min(rvect,0.999); rvect= max(rvect,0.001);
cdn= prod(1-rvect)*esf(rvect./(1-rvect));
eap= (0:(length(cdn)-1))*cdn(:);
var= (0:(length(cdn)-1)).^2*cdn(:)-((0:(length(cdn)-1))*cdn(:))^2;
if ~strcmp(filter.run_flag,'silence')
    disp([' time= ',num2str(k),...
        ' #eap cdn=' num2str(eap),...
        ' #var cdn=' num2str(var,4),...
        ' #est card=' num2str(est.N(k),4),...
        ' #trax pred=' num2str(T_predict,4),...
        ' #trax post=' num2str(T_posterior,4),...
        ' #trax updt=',num2str(T_clean,4)   ]);
end
end


function clipped_r= limit_range(r)
    r(r>=1-eps)=1-eps;
    r(r<eps)=eps;
    clipped_r= r;
end

function est =  init_est_re(model)
    est.X= cell(model.K,1);
    est.N= zeros(model.K,1);
    est.L= cell(model.K,1);
    est.P= cell(model.K,1);
    est.T= {};
    est.M= 0;
    est.J= []; est.H= {};
    est.strategy = 'STE-LMB';
end

function est=extract_lmb_estimates_smooth(model,tt_lmb,est,meas)
%extract estimates via recursive estimator, where  
%trajectories are extracted via association history, and
%track continuity is guaranteed with a non-trivial estimator
K = model.K;
prune_flag = 0;
prune_thres = 3;
%extract estimates via MAP cardinality and corresponding tracks
rvect= get_rvals(tt_lmb); rvect= min(rvect,0.999); rvect= max(rvect,0.001);
cdn= prod(1-rvect)*esf(rvect./(1-rvect));
[~,mode] = max(cdn);
M = min(length(rvect),mode-1);
[~,idxcmp]= sort(rvect,'descend');
% prepare the association history and its labels
T= cell(M,1);
J= zeros(2,M);
for m=1:M
    T{m,1}= tt_lmb{idxcmp(m)}.ah;
    J(:,m)= tt_lmb{idxcmp(m)}.l;
end

H= cell(M,1);
for m=1:M
   H{m}= [num2str(J(1,m)),'.',num2str(J(2,m))]; 
end

% compute dead & updated & new tracks
[~,io,is]= intersect(est.H,H);
[~,id,in]= setxor(est.H,H);

est.M= M;
est.T= cat(1,est.T(id),T(is),T(in));
est.J= cat(2,est.J(:,id),J(:,is),J(:,in));
est.H= cat(1,est.H(id),H(is),H(in));

%write out estimates in standard format
est.N= zeros(meas.K,1);
est.X= cell(meas.K,1);
est.L= cell(meas.K,1);
for t=1:length(est.T)
    tah= est.T{t};
    traj_length = length(tah) ; 
    ks= est.J(1,t);
    k_end = ks+traj_length-1 ;
    if prune_flag
        if k_end<K && traj_length<=prune_thres
            continue
        end
    end
    bidx= est.J(2,t);
    
    m_cell = cell(length(tah),1) ; 
    P_cell = cell(length(tah),1) ; 
    w_cell = cell(length(tah),1) ; 

    % assign birth
    m_cell{1} = model.m_birth{bidx};
    P_cell{1} = model.P_birth{bidx};
    w_cell{1} = model.w_birth{bidx};
    % forward filtering of trajectory
    for u=2:length(tah)
        
        [m_cell{u},P_cell{u}] = model.filter.predict_func(model,m_cell{u-1},P_cell{u-1});
        w_cell{u} = w_cell{u-1} ; 
        k= ks+u-1;
        emm= tah(u);
        if emm > 0
            [qz,m_cell{u},P_cell{u}] = model.filter.update_func(meas.Z{k}(:,emm),model,m_cell{u},P_cell{u});
            w_cell{u}= qz.*w_cell{u}+eps;
            w_cell{u}= w_cell{u}/sum(w_cell{u});
        end

    end
    % backward smoothing of trajectory 
    [~,idxtrk] = max(w_cell{traj_length}) ; 
    est.N(k_end)= est.N(k_end)+1;
    est.X{k_end}= cat(2,est.X{k_end},m_cell{traj_length}(:,idxtrk));
    est.L{k_end}= cat(2,est.L{k_end},est.J(:,t));
    for u= length(tah)-1 : -1 : 1
        m = model.filter.smooth_func(model,w_cell{u},m_cell{u},P_cell{u}, w_cell{u+1},m_cell{u+1}) ;
        k= ks+u-1;
        est.N(k)= est.N(k)+1;
        est.X{k}= cat(2,est.X{k},m);
        est.L{k}= cat(2,est.L{k},est.J(:,t));
    end
end
end

function est = extract_lmb_estimates_recursive(model,tt_lmb,est,meas)
    % Extract the LMB estimates recursively, which change the estimate
    % history from when each object is born until time k_in

    % make code look nicer
    meas_type = model.source_info.meas_type;
    R = model.R;
    

    %extract estimates via MAP cardinality and corresponding tracks
    rvect= get_rvals(tt_lmb); rvect= min(rvect,0.999); rvect= max(rvect,0.001);
    cdn= prod(1-rvect)*esf(rvect./(1-rvect));
    [~,mode] = max(cdn);
    M = min(length(rvect),mode-1);
    [~,idxcmp]= sort(rvect,'descend');

    % prepare the association history and its labels
    T= cell(M,1);
    J= zeros(2,M);
    for m=1:M
        T{m,1}= tt_lmb{idxcmp(m)}.ah;
        J(:,m)= tt_lmb{idxcmp(m)}.l;
    end

    H= cell(M,1);
    for m=1:M
       H{m}= [num2str(J(1,m)),'.',num2str(J(2,m))]; 
    end

    % compute dead & updated & new tracks
    [~,io,is]= intersect(est.H,H);
    [~,id,in]= setxor(est.H,H);

    est.M= M;
    est.T= cat(1,est.T(id),T(is),T(in));
    est.J= cat(2,est.J(:,id),J(:,is),J(:,in));
    est.H= cat(1,est.H(id),H(is),H(in));

    % write out estimates in standard format
    est.N= zeros(model.K,1);
    est.X= cell(model.K,1);
    est.P= cell(model.K,1);
    est.L= cell(model.K,1);
    
    %% main program
    for t=1:length(est.T)
        ks= est.J(1,t);
        bidx= est.J(2,t); % birth idx from measurements
        
        tah= est.T{t};                      % extract association history
        
        % compute initial birth value at time ks
        w = 1;
        if ks>1 && ~isempty(meas.Z{ks}) && model.abp.enable
            try
                [m,P] = convert_z_to_x(meas.Z{ks-1}(:,bidx),meas_type,model.source_info.source_pos,R);
                if size(m,1)<model.x_dim
                    m = [m;0];
                    P = diag([diag(P);model.P_birth{1}(5,5)]);
                end
            catch err
                disp(err.message);
                m= model.m_birth{bidx};
                P= model.P_birth{bidx};
                w= model.w_birth{bidx};
            end
        else
            m= model.m_birth{bidx};
            P= model.P_birth{bidx};
            w= model.w_birth{bidx};
        end

        % recompute the whole trajectory
        tah_len = length(tah);
%         for k=ks:1:min(k_in,find(tah >= 0 , 1, 'last' ))
        for i = 1 : tah_len
            k = ks + i - 1; 
            try
            [m,P] = model.filter.predict_func(model,m,P);
            catch err
                disp(err.message);
            end
            emm= tah(i);
            if emm > 0
                try
                [qz,m,P] = model.filter.update_func(meas.Z{k}(:,emm),model,m,P);
                catch err
                   disp(err.message); 
                end
                w= qz.*w+eps;
                w= w/sum(w);
            end

            [~,idxtrk]= max(w);
            try
            est.N(k)= est.N(k)+1;
            est.X{k}= cat(2,est.X{k},m(:,idxtrk));
            est.P{k}= cat(2,est.P{k},diag(P(:,:,idxtrk)));
            est.L{k}= cat(2,est.L{k},est.J(:,t));
            catch err
                disp(err.message);
            end
        end
    end
    
end

function [x,P] = convert_z_to_x(z,meas_type,source_pos,R)
    switch meas_type
        case 'pos'
            if size(R,1) == 2
                x = [z(1);0;z(2);0]; 
                P = eye(4);
                P(1,1) = 9*R(1,1);
                P(2,2) = 4*R(1,1);
                P(3,3) = 9*R(2,2);
                P(4,4) = 4*R(2,2);
            elseif  size(R,1) == 3
                x = [z(1);0;z(2);0;z(3);0]; 
                P = eye(6);
                P(1,1) = 9*R(1,1);
                P(2,2) = 4*R(1,1);
                P(3,3) = 9*R(2,2);
                P(4,4) = 4*R(2,2);
                P(5,5) = 9*R(3,3);
                P(6,6) = 4*R(3,3);
            end
        case 'brg_rng'
            x = zeros(4,1);
            x(1) = source_pos(1) + z(2) * sin(z(1));
            x(3) = source_pos(2) + z(2) * cos(z(1));
            
            P = eye(4);
            noise_range = sqrt(R(1,1)); 
            noise_azimuth = sqrt(R(2,2));
            P(1,1) = (3*noise_range * sin(noise_azimuth))^2;
            P(2,2) = (2*noise_range * sin(noise_azimuth))^2;
            P(3,3) = (3*noise_range * cos(noise_azimuth))^2;
            P(4,4) = (2*noise_range * cos(noise_azimuth))^2;
        case 'rng_rngrt_brg'
            x = zeros(4,1);
            x(1) = source_pos(1) + z(1) * sin(z(3));
            x(3) = source_pos(2) + z(1) * cos(z(3));
            
            P = eye(4);
            noise_range = sqrt(R(1,1)); 
            
            noise_azimuth = sqrt(R(3,3));
            P(1,1) = (3*noise_range * sin(noise_azimuth))^2;
            P(2,2) = 4*R(2,2);
            P(3,3) = (3*noise_range * cos(noise_azimuth))^2;
            P(4,4) = 4*R(2,2);
        case 'radar3d'
            x = zeros(6,1); % https://au.mathworks.com/help/matlab/ref/sph2cart.html
            
            x(1) = source_pos(1) + z(1) * cos(z(3)) * sin(z(4));
            x(3) = source_pos(2) + z(1) * cos(z(3)) * cos(z(4));
            x(5) = source_pos(3) + z(1) * sin(z(3));
            
            P = eye(6);
            noise_range = sqrt(R(1,1)); 
            
            noise_elevation = sqrt(R(3,3));
            noise_azimuth = sqrt(R(4,4));
            P(1,1) = (3*noise_range * cos(noise_elevation) * sin(noise_azimuth))^2;
            P(2,2) = 4*R(2,2);
            P(3,3) = (3*noise_range * cos(noise_elevation) * cos(noise_azimuth))^2;
            P(4,4) = 4*R(2,2);
            P(5,5) = (3*noise_range * sin(noise_elevation))^2;
            P(6,6) = 4*R(2,2);            
    end
end

