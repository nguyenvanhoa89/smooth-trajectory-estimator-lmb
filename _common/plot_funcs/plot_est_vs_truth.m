function plot_est_vs_truth(model, truth, est,varargin)
    % Plot estimates versus ground truth using OSPA2 to match the colors of estimates to colors of ground truth
    % via model, truth and a variable-length input argument list, i.e., a pair of (property, value).

    % ---Input Parser
    p = inputParser;
    % Setup parsing schema
    addParameter(p,'K', truth.K, @isnumeric);
    addParameter(p,'Transparency', 1, @isnumeric);
    addParameter(p,'colorarray',[]);
    parse(p, varargin{:});
    
    % --- Make code look nicer
    Transparency = p.Results.Transparency;

    colorarray = p.Results.colorarray;
    if isempty(colorarray)
       colorarray= makecolorarray(100); 
    end
    LineWidth = 1;
    text_offset = 50;
    font_name = 'Times New Roman';
    font_size = 14;
    
    K = p.Results.K;
    model.pos_idx = [1 3];
    % --- plot truth
    color_list = color_vector(10000)';
    [X_track,k_birth,k_death]= extract_tracks(truth.X,truth.track_list,truth.total_tracks);
    ntarget = truth.total_tracks;
    for i = 1 : ntarget                                                                         % assign truth label id first
       [~,colorarray]= assigncolor(i,colorarray);
    end
    figure(); hold on;

    
    htruth = cell(ntarget,1);
    for i=1:ntarget
        k_b_temp = k_birth(i); k_b_temp = k_b_temp(k_b_temp<=K);                                % update birth time
        k_d_temp = k_death(i); k_d_temp = min(k_d_temp,K);                                      % update death time
        life_temp = k_b_temp : k_d_temp;
        pos_temp = X_track(model.pos_idx,:,i);
        cur_color = 'k';
        Transparency_temp = Transparency;                                                       % to make the whole truth blure
        if ~isempty(k_b_temp)
                htruth{i} = plot(pos_temp(1,life_temp),pos_temp(2,life_temp),'LineWidth',LineWidth, 'LineStyle','-','Color' , cur_color);
                htruth{i}.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("k",life_temp);
                htruth{i}.Color(4) = Transparency_temp;
                plot_temp = scatter(pos_temp(1,k_b_temp),pos_temp(2,k_b_temp),100, 'LineWidth',1,...
                    'Marker' , 'o', 'MarkerFaceColor', cur_color, 'MarkerEdgeColor', 'black');
                plot_temp.MarkerEdgeAlpha =  Transparency_temp; plot_temp.MarkerFaceAlpha =  Transparency_temp;
                if K >= k_death(i)
                    plot_temp = scatter(pos_temp(1,k_d_temp),pos_temp(2,k_d_temp),100, 'LineWidth',1,...
                        'Marker' , 's', 'MarkerFaceColor', cur_color, 'MarkerEdgeColor', 'black');
                end
                plot_temp.MarkerEdgeAlpha =  Transparency_temp; plot_temp.MarkerFaceAlpha =  Transparency_temp;
        end
    end
    

    
    % --- plot est
    [Y_track,l_list,ke_birth,ke_death]= extract_tracks_with_labels(model,est,1,K);
    n_est = size(l_list,2);

    % For color matching between truth and est

    [~,allcostm]=compute_ospa2(X_track([1 3],:,:),Y_track([1 3],:,:),model.ospa.c,model.ospa.p,K);
    if size(allcostm,2) ~= n_est
        allcostm = allcostm';
    end
    % --- List of matched and unmatched tracks
    Matching = Hungarian(allcostm);
    cost_check = allcostm < model.ospa.c;
    Matching = Matching .* cost_check;

    l1_idx = (1 : ntarget)';
    l2_idx = (1 : n_est)';
    L1_idx =  Matching * l2_idx;
    Q = [l1_idx,L1_idx];
    Q_check = prod(Q>0,2)>0;
    Q = Q(Q_check,:); % List of all matched track
    
    
    hest = cell(n_est,1);
    for i = 1 : n_est
        pos_temp = Y_track(model.pos_idx,:,i);
        k_b_temp = ke_birth(i); k_b_temp = k_b_temp(k_b_temp<=K);                                   % update birth time
        k_d_temp = ke_death(i); k_d_temp = min(k_d_temp,K);                                         % update death time
        life_temp = k_b_temp : k_d_temp;
        cur_l = convert_object_label_to_string(l_list(:,i));
        
        if ismember(i,Q(:,2))
            truth_idx = Q(i == Q(:,2),1);
            cur_color = colorarray.rgb(assigncolor(truth_idx,colorarray),:)' ;
        else
            cur_color =  color_list(:,i+n_est);
        end
        if K > k_d_temp, Transparency_temp = Transparency; else, Transparency_temp = 1; end
        if ~isempty(k_b_temp)
            cur_l_list = repmat(convert_object_label_to_string(l_list(:,i)),1,length(life_temp));
            hest{i} = plot(pos_temp(1,life_temp),pos_temp(2,life_temp), '.','Color', cur_color, 'LineWidth',2,'MarkerSize',15); hold on;
            hest{i}.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("ID",repmat(convert_object_label_to_string(l_list(:,i)),1,length(life_temp)));
            hest{i}.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("k",life_temp);
            hest{i}.Color(4) = Transparency_temp;
        end
    end
    
    % --- Format plot
    title(['Truth vs Estimated - ',est.strategy]);
    xlabel('x-coordinate (m)', 'FontSize', font_size);
    ylabel('y-coordinate (m)', 'FontSize', font_size);
    axis equal
    xlim(model.limit(1,:));
    ylim(model.limit(2,:) );     
    set(gcf,'color','w');
    set(gca, 'FontSize', font_size, 'FontName', font_name);
    grid on;

    
end

function [idx,colorarray]= assigncolor(label,colorarray)
        str= sprintf('%i*',label);
        tmp= strcmp(str,colorarray.lab);
        if any(tmp)
            idx= find(tmp);
        else
            colorarray.cnt= colorarray.cnt + 1;
            colorarray.lab{colorarray.cnt}= str;
            idx= colorarray.cnt;
        end
    end

function [X_track,k_birth,k_death, P_track]= extract_tracks(X,track_list,total_tracks,varargin)
    % Extract tracks from the estimated tracks
    
    % --- Input Parser
    pRaw = inputParser;
    addParameter(pRaw,'P',[]);      
    parse(pRaw, varargin{:});     
    P = pRaw.Results.P;
    
    K= size(X,1);
    x_dim= size(X{K},1);
    k=K-1; while x_dim==0, x_dim= size(X{k},1); k= k-1; end
    X_track= NaN(x_dim,K,total_tracks);
    k_birth= zeros(total_tracks,1);
    k_death= zeros(total_tracks,1);
    P_track = [];
    if nargout == 4 && ~isempty(P)
        P_track= NaN(x_dim,K,total_tracks);
    end
    
    max_idx= 0;
    for k=1:K
        if ~isempty(X{k})
            X_track(:,k,track_list{k})= X{k};
            if nargout == 4 && ~isempty(P)
                P_track(:,k,track_list{k})= P{k};
            end
        end
        if max(track_list{k})> max_idx %new target born?
            idx= find(track_list{k}> max_idx);
            k_birth(track_list{k}(idx))= k;
        end
        if ~isempty(track_list{k}), max_idx= max([track_list{k},max_idx]); end
        k_death(track_list{k})= k;
    end
end

function [Y_track,l_list,k_birth,k_death,Y_P_track]= extract_tracks_with_labels(model,est,cur_j,cur_k,varargin)
    % extract_tracks_with_labels
    % % --- Input Parser
    p = inputParser;
    addParameter(p,'use_cov',false);                                                     % use covariance or not
    parse(p, varargin{:});
    use_cov = p.Results.use_cov;
    try
        sel_idx = cur_j : cur_k;
        l_list = unique([cell2mat(est.L(sel_idx)')]','rows','stable')';
        if isempty(l_list)
           l_list = zeros(2,0); 
        end
        labelcount= size(l_list,2);
        l_idx = 1 : labelcount;
        
        est.total_tracks= labelcount;
        est.track_list= cell(length(sel_idx),1);
        idx = 0;
        for k=cur_j:cur_k
            if ~isempty(l_idx)
                idx = idx + 1;
                n_l = size(est.L{k},2);
                for j = 1 : n_l
                    cur_l = est.L{k}(:,j);
                    cur_l_idx = l_idx(ismember(l_list',cur_l','rows')');
                    est.track_list{idx} = [est.track_list{idx},cur_l_idx];
                end
            end
        end
        if labelcount == 0
            Y_track = zeros(model.x_dim,cur_k-cur_j+1,0);
            k_birth = [];
            k_death = [];
            Y_P_track = zeros(model.x_dim,cur_k-cur_j+1,0);
        elseif use_cov
            [Y_track,k_birth,k_death, Y_P_track]= extract_tracks(est.X(sel_idx),est.track_list,est.total_tracks,'P', est.P(sel_idx));
        else
            [Y_track,k_birth,k_death, Y_P_track]= extract_tracks(est.X(sel_idx),est.track_list,est.total_tracks);
        end
    catch err
        error(err.message);
    end
end

function [result,trk_dist] = compute_ospa2(X,Y,c,p,wl)

% This is the MATLAB code for the implementation OSPA(2) metric proposed in
% M. Beard, B.-T. Vo, and B.-N. Vo, "Performance Evaluation for Large-Scale Multi-Target Tracking Algorithms," Proc. 21st IEEE Intl. Conf. Information Fusion, Jul. 2018, Cambridge, UK.
% http://ba-ngu.vo-au.com/vo/BVV_OSPA2_FUSION18.pdf
% and
% M. Beard, B.-T. Vo, and B.-N. Vo, "A Solution for Large-Scale Multi-Object Tracking," arXiv preprint, arXiv:1804.06622, Apr. 2018.
% https://arxiv.org/abs/1804.06622
% based on the OSPA metric proposed in
% D. Schuhmacher, B.-T. Vo, and B.-N. Vo, "A consistent metric for performance evaluation in multi-object filtering," IEEE Trans. Signal Processing, Vol. 56, No. 8 Part 1, pp. 3447� 3457, 2008.
% http://ba-ngu.vo-au.com/vo/SVV08_OSPA.pdf
%
% ---BibTeX entry
% @inproceedings{OSPA2,
% author={M. Beard and B.-T. Vo and B.-N. Vo},
% booktitle = {Proc. 21st IEEE Intl. Conf. Information Fusion},
% title={Performance Evaluation for Large-Scale Multi-Target Tracking Algorithms},
% month= {Jul},
% year={2018},
% location= {Cambridge, UK}}
%
% @ARTICLE{LST,
% author={M. Beard and B.-T. Vo and B.-N. Vo},
% title = "{A Solution for Large-scale Multi-object Tracking}",
% journal = {ArXiv e-prints},
% archivePrefix = "arXiv",
% eprint = {1804.06622},
% year = 2018,
% month = apr}
%
% @ARTICLE{OSPA,
% author={D. Schuhmacher and B.-T. Vo and B.-N. Vo},
% journal={IEEE Transactions on Signal Processing},
% title={A Consistent Metric for Performance Evaluation of Multi-Object Filters},
% year={2008},
% month={Aug},
% volume={56},
% number={8},
% pages={3447-3457}}  
% ---
%
% Inputs:
% X    DxTxN array, where D is the state dimension, T is the number of
%      time steps, and N is the number of objects. NaN is used to
%      indicate when an object is not present.
%
% Y    DxTxM array, where D and T must match the dimensions of X, and M
%      is the number of tracks.
%
% c    Cutoff parameter (used for both inner and outer OSPA)
%
% p    Order parameter (used for both inner and outer OSPA)
%
% wl   Size of the moving window. For example, at time t, the metric
%      will be computed over the window [t-wl+1:t].
%
% Output:
% result   3xT array, where result(1,t) is the OSPA(2) at time t,
%          result(2,t) is the localisation component at time t, and
%          result(3,t) is the cardinality component at time t.
%
%
  
  if (size(X,1) ~= size(Y,1)) || (size(X,2) ~= size(Y,2))
    error('Dimensions of X and Y are inconsistent');
  end

  if ~isnumeric(c) || ~isscalar(c) || (c <= 0)
    error('c must be a positive scalar');
  end

  if ~isnumeric(p) || ~isscalar(p) || (p <= 0)
    error('p must be a positive scalar');
  end

  wl = floor(wl);
  if ~isnumeric(wl) || ~isscalar(wl) || (wl <= 0)
    error('Window length must be a positive integer');
  end

  eval_idx = 1:size(X,2);
  truncated = true; 
  win_off = (-wl+1):0;
  
  num_x = size(X,3);
  num_y = size(Y,3);
  num_step = size(X,2);
  num_eval = length(eval_idx);
  
  result = zeros(3,num_eval);
  
  % First, for each time index compute the matrix of inter-point
  % distances and the track existence flags
  
  distances = zeros(num_x,num_y,num_step);
  x_exists = false(num_x,num_step);
  y_exists = false(num_y,num_step);
  
  for i = 1:num_step
      
    % Compute distance between every pair of points
    x = permute(X(:,i,:),[1 3 2]);
    y = permute(Y(:,i,:),[1 3 2]);
    d = permute(sum(abs(bsxfun(@minus,permute(y,[1 3 2]),x)).^p,1),[2 3 1]);
    
    % Compute track existence flags
    x_exists(:,i) = ~isnan(x(1,:));
    y_exists(:,i) = ~isnan(y(1,:));
        
    % Distance between an empty and non-empty state
    one_exists = bsxfun(@xor,x_exists(:,i),y_exists(:,i)');
    d(one_exists) = c^p;     
    
    % Find times when neither object exists
    neither_exists = bsxfun(@and,~x_exists(:,i),~y_exists(:,i)');
    if truncated
      % Truncated window, exclude times when neither objects exists
      d(neither_exists) = NaN;
    else
      % Full window, distance between empty states is zero
      d(neither_exists) = 0;
    end
    
    % Store the distance matrix for this step
    distances(:,:,i) = d;
    
  end
  
  % Cap all inter-point distances at c^p
  if truncated
    % Truncated window
    distances = min(c^p,distances,'includenan');
  else
    % Full window
    distances = min(c^p,distances);
  end
  
  % Compute the OSPA(2) at each evaluation point
  for i = 1:num_eval
    
    % Window indices
    win_idx = eval_idx(i) + win_off;
    idx_val = (win_idx > 0) & (win_idx <= num_step);
    win_idx = win_idx(idx_val);
    
    % Compute the matrix of weighted time-averaged
    % OSPA distances between tracks
    trk_dist = mean(distances(:,:,win_idx),3,'omitnan');
    trk_dist(isnan(trk_dist)) = 0;
    
    % Get the number of objects in X and Y that exist
    % at any time inside the current window
    valid_rows = any(x_exists(:,win_idx),2);
    valid_cols = any(y_exists(:,win_idx),2);
    m = sum(valid_rows);
    n = sum(valid_cols);
    
    % Solve the optimal assignment problem
    trk_dist = trk_dist(valid_rows,valid_cols);
    if isempty(trk_dist)
        cost = 0;
    else
      if m > n
          trk_dist = trk_dist';
      end
      [~,cost] = lapjv(trk_dist);
    end
    
    % Compute the OSPA(2) distance for the current time index
    if max(m,n) == 0
      result(:,i) = 0;
    else
      result(1,i) = ( ( c^p * abs(m-n) + cost ) / max(m,n) ) ^ (1/p);
      result(2,i) = ( cost / max(m,n) ) ^ (1/p);
      result(3,i) = ( c^p * abs(m-n) / max(m,n) ) ^ (1/p);
    end
    
  end
  
end

function l_str = convert_object_label_to_string(l)
    if isempty(l)
        l_str = cell(0,1);
    else
        n_l = size(l,2);
        l_str = cell(1,n_l);
        for i = 1 : n_l
            l_str{i} = [num2str(l(1,i)),'-',num2str(l(2,i))];
        end
    end
end