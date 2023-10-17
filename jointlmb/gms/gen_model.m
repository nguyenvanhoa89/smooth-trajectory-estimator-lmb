function model= gen_model

% basic parameters
model.x_dim= 4;   %dimension of state vector
model.z_dim= 2;   %dimension of observation vector

model.K = 100; % Total number of time steps.

% dynamical model parameters (CV model)
model.T= 1;                                     %sampling period
model.A0= [ 1 model.T; 0 1 ];                         %transition matrix                     
model.F= [ model.A0 zeros(2,2); zeros(2,2) model.A0 ];
model.B0= [ (model.T^2)/2; model.T ];
model.B= [ model.B0 zeros(2,1); zeros(2,1) model.B0 ];
model.sigma_v = 5;
model.Q= (model.sigma_v)^2* model.B*model.B';   %process noise covariance

% survival/death parameters
model.P_S= .99;
model.Q_S= 1-model.P_S;

% birth parameters (LMB birth model, single component only)
model.T_birth= 4;         %no. of LMB birth terms
model.L_birth= zeros(model.T_birth,1);                                          %no of Gaussians in each LMB birth term
model.r_birth= zeros(model.T_birth,1);                                          %prob of birth for each LMB birth term
model.w_birth= cell(model.T_birth,1);                                           %weights of GM for each LMB birth term
model.m_birth= cell(model.T_birth,1);                                           %means of GM for each LMB birth term
model.B_birth= cell(model.T_birth,1);                                           %std of GM for each LMB birth term
model.P_birth= cell(model.T_birth,1);                                           %cov of GM for each LMB birth term

model.L_birth(1)=1;                                                             %no of Gaussians in birth term 1
model.r_birth(1)=0.05;                                                          %prob of birth
model.w_birth{1}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
model.m_birth{1}(:,1)= [ 0.1; 0; 0.1; 0 ];                                      %mean of Gaussians
model.B_birth{1}(:,:,1)= diag([ 10; 10; 10; 10 ]);                              %std of Gaussians
model.P_birth{1}(:,:,1)= model.B_birth{1}(:,:,1)*model.B_birth{1}(:,:,1)';      %cov of Gaussians

model.L_birth(2)=1;                                                             %no of Gaussians in birth term 2
model.r_birth(2)=0.05;                                                          %prob of birth
model.w_birth{2}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
model.m_birth{2}(:,1)= [ 400; 0; -600; 0 ];                                     %mean of Gaussians
model.B_birth{2}(:,:,1)= diag([ 10; 10; 10; 10 ]);                              %std of Gaussians
model.P_birth{2}(:,:,1)= model.B_birth{1}(:,:,1)*model.B_birth{1}(:,:,1)';      %cov of Gaussians

model.L_birth(3)=1;                                                             %no of Gaussians in birth term 3
model.r_birth(3)=0.05;                                                          %prob of birth
model.w_birth{3}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
model.m_birth{3}(:,1)= [ -800; 0; -200; 0 ];                                    %mean of Gaussians
model.B_birth{3}(:,:,1)= diag([ 10; 10; 10; 10 ]);                              %std of Gaussians
model.P_birth{3}(:,:,1)= model.B_birth{1}(:,:,1)*model.B_birth{1}(:,:,1)';      %cov of Gaussians

model.L_birth(4)=1;                                                             %no of Gaussians in birth term 4
model.r_birth(4)=0.05;                                                          %prob of birth
model.w_birth{4}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
model.m_birth{4}(:,1)= [ -200; 0; 800; 0 ];                                     %mean of Gaussians
model.B_birth{4}(:,:,1)= diag([ 10; 10; 10; 10 ]);                              %std of Gaussians
model.P_birth{4}(:,:,1)= model.B_birth{1}(:,:,1)*model.B_birth{1}(:,:,1)';      %cov of Gaussians

% observation model parameters (noisy x/y only)
model.H= [ 1 0 0 0 ; 0 0 1 0 ];    %observation matrix
model.D= diag([ 10; 10 ]); 
model.R= model.D*model.D';              %observation noise covariance

% detection parameters
model.P_D= .88;   %probability of detection in measurements
model.Q_D= 1-model.P_D; %probability of missed detection in measurements

% clutter parameters
model.lambda_c= 66;                             %poisson average rate of uniform clutter (per scan)
model.range_c= [ -1000 1000; -1000 1000 ];      %uniform clutter region
model.limit=model.range_c;
model.pdf_c= 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density

% Single-object filter type
model.source_info.meas_type = 'pos';
model.source_info.source_pos = [0,0];
model.filter.type = 'kf';
switch model.filter.type
    case 'kf'
        model.filter.predict_func = @kalman_predict_multiple;
        model.filter.gate_func    = @gate_meas_gms_idx;
        model.filter.update_func  = @kalman_update_multiple; 
        model.filter.smooth_func  = @rts_smooth_estimate; 
    case 'ekf'
        model.filter.predict_func = @ekf_predict_multiple;
        model.filter.gate_func    = @gate_meas_ekf_idx;
        model.filter.update_func  = @ekf_update_multiple; 
    case 'ukf'
        model.filter.predict_func = @ukf_predict_multiple;
        model.filter.gate_func    = @gate_meas_ukf_idx;
        model.filter.update_func  = @ukf_update_multiple; 
        model.filter.smooth_func  = @urts_smooth_estimate; 
        model.filter.ukf_alpha= 1;                %scale parameter for UKF - choose alpha=1 ensuring lambda=beta and offset of first cov weight is beta for numerical stability
        model.filter.ukf_beta= 2;                 %scale parameter for UKF
        model.filter.ukf_kappa= 2; %2;                %scale parameter for UKF (alpha=1 preferred for stability, giving lambda=1, offset of beta for first cov weight)
        kappa = model.filter.ukf_kappa; nx = model.x_dim;
        model.filter.ukf_W=[kappa/(kappa+nx) 0.5/(kappa+nx)+zeros(1,2*nx)];
end
model.ospa.c = 100;
model.ospa.p = 1; 
model.ospa.win_len = 10;
model.abp.enable = false;
model.file_name = [model.filter.type,'_pD_',num2str(model.P_D),...
                '_ldc_',num2str(model.lambda_c),...
                '.mat'];
