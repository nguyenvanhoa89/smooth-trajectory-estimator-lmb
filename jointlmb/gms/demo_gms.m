%% Info
% This is the demo script for the The Smooth Trajectory Estimator for LMB Filters proposed in
% H. V. Nguyen, T. T. D. Nguyen, C. Shim, M. Anuar, "The Smooth Trajectory Estimator for LMB Filters", 
% in Proc of 2023 ICCAIS, Nov 2023, pp. 1--6.

% The implemtation was developed based on the following codes provided in the following papers

% [1]  T. T. D. Nguyen and D. Y. Kim, "GLMB Tracker with Partial Smoothing," Sensors, Vol. 19, No. 20, pp. 4419, 2019.  
% https://https://www.mdpi.com/1424-8220/19/20/4419

% [2] B.-T. Vo, and B.-N. Vo, "An Efficient Implementation of the Generalized Labeled Multi-Bernoulli Filter," 
% IEEE Trans. Signal Processing, Vol. 65, No. 8, pp. 1975-1987, 2017.
% http://ba-ngu.vo-au.com/vo/VVH_FastGLMB_TSP17.pdf

% Please cite our paper if you are using our codes as follows: 
%
% @inproceedings{nguyen2023the,
%   title={The Smooth Trajectory Estimator for LMB Filters},
%   author={Nguyen, Hoa Van and Nguyen, Tran Thien Dat and Shim, Changbeom and Anuar, Marzhar},
%   booktitle={Proceedings of the 12th International Conference on Control, Automation and Information Sciences (ICCAIS)},
%   year={2023},
%   pages={1--6},
%   organization={IEEE}}
%
%---
%% Workspace preparation
rng(10); dbstop if error
close all; clc; clearvars;
addpath(genpath('../../_common'));

%% Main program
model= gen_model;
truth= gen_truth(model);
meas=  gen_meas(model,truth);
[est,est_re]=   run_filter(model,meas);

%% Plot the results
colorarray = load('colorarray.mat');
plot_est_vs_truth(model, truth, est,'colorarray',colorarray.colorarray);set(gca,'FontSize',20); 
plot_est_vs_truth(model, truth, est_re,'colorarray',colorarray.colorarray);set(gca,'FontSize',20); 

[~,ospa_vals_lmb,ospa2_vals_lmb]= plot_results(model,truth,meas,est,'plot_flag',false,'save_plot_flag',false);
[~,ospa_vals_lmbre,ospa2_vals_lmbre]= plot_results(model,truth,meas,est_re,'plot_flag',false,'save_plot_flag',false);

temp_txt = 'LMB Results - Linear Scenario:';
disp(temp_txt);
temp_txt = sprintf('OSPA  Dist:  %s [m]',jsonencode(round(mean(ospa_vals_lmb,1),2)));
fprintf([temp_txt,'\n']);
temp_txt=sprintf('OSPA2 Dist:  %s [m]',jsonencode(round(ospa2_vals_lmb(end,:),2)));
fprintf([temp_txt,'\n']);

temp_txt = 'LMB-RE Results - Linear Scenario:';
disp(temp_txt);
temp_txt = sprintf('OSPA  Dist:  %s [m]',jsonencode(round(mean(ospa_vals_lmbre,1),2)));
fprintf([temp_txt,'\n']);
temp_txt=sprintf('OSPA2 Dist:  %s [m]',jsonencode(round(ospa2_vals_lmbre(end,:),2)));
fprintf([temp_txt,'\n']);
