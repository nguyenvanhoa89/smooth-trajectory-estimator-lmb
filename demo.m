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


%% Prepare workspace
clc; clearvars; close all;
addpath(genpath(pwd))
%% Run experiment
% demo_gms % linear dynamic scenario (uncomment to run)
demo_ukf % non-linear dynamic scenario (uncomment to run)