clear all
clc
%%
% Author: Yun Zhao (Monash University, Australia)
% Email: yun.zhao1@monash.edu
% -------------------------------------------------------------------------
% To run this demo, one needs to install "signal processing toolbox", 
% "statistics and machine learning toolbox", "parallel computing toolbox" 
% in the MATLAB.
%
% Fieldtrip version equal to or later than 20200828 is also required and 
% available to download from https://www.fieldtriptoolbox.org.
% -------------------------------------------------------------------------
% Here we provide two showcases to enable reviewers and readers to get
% in touch with the inference-based whole-brain imaging framework.
% 
% Case 1: Neurophysiological variables inference by AKF
% To run Case 1, one needs to put "AKF estimation" in the pipeline variable
% below.
%
% Case 2: Contrast imaging between strong and weak occipital alpha 
% oscillation. 
% To run Case 2, one needs to put "Contrast imaging" in the pipeline 
% variable below.
%
% Note: Case 1 should be run prior to Case 2.
% -------------------------------------------------------------------------
% Here we show basic time statistics of running this framework using a
% MacBook Pro (Processor 2.9 GHz 6-Core Intel Core i9, Memory 32 GB 2400 
% MHz DDR4).
%
% "AKF estimation": 3.16 hr on 6 cpus (14.5 seconds each channel)
% "Contrast imaging": 1 hr on 6 cpus (with corrections for multiple comparisons)
%                     0.14 hr on 6 cpus (without corrections for multiple comparisons)
% -------------------------------------------------------------------------
%
%% add source to MATLAB path
addpath(genpath('./src'))

%% users needs to define following critical parameters
pipeline = ["Contrast imaging"];                                            % which showcase to run. two showcases: "AKF estimation", "Contrast imaging".
ncpu = feature('numcores');                                                 % number of cpus you have
multiple_comparison_correction = 1;                                         % whether to perform corrections for multiple comparison problem (0 or 1)

%% define some other parameters and save them in one variable 'tasks'
tasks.pipeline = pipeline;
tasks.ncpu = ncpu;  
tasks.multiple_comparison_correction = multiple_comparison_correction;  
tasks.fs = 400;                                                             % data sampling rate
tasks.Fs = 400;                                                             % sampling rate for estimation; increase Fs for precise tracking model states
tasks.FS = 150;                                                             % saved estimates sampling rate; reduce FS to save space
tasks.size_cut = 89000;                                                     % if there is more than one subject, truncate recordings to be the same size

tic
try 
    main(tasks)
catch ME
    disp(ME)
    
    % shut down parallel pool if it exists
    poolobj = gcp('nocreate');
    delete(poolobj)
end
toc
