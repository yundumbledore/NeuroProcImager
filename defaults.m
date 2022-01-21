%% create some directories if they do not exist
% -------------------------------------------------------------------------
% After downloading the package, you should execute the defaults function 
% (formerly called defaults.m), which sets the defaults and  configures up 
% the minimal required path settings.
% -------------------------------------------------------------------------
%%
contrast_image_dir = './demo_cases/contrast_imaging/visualisation_outputs'; % directory to save contrast imaging
output_dir = './output'; % directory to save variable estimates
data_dir = './data'; % directory to put data

if ~exist(contrast_image_dir, 'dir')
    mkdir(contrast_image_dir)
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir)
end

if ~exist(data_dir, 'dir')
    mkdir(data_dir)
end