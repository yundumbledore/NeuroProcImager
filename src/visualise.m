function [] = visualise(sub_ind, tstat, corrected, parameter_name, output_dir)
    load('./files/template_6mmROI.mat', 'template_sourcemodel'); % load template source model
    load('./files/mri_temp.mat', 'mri_temp'); % load normalised segmented mri
    load("./files/all_roi_points.mat", 'all_roi_points'); % load upsampled roi points coordinates
    load('./files/cerebellum_MNI_pos.mat', 'cerebellum_MNI_pos'); % load cerebellum source points coordinates
    
    % prepare functional data (statistics) for both cerebral and cerebellum
    all_source_points_activation = [tstat(1:2067); zeros(70,1); tstat(2068:2153); zeros(19,1);tstat(2154:3998); zeros(258,1); tstat(3999:4714)];
    activation = all_source_points_activation;
    ROI_pos_MNI = template_sourcemodel.pos;
    ROI_pos_MNI = [ROI_pos_MNI; cerebellum_MNI_pos];
    activation = [activation; zeros(length(cerebellum_MNI_pos),1)];
    
    % upsample to make high resolution image
    vq2 = griddatan(ROI_pos_MNI,activation,all_roi_points);
    ROI_pos_MNI = all_roi_points;
    activation = vq2;
    
    % make a variable to include source points coordinates and activation
    alpha_stats.pos = ROI_pos_MNI;
    alpha_stats.pow = activation;
    alpha_stats.inside = logical(ones(length(activation),1));
    ind_nan = find(isnan(alpha_stats.pow(:)));
    alpha_stats.pow(ind_nan) = 0;

    % interpolate the functional data to the anatomical data
    cfg = [];
    cfg.parameter = 'pow';
    cfg.interpmethod  = 'nearest';
    alpha_stats_int = ft_sourceinterpolate(cfg,alpha_stats,mri_temp);
    ind_nan = find(isnan(alpha_stats_int.pow(:)));
    alpha_stats_int.pow(ind_nan) = 0;
    
    % generate images (D1 lateral view, D2 rear view, D3 dorsal view)
    if ~corrected % generate slice images for corrected statistics
        cfg = [];
        cfg.method = 'slice';
        cfg.slicerange = [30 150]; 
        cfg.slicedim = 1;
        cfg.nslices = 20;
        cfg.funparameter = 'pow';
        cfg.maskparameter  = cfg.funparameter;
        ft_sourceplot(cfg, alpha_stats_int);
        saveas(gcf,[output_dir 'sub_' num2str(sub_ind) '_D1_' parameter_name '_raw.jpg'])
        close all
        
        cfg = [];
        cfg.method = 'slice';
        cfg.slicerange = [30 150]; 
        cfg.slicedim = 2;
        cfg.nslices = 20;
        cfg.funparameter = 'pow';
        cfg.maskparameter  = cfg.funparameter;
        ft_sourceplot(cfg, alpha_stats_int);
        saveas(gcf,[output_dir 'sub_' num2str(sub_ind) '_D2_' parameter_name '_raw.jpg'])
        close all
        
        cfg = [];
        cfg.method = 'slice';
        cfg.slicerange = [30 150]; 
        cfg.slicedim = 3;
        cfg.nslices = 20;
        cfg.funparameter = 'pow';
        cfg.maskparameter  = cfg.funparameter;
        ft_sourceplot(cfg, alpha_stats_int);
        saveas(gcf,[output_dir 'sub_' num2str(sub_ind) '_D3_' parameter_name '_raw.jpg'])
        close all
    else % generate slice images for uncorrected statistics
        cfg = [];
        cfg.method = 'slice';
        cfg.slicerange = [30 150]; 
        cfg.slicedim = 1;
        cfg.nslices = 20;
        cfg.funparameter = 'pow';
        cfg.maskparameter  = cfg.funparameter;
        ft_sourceplot(cfg, alpha_stats_int);
        saveas(gcf,[output_dir 'sub_' num2str(sub_ind) '_D1_' parameter_name '_corrected.jpg'])
        close all
        
        cfg = [];
        cfg.method = 'slice';
        cfg.slicerange = [30 150]; 
        cfg.slicedim = 2;
        cfg.nslices = 20;
        cfg.funparameter = 'pow';
        cfg.maskparameter  = cfg.funparameter;
        ft_sourceplot(cfg, alpha_stats_int);
        saveas(gcf,[output_dir 'sub_' num2str(sub_ind) '_D2_' parameter_name '_corrected.jpg'])
        close all
        
        cfg = [];
        cfg.method = 'slice';
        cfg.slicerange = [30 150]; 
        cfg.slicedim = 3;
        cfg.nslices = 20;
        cfg.funparameter = 'pow';
        cfg.maskparameter  = cfg.funparameter;
        ft_sourceplot(cfg, alpha_stats_int);
        saveas(gcf,[output_dir 'sub_' num2str(sub_ind) '_D3_' parameter_name '_corrected.jpg'])
        close all
    end
end