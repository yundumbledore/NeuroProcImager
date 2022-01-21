function [] = main(tasks)
    %% create working directories if not exist
    contrast_image_dir = './demo_cases/contrast_imaging/visualisation_outputs'; % directory to put contrast imaging
    output_dir = './output'; % directory to put variable estimates

    if ~exist(contrast_image_dir, 'dir')
        mkdir(contrast_image_dir)
    end
    
    if ~exist(output_dir, 'dir')
        mkdir(output_dir)
    end
    
    %% define a dictionary to put variable names
    name_list{1} = 'aEP';
    name_list{2} = 'aIP';
    name_list{3} = 'aPE';
    name_list{4} = 'aPI';
    name_list{5} = 'input';
    name_list{6} = 'v_es';
    name_list{7} = 'v_ii';
    name_list{8} = 'v_pyr';
    name_list{9} = 'alpha_pw';
    
    %% run showcases below
    pipeline = tasks.pipeline; % get which showcase to run
    
    listing = dir('./data'); % check data folder
    name = {listing.name}; % get the data file name
    name = name(~strncmp(name, '.', 1)); % get valid file names by excluding '.' and '..'
    
    if ismember("AKF estimation", pipeline) % run AKF model inversion showcase
        disp('Now running the 1st demo case: AKF estimation...')
        for i = 1:length(name) % iterate over all data files if there is more than one to run
            tmp = regexp(name{i}, '[_.]', 'split');
            sub_ind = tmp{2}; % get subject index
            
            task.sub_ind = str2num(sub_ind);
            task.ncpu = tasks.ncpu;
            task.fs = tasks.fs;
            task.Fs = tasks.Fs;
            task.FS = tasks.FS;
            task.size_cut = tasks.size_cut;
            [N_channels, N_timesteps] = run_akf(task); % AKF estimation starts here
            
            task.N_channels = N_channels;
            task.N_timesteps = N_timesteps;
            fetch_estimates(task) % fetch estimates from individual files and save into a single .mat file
        end
    end  
    
    if ismember("Contrast imaging", pipeline) % run contrast imaging showcase
        disp('Now running the 2nd demo case: Contrast imaging...')
        cmc = tasks.multiple_comparison_correction; % whether to perform corrections for multiple comparisons
        
        for i = 1:length(name) % iterate over all data files if there is more than one to run
            tmp = regexp(name{i}, '[_.]', 'split');
            sub_ind = tmp{2}; % get subject index
            
            task.sub_ind = str2num(sub_ind);
            task.size_cut = tasks.size_cut;
            task.fs = tasks.fs;
            task.FS = tasks.FS;
             
            [upper_index, lower_index, downsampled_alpha_pw] = find_strong_weak_alpha(task); % get alpha power from the signal and get occipital strong and weak alpha time indices
          
            load(['./output/variable_estimates_' num2str(sub_ind) '.mat'],'compact_list') % load variable estimates that have been obtained in showcase 1
            compact_list{9} = downsampled_alpha_pw; % add it into compact_list for parfor loop
            clear downsampled_alpha_pw
           
            % Run two sample t test to compare sample mean of variable estimates and alpha power between strong and weak occipital alpha power
            poolobj = parpool(tasks.ncpu);
            parfor (par_ind = 1:9, tasks.ncpu)
                upper_sample = compact_list{par_ind}(:,upper_index); % retrieve variable estimates with respect to strong occipital alpha
                lower_sample = compact_list{par_ind}(:,lower_index); % retrieve variable estimates with respect to weak occipital alpha
                parameter_name = name_list{par_ind}; % get the variable name
                twosample_ttest(str2num(sub_ind), parameter_name, upper_sample, lower_sample, cmc);
            end
            delete(poolobj)
 
            % whole brain visulisation starts here
            task.name_list = name_list;
            task.ncpu = tasks.ncpu;
            task.cmc = cmc;
            visualise_contrast_imaging(task); 
        end
    end
end
