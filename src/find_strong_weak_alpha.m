function [upper_index, lower_index, downsampled_alpha_pw] = find_strong_weak_alpha(task)
    size_cut = task.size_cut; % data truncate size
    sub_ind = task.sub_ind; % subject index
    fs = task.fs; % data sampling rate
    FS = task.FS; % saved estimates sampling rate;
    
    roi_info = './files/roi_info.mat';
    load(roi_info, 'occipital','selected_roi_index'); % load roi tissue index
    
    virtual_MEG = ['./data/data_' num2str(sub_ind) '.mat'];
    load(virtual_MEG, 'virtualdata_timeseries'); % load MEG data file
    
    % scale the signal to the same amplitude of the NMM output
    % this corresponds to data pre-processing in AKF estimation
    load('./files/model_data_std.mat','model_std','grand_max_std') % load saved standard deviation of the signal
    y = (model_std/grand_max_std)*virtualdata_timeseries; % scale the signal
    y = y(:,1:size_cut);

    downsampled_y = resample(y', FS, fs)'; % downsample the signal to match the sampling rate of the saved variable estimates
    alpha_signal = bandpass(downsampled_y, [8,12], FS); % extract alpha component from the signal
    
    [yupper,~] = envelope(alpha_signal,FS/3,'rms'); % alpha power envelope
    downsampled_alpha_pw = yupper;
    save(['./demo_cases/contrast_imaging/downsampled_alpha_pw_' num2str(sub_ind) '.mat'], 'downsampled_alpha_pw','-v7.3');
   
    % aggregate alpha power from all occipital source points
    roi_index = unique(cell2mat(selected_roi_index(1,:)));
    average_pw_list = [];
    for i = roi_index 
        indice = cell2mat(selected_roi_index(1,:))==i;
        k = selected_roi_index(3,indice);
        average_pw = mean(downsampled_alpha_pw(cell2mat(k),:),1);
        average_pw_list(i,:) = average_pw;
    end
    
    % compute occipital power over time
    occipital_alpha_pw = mean(average_pw_list(occipital,:),1); 
    
    % get strong and weak occipital alpha time indices
    lower = prctile(occipital_alpha_pw,25); % 25th percentile; power smaller than 25th regarded as weak alpha
    upper = prctile(occipital_alpha_pw,75); % 75th percentile; power greater than 75th regarded as strong alpha
    occipital_alpha_pw(2,:)=occipital_alpha_pw(1,:)>=upper; % the 2nd row in occipital_alpha_pw denotes whether the pw is greater than 75th percentile
    occipital_alpha_pw(3,:)=occipital_alpha_pw(1,:)<=lower; % the 3rd row in occipital_alpha_pw denotes whether the pw is smaller than 25th percentile
    upper_index = find(occipital_alpha_pw(2,:)==1); % get the time indices of strong alpha
    lower_index = find(occipital_alpha_pw(3,:)==1); % get the time indices of weak alpha
    
    save(['./demo_cases/contrast_imaging/upper_index_' num2str(sub_ind) '.mat'], 'upper_index','-v7.3');
    save(['./demo_cases/contrast_imaging/lower_index_' num2str(sub_ind) '.mat'], 'lower_index','-v7.3');
end