function [N_channels, N_timesteps] = run_akf(task)
    %% hyperparameters
    ncpu = task.ncpu;
    fs = task.fs;
    Fs = task.Fs;
    FS = task.FS;
    size_cut = task.size_cut;
    sub_ind = task.sub_ind;

    %% make output directory if not exists
    saving_dir = ['./output/' num2str(sub_ind)];
    if ~exist(saving_dir, 'dir')
        mkdir(saving_dir)
    end
    
    %% load data
    load(['./data/data_' num2str(sub_ind) '.mat'],'virtualdata_timeseries'); % change this path to load alternative data
    virtualdata_timeseries = virtualdata_timeseries(:,1:size_cut);
    [N_channels, N_timesteps]  = size(virtualdata_timeseries);
    N_timesteps = N_timesteps*(FS/fs);

    %% set filter
    [A,B,C,N_states,N_syn,N_inputs,N_samples,xi,v0,varsigma,Q,R,H,s_y] = set_params(fs); % set parameters for the filter

    model_std = std(s_y); % find model output standard deviation
    grand_max_std = max(std(virtualdata_timeseries, [], 2)); % find the maximum standard deviation across all source points of the subject
    save('./files/model_data_std.mat','model_std','grand_max_std','-v7.3')
    
    y = (model_std/grand_max_std)*virtualdata_timeseries; % scale data to match the amplitude of the NMM output
    y = resample(y',Fs,fs)'; % upsample data to make smoother estimation/downsample data to save time
    clear virtualdata_timeseries

    % set model matrices
    nmm.A = A;
    nmm.B = B;
    nmm.C = C;

    f_ = @(x,P)nmm_run(nmm, x, P); % this is the analytic time update equation

    % set filter initial state and covariance
    m0 = mean(xi(:,N_samples/2:end),2);
    nmm.x0 = m0;
    P_hat_init = 10*cov(xi(:,N_samples/2:end)');
    P_hat_init(2*N_syn+1:end,2*N_syn+1:end) = eye(N_syn+N_inputs)*10e-2;
    nmm.P0 = P_hat_init;

    %% apply filter
    poolobj = parpool(ncpu);
    parfor (ich = 1:4714, ncpu)
        y_ich = y(ich,:);
        [xi_hat, er] = analytic_kalman_filter(y_ich,f_,nmm,H,Q,R);
        if ~er % save estimates if no error
            save_estimates(xi_hat, sub_ind, ich, FS, Fs)
        end
    end
    delete(poolobj)
end

function [] = save_estimates(xi_hat, sub_ind, ich, FS, Fs)
    v_pyr_hat = xi_hat(1,:) + xi_hat(7,:) + xi_hat(9,:); % pyramidal membrane potential
    v_pyr_hat_ = v_pyr_hat/50;

    v_es_hat = xi_hat(5,:); % excitatory stellate membrane potential
    v_ii_hat = xi_hat(3,:); % inhibitory interneuron membrane potential
    theta_hat = xi_hat(9:end,:);

    input_hat_ = theta_hat(1,:);
    aIP_hat_ = theta_hat(2,:);
    aPI_hat_ = theta_hat(3,:);
    aPE_hat_ = theta_hat(4,:);
    aEP_hat_ = theta_hat(5,:);

    % downsample to the specified sampling rate to save space
    v_pyr_hat = single(resample(v_pyr_hat_,FS,Fs));
    v_es_hat = single(resample(v_es_hat,FS,Fs));
    v_ii_hat = single(resample(v_ii_hat,FS,Fs));
    input_hat = single(resample(input_hat_,FS,Fs));
    aIP_hat = single(resample(aIP_hat_,FS,Fs));
    aPI_hat = single(resample(aPI_hat_,FS,Fs));
    aPE_hat = single(resample(aPE_hat_,FS,Fs));
    aEP_hat = single(resample(aEP_hat_,FS,Fs));
    save(['./output/' num2str(sub_ind) '/est_Ch_' num2str(ich) '.mat'],...
        'v_pyr_hat','v_es_hat','v_ii_hat','input_hat','aIP_hat','aPI_hat','aPE_hat','aEP_hat','-v7.3');
end