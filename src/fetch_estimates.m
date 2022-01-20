function [] = fetch_estimates(task)
    %% initial some variables
    sub_ind = task.sub_ind;
    N_channels = task.N_channels;
    N_timesteps = task.N_timesteps;
    
    chan_error = zeros(N_channels,1);
    aEP = zeros(N_channels, N_timesteps);
    aIP = zeros(N_channels, N_timesteps);
    aPE = zeros(N_channels, N_timesteps);
    aPI = zeros(N_channels, N_timesteps);
    input = zeros(N_channels, N_timesteps);
    v_es = zeros(N_channels, N_timesteps);
    v_ii = zeros(N_channels, N_timesteps);
    v_pyr = zeros(N_channels, N_timesteps);
    
    %% go to each file to fetch estimates
    for iCh = 1:N_channels
        try
            load_file = ['./output/' num2str(sub_ind) '/est_Ch_' num2str(iCh) '.mat'];
            load(load_file)
            aEP(iCh,:) = aEP_hat;
            aIP(iCh,:) = aIP_hat;
            aPE(iCh,:) = aPE_hat;
            aPI(iCh,:) = aPI_hat;
            input(iCh,:) = input_hat;
            v_es(iCh,:) = v_es_hat;
            v_ii(iCh,:) = v_ii_hat;
            v_pyr(iCh,:) = v_pyr_hat; 
        catch me
            chan_error(iCh,1) = 1;
        end
    end
    
    %% put estimates into a dictionary
    compact_list{1} = aEP;
    compact_list{2} = aIP;
    compact_list{3} = aPE;
    compact_list{4} = aPI;
    compact_list{5} = input;
    compact_list{6} = v_es;
    compact_list{7} = v_ii;
    compact_list{8} = v_pyr;
    
    save(['./output/variable_estimates_' num2str(sub_ind) '.mat'],'compact_list', '-v7.3');
    save(['./output/variable_estimates_errorReport_' num2str(sub_ind) '.mat'],'chan_error','-v7.3');
end