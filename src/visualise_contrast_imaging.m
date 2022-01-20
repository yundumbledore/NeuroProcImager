function [] = visualise_contrast_imaging(task)
    name_list = task.name_list;
    sub_ind = task.sub_ind;
    ncpu = task.ncpu;
    cmc = task.cmc;
    
    % only need 9 cpus
    if ncpu > 9
        ncpu = 9;
    end
    
    % for each neurophysiological variable perform corrections for multiple
    % comparison and visualise on MRI image
    poolobj = parpool(ncpu);
    parfor (par_ind = 1:9, ncpu)
        parameter_name = name_list{par_ind};
        stats_correction(parameter_name, sub_ind, cmc);
    end
    delete(poolobj)
end

function [] = stats_correction(parameter_name, sub_ind, cmc)
    load(['./demo_cases/contrast_imaging/permuted_sub_' num2str(sub_ind) '_tstats_' parameter_name '.mat'],'eavg_across');
    tstat = eavg_across.t_stat_original';
    
    pmax = prctile(eavg_across.mx,(100-0.5)); % critical upper threshold
    pmin = prctile(eavg_across.mn,(0.5)); % critical lower threshold
    
    output_dir = './demo_cases/contrast_imaging/visualisation_outputs/';
    % output uncorrected imaging
    corrected = 0;
    evalc('visualise(sub_ind, tstat, corrected, parameter_name, output_dir)');
    
    % output uncorrected imaging
    if cmc
        corrected = 1;
        ind_voxel = find((tstat(:,1) <= pmax) & (tstat(:,1) >= pmin)); % find insignificant voxels (within the critical thresholds)

        if length(ind_voxel) ~= 4714
            tstat(ind_voxel) = 0; % set insignificant voxels to zero
            evalc('visualise(sub_ind, tstat, corrected, parameter_name, output_dir)');
        end
    end
end