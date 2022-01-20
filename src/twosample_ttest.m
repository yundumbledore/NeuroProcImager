function [] = twosample_ttest(sub_ind, parameter_name, upper_sample, lower_sample, cmc)
    % if correction for multiple comparison is true, permute the data 500 times
    if cmc
        nperm = 500;
    else
        nperm = 1;
    end
    
    % run two-sample t-test with or without nonparametric permutation test
    [eavg_across] = np_pm_test_ttest(upper_sample', lower_sample', nperm); 
    
    save(['./demo_cases/contrast_imaging/permuted_sub_' num2str(sub_ind) '_tstats_' parameter_name '.mat'],'eavg_across','-v7.3');
end