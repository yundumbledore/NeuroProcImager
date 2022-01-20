%% Nonparametric permutation test for t-test. 
%% Note that this code should be used for an individual imaging. This is not for group-level imaging.
% Nichols, T.E. & Holmes, A.P. Nonparametric permutation tests for functional neuroimaging: a primer with examples. Human brain mapping 15, 1-25 (2002).

function [eavg_across] = np_pm_test_ttest(x_, y_, nperm)
    eavg_across.mx = zeros(1,nperm);
    eavg_across.mn = eavg_across.mx;
    sc = 1;%sqrt(size(x_,1));
    data_length = size(x_, 1);
    
    for i = 1:nperm
        if i == 1
            [~,~,~,stats] = ttest2(x_,y_); % calculate t-score with unpermuted data
            tscore = stats.tstat;
            eavg_across.t_stat_original = tscore*sc; % remember the unpermuted result
        else
            x = x_; % make a copy
            y = y_;
            tmp = rand(data_length,1); 
            ind = tmp >= 0.5;
            idx = find(ind == 1); % row index to swap elements between x and y

            y_copy = y; % swap here
            y(idx,:) = x(idx,:);
            x(idx,:) = y_copy(idx,:);
        
            [~,~,~,stats] = ttest2(x,y); % calculate t-score with permuted data
            tscore = stats.tstat;
        end
        
        % find the max and min t-score from the whole brain source points
        eavg_across.mx(i) = max(tscore)*sc;
        eavg_across.mn(i) = min(tscore)*sc;
    end
end