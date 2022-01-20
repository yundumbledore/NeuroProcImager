function [x_hat, er] = analytic_kalman_filter(y,f_,nmm,H,Q,R)
    m0 = nmm.x0; % initial model state
    P0 = nmm.P0; % initial state covariance
    
    NSamples = length(y); % Number of samples. Should be equal to Ns, but setting it to lenght of the observed EEG
    NStates = length(m0); % Number of states

    % Set inital state and covariance for AKF
    x_hat = zeros(NStates, NSamples);
    P_hat = zeros(NStates, NStates, NSamples);
    x_hat(:,1) = m0;
    P_hat(:,:,1) = P0;
    
    % run AKF to iteratively estimate neurophysiological variables
    for n = 1:NSamples-1
        [x_hat(:,n+1), P_hat(:,:,n+1), er] = akf_iterate(x_hat(:,n), P_hat(:,:,n), Q, R, y(n+1), H, f_);

        if er
            x_hat = zeros(NStates, NSamples);
            return
        end
    end
end