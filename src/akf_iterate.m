function [x_hat, P_hat, er] = akf_iterate(x, P, Q, R, y, H, f_)
er = 0;
%% prediction step
[x_hat, P_hat] = f_(x, P); % [x_hat, P_hat] = f_(x, zeros(size(P0)));
P_hat = P_hat + Q;
%% update step
    % Update step
    K = P_hat*H' / ((H*P_hat*H' + R)); 
    x_hat = x_hat + K*(y(:)-H*x_hat);
    P_hat = (eye(length(x_hat))-K*H)*P_hat; 

    % Force symmetry on P_hat
    P_hat = (P_hat + P_hat')/2;
    % Check eigenvalues
    [~,flag] = chol(P_hat);
    if flag
        % If any is negative, find the nearest Semipositve Definite matrix
        try
            [P_hat, k]= find_nearest_spd(P_hat);
        catch ME
            if strcmp('MATLAB:svd:matrixWithNaNInf', ME.identifier)
                disp(['Error during iteration: ?']);% num2str(n)]);
            end
            rethrow(ME);
        end
        if k == -1
            % Infinite loop in the nearestSPD script. No SPD matrix found
            x_hat = zeros(size(x));
            er = 1;
            return
        end
    end
end % End function