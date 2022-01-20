function varargout = nmm_run(nmm, x, P)
    % Model
    A       = nmm.A;
    B       = nmm.B;
    C       = nmm.C;

    % Indexes
    v_idx       = [1 3 5 7];
    z_idx       = [2 4 6 8];
    u_idx       = 9;
    alpha_idx   = [10 11 12 13];

    NStates = length([v_idx z_idx u_idx alpha_idx]);
    NSynapses = length(v_idx);

    % The parameters
    r           = 3.033;         % varsigma
    v0          = 6;        % Threshold

    % The numeric solution of P6
    w_fast = [0.1713244923791705 0.3607615730481384  0.4679139345726904 0.1713244923791705, 0.3607615730481384 0.4679139345726904];
    y_fast = [0.033765242898424 0.169395306766868 0.380690406958402 0.966234757101576, 0.830604693233132 0.619309593041599];

    C_inhibit = C; % Auxiliary matrix C_inhibit, same as C but the inhibitory element is negative to include it in the nonlinearity as inhibition

    CPC = C(z_idx,:)*P*C(z_idx,:)';
    dCPB = diag(C(z_idx,:)*P*B(z_idx,:)');
    Bxi = B(z_idx,alpha_idx)*x(alpha_idx);

    gamma = 1./sqrt(2*(diag(CPC) + r^2));
    beta = (C_inhibit(z_idx,[v_idx u_idx])*x([v_idx u_idx]) - v0).*gamma;       % <- v_0
    Xi = (erf(beta) + 1)/2;
    Upsilon = exp(-(beta.^2)).*gamma./sqrt(pi);
    psi = Bxi.*Xi + dCPB.*Upsilon;                              % E[ Bxi_t o g(Cxi_t) ]

    %% Analytic mean
    analytic_mean = A*x;
    analytic_mean(2:2:2*NSynapses) = analytic_mean(2:2:2*NSynapses) + psi;

    %% Analytic Covariance (Pip)
    % cov part 1 (Phi)
    q2 = Upsilon.*(Bxi - dCPB.*beta.*gamma.^2 *2);
    AP = A*P;        
    Phi = ones(NStates,1)*q2'.*(AP*C_inhibit(z_idx,:)') + ones(NStates,1)*Xi'.*(AP*B(z_idx,:)');

    CPCgammagamma = asin(CPC.*(gamma*gamma')*2);
    CPCgammagamma = CPCgammagamma(:);                                           % change to a vector
    CPCgammagamma_y = CPCgammagamma*y_fast;                                     % a matrix of row vectors

    betabeta = (beta*beta')*2;
    betabeta_mat = betabeta(:)*ones(1,length(w_fast));                          % change to a vector and repmat

    beta2mat = (beta.^2)*ones(1,NSynapses);                                         % sq and repmat
    beta2matT = beta2mat';                                                      % the transpose allow for permutations when we sum below
    beta2_plus_beta2T = (beta2mat(:) + beta2matT(:))*ones(1,length(w_fast));    % change to a vectors, add, and repmat

    % put it together
    %
    Psi = reshape(sum(CPCgammagamma*w_fast.*exp(-(beta2_plus_beta2T - betabeta_mat.*sin(CPCgammagamma_y))./cos(CPCgammagamma_y).^2),2),NSynapses,NSynapses)/(4*pi);
    %                 ~~~~~~~~~~~~~~~~~~~~  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %             ~~~
    Omega = ((Xi*Xi') + Psi).*(Bxi*Bxi' + P(alpha_idx,alpha_idx));
    %        ~~~~~~~~~~~~~~    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %         E[g^2(Cxi)]           (alpha^2 + sigma_alpha^2)


    %%
    % here we construct the cov mat.
    %
    P_1m = AP*A';% + Q; % Q added outside this function
    P_1m(z_idx,z_idx) = P_1m(z_idx,z_idx)  + Omega - psi*psi';
    P_1m(:,z_idx) = P_1m(:,z_idx) + Phi;
    P_1m(z_idx,:) = P_1m(z_idx,:) + Phi';
    analytic_cov = P_1m;


    varargout{1} = analytic_mean; % E_t
    varargout{2} = analytic_cov; % P_t            
end % End function