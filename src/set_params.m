function [A,B,C,N_states,N_syn,N_inputs,N_samples,xi,v0,varsigma,Q,R,H,s_y] = set_params(fs)
    input = 300; % cortical input

    % set_params.m

    % set parameters - remove redundant states

    scale = 50;                         % this is to get states and derivatives on the same order of magnitude

    % set units of membrane potentials (1 for mV, 0 for V)
    %
    mV = 1;
    V2mVfactor = 1e3;

    % time parameters
    %
    Fs = fs;                           % sampling rate per second
    dt = 1/Fs;                         % time step (s)
    TimeOfSim = 50;                    % time of simulation (s)
    N_samples = round(TimeOfSim/dt);   % number of time points to simulate

    N_inputs = 1;
    N_syn = 4;                            % number of synapses
    N_states = 3*N_syn + N_inputs;        % u, z, alpha plus an input

    %% define the disturbance covariance matrix
    %
    sigma_all = 5e-8;                               % something small for all states
    sigma_input = 5e-4;                             % for input
    sigma_params = 5e-5;%sigma_all;
    sigma_offset = 5e-6;
    sigma_R = 1e-3;

    Q = eye(N_states)*(scale*sqrt(dt)*sigma_all)^2;                             % add a tiny bit of noise to all states (added for numerical stability)
    Q(2*N_syn+1:end,2*N_syn+1:end) = eye(N_syn+N_inputs)*(scale*sqrt(dt)*sigma_params)^2;
    % **** HARDCODED NUMBERS HERE
    Q(2*N_syn+1,2*N_syn+1) = (scale*sqrt(dt)*sigma_input)^2;
    if N_inputs > 1
        Q(2*N_syn+2,2*N_syn+2) = (scale*sqrt(dt)*sigma_offset)^2;
    end

    % measurement disturbance covariance
    %
    R = sigma_R^2;

    %% General parameters from J&R
    %
    % sigmoid bits
    %
    f_max = 2.5;  % maximum firing rate (spikes/s)
    r = 560; 
    varsigma = 1.699/r;     % slope of the sigmoid                                                     % (spikes/(Vs))
    varsigma_sq = varsigma^2;
    v0 = 0.006;    % mean firing threshold                                                             % (V)

    % synaptic gains
    %
    alpha_e = 3.25e-3;                                                          % gain of excitatory synapses (V)
    alpha_i = -22e-3;                                                           % gain of inhibitory synapses (V)

    % synaptic kernel time constants
    %
    ex_tau = 0.010;                     % excitatory synaptic time constant (s)
    in_tau = 0.020;                     % inhibitory synaptic time constant (s)

    % input to py population
    %
    input = input*scale;                   % ! scale here ! the input is a membrane potential input the pyramidal population this is similar to setting it at 270
    input = input * alpha_e/ex_tau * ex_tau^2; % transformed input
    %       ~~~~~   ~~~~~~~~~~~~~~   ~~~~~~~~~
    %       input   synaptic gain    integral of kernel

    if mV == 1   % adjust units
        Q = V2mVfactor^2 * Q;
        R = V2mVfactor^2 * R;

        r = r/V2mVfactor;
        varsigma = 1.699/r;                                                     % (spikes/(Vs))
        varsigma_sq = varsigma^2;
        v0 = v0*V2mVfactor;
        alpha_e = alpha_e*V2mVfactor;                                           % gain of excitatory synapses (V)
        alpha_i = alpha_i*V2mVfactor;                                           % gain of inhibitory synapses (V)

        input= input*V2mVfactor;
    end

    % conectivity constants to relate to Jansen and Rit 1995 model
    %
    ConnectivityConst = 270;                            % Jansen and Rit connectivity parameters. Either 135, 270 or 675
    C1 = ConnectivityConst;
    C2 = 0.8*ConnectivityConst;
    C3 = 0.25*ConnectivityConst;
    C4 = 0.25*ConnectivityConst;

    %% model structure
    % ~~~~~~~~~~~~~~~~
    %           X
    %       __  |  __
    %      /  \ | /  \
    %     /  04 | 01  \
    %     |     P     |
    %  ^  |     | |   |  ^
    %  |  E     | v   I  |  direction of info
    %     03   /|\   02
    %     |   / | \   |
    %      \_/  |  \_/
    %           v
    % population types: E, P, I, X
    % synapses: 01 (IP), 02 (PI), 03 (PE), 04 (EP)

    %% this is the observation function.
    %
    H = zeros(1,N_states);        %Initialize to zeros and later add 1s to states that contribute to EEG

    % initialize adjancy matrix
    %
    Gamma = zeros(2*N_syn+N_inputs,2*N_syn+N_inputs);   %  - plus 1 for input

    % specify synapses
    %
    syn_index = 0;

    % syn1, connection from I to P  % x = [ve ze vp1 zp1 vp2 zp2 vi zi];
    %
    syn_index = syn_index + 1;
    tau(syn_index) = in_tau;
    alpha(syn_index) = alpha_i*2*f_max*C4*dt / tau(syn_index);          % note the time constant and time step are in the gains
    presyn_inputs = 2;                                                  % the presynaptic population is getting inputs from synapses 2
    if ~isempty(presyn_inputs)
        Gamma(2*syn_index,2*presyn_inputs-1) = 1;                       % set the entries of Gamma corresponding to indices of presynaptic inputs to 1
    end
    H(2*syn_index-1) = 1;

    % syn2, connection from P to I
    %
    syn_index = syn_index + 1;
    tau(syn_index) = ex_tau;
    alpha(syn_index) = alpha_e*2*f_max*C3*dt / tau(syn_index);          % note the time constsnt and time step are in the gains
    presyn_inputs = [1 4 5];                                            % the presynaptic population is getting inputs from synapses 1, 4, 5
    if ~isempty(presyn_inputs)
        Gamma(2*syn_index,2*presyn_inputs-1) = 1;
    end
    H(2*syn_index-1) = 0;                                               % set to one if it contributes to the EEG (i.e. if the synapse is to Py cells)

    % syn3, connection from P to E
    %
    syn_index = syn_index + 1;
    tau(syn_index) = ex_tau;
    alpha(syn_index) = alpha_e*2*f_max*C1*dt / tau(syn_index);
    presyn_inputs = [1 4 5];                                         	% the presynaptic population is getting inputs from no other synapses (in the model)
    if ~isempty(presyn_inputs)
        Gamma(2*syn_index,2*presyn_inputs-1) = 1;
    end
    H(2*syn_index-1) = 0;

    % syn4, connection from E to P
    %
    syn_index = syn_index + 1;
    tau(syn_index) = ex_tau;
    alpha(syn_index) = alpha_e*2*f_max*C2*dt / tau(syn_index);          % note the time constsnt and time step are in the gains
    presyn_inputs = 3;                                                  % the presynaptic population is getting inputs from synapse 3
    if ~isempty(presyn_inputs)
        Gamma(2*syn_index,2*presyn_inputs-1) = 1;
    end
    H(2*syn_index-1) = 1;

    % for input
    %
    syn_index = syn_index + 1;
    H(2*syn_index-1) = 1;           % the input contributes to the observation function

    if N_inputs > 1
        % offset term
        H(2*syn_index) = 1;            % offset contributes to the observation function
    end

    % rescale H
    %
    H = H/scale;                                                        % !scale! this is help deal with our numerical issues.

    % set_ABC.m

    % dean freestone

    % this script take the parameters from set_params and creates the system
    % matrices for the neural mass model (recurrent NN)

    %% define A
    % A encodes the dynamics induced by the membrane time constants.
    % A is made up of the submatrices Psi in a block diagonal structure.
    % There is a Psi for each connection in the model. This is where all the
    % synaptic time constants enter the system. Further, the scale paramter
    % enters here (and with C (multiplicative factor) and with the H (divisor).

    Psi = zeros(2*N_syn,2*N_syn);               % initialise Psi, the component of A for fast states
    for n=1:N_syn                               % build block diagonal structure
        index = 2*(n-1)+1;
        Psi(index:index+1,index:index+1) = [0 scale ; -1/(scale*tau(n)^2) -2/(tau(n))];
    end

    % A = [1 + dt*Psi, 0
    %      0           1]
    % where 1 is the identity matrix of the appropriate size
    %
    A = [(eye(2*N_syn) + dt*Psi) zeros(2*N_syn,N_syn+N_inputs) ; zeros(N_syn+N_inputs,2*N_syn) eye(N_syn+N_inputs)]; % top left part is discrete time versions.

    %% define B (syanptic gain selection)
    %
    % Theta = [0 0 ... 0
    %          1 0 ... 0
    %          0 0 ... 0
    %          0 1 ... 0
    %          ...
    %          0 0 ... 0
    %          0 0 ... 1]
    % B = [0 Theta ; 0 0]

    Theta = zeros(2*N_syn,N_syn);                   % THETA IS USED TWICE
    for n=1:N_syn
        index = 2*(n-1)+1;
        Theta(index:index+1,n) = [0 ; 1];
    end
    B = [zeros(2*N_syn,2*N_syn+N_inputs) Theta ; zeros(N_syn+N_inputs,3*N_syn+N_inputs)];

    %% define C (adjacency matrix)
    %
    C = [Gamma/scale zeros(2*N_syn+N_inputs,N_syn) ; zeros(N_syn,3*N_syn+N_inputs)];

    %%
    %%
    % WE ARE NOW BORN TO RUN
    % ~~~~~~~~~~~~~~~~~~~~~~~
    % set up the forward simulation
    %

    % define initial conditions
    %
    ve = 0; ze = 0; vp1 = 0; zp1 = 0; vp2 = 0; zp2 = 0; vp3 = 0; zp3 = 0; vi = 0; zi = 0;
    x = [ve ze vp1 zp1 vp2 zp2 vi zi];

    xi = zeros(N_states,N_samples);
    xi(:,1) = [x input alpha]';                 % note: input and alpha are set in set_params

    % rng(1)
    w = mvnrnd(zeros(1,N_states),Q,N_samples)';

    % phi_p = zeros(1,N_samples);

    % run the model forward !!!This is the model
    %
    for n=1:N_samples-1

        phi = g(C*xi(:,n),v0,varsigma); % g is error function sigmoid
    %     phi_p(n) = phi(4);
        xi(:,n+1) = A*xi(:,n) + B*xi(:,n).*phi;% + w(:,n);

    end

    % rng(1)
    v = sqrt(R)*randn(1,N_samples); % variance
    s_y = H*xi; %+ v;                           % simulated EEG signal
end

% error function sigmoid    
function out = g(v,v0,varsigma)
    out = 0.5*erf((v - v0) / (sqrt(2)*varsigma)) + 0.5;
end