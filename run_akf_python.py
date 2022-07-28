import numpy as np
import numpy.matlib
import random
import math
import mat73
from tqdm import tqdm
import os
from scipy import signal
from scipy.io import savemat
import pyarrow as pa
import pyarrow.parquet as pq
import warnings
warnings.filterwarnings('ignore')

### Define some functions
def isPD(B):
    """Returns true when input is positive-definite, via Cholesky"""
    try:
        _ = np.linalg.cholesky(B)
        return True
    except np.linalg.LinAlgError:
        return False

def nearest_spd(A):
    k = 1
    B = (A + A.T)/2
    _, sigma, V = np.linalg.svd(B)
    H = np.dot(V.T, np.dot(np.diag(sigma),V))
    Ahat = (B+H)/2
    Ahat = (Ahat + Ahat.T)/2
    if isPD(Ahat):
        return Ahat, k
    spacing = np.spacing(np.linalg.norm(A))
    I = np.eye(A.shape[0])
    k = 1
    while not isPD(Ahat):
        mineig = np.min(np.real(np.linalg.eigvals(Ahat)))
        Ahat += I * (-mineig * k ** 2 + spacing)
        k += 1
        if k >= 1e5:
            k = -1
            return [Ahat, k]
    return [Ahat, k]

# NMM parameter initialization
def set_params(ext_input, input_offset, TimeOfSim, Fs, sigma_R = 1e-3):
    """
    set_params, migrated from the MATLAB version by Pip Karoly
    Set parameters for the neural mass model

    Inputs:
        ext_input - the input to the model
        input_offset - value of the offset (to compensate for a DC offset if required - i.e. there is a DC shift in the model but data may not be recorded with DC component)
        TimeOfSim - length of time to simulate data for
        Fs - sampling frequency (Hz)
        sigma_R (optional) - Default value: 1e-3

    Outputs:
        A,B,C,H: model and observation matrices (defined in Karoly et al 2018)
        N_states,N_syn,N_inputs,N_samples: model dimensions
        xi, y: simulated data (xi = state vector, y = measurement)
        v0,varsigma: model constants
        Q,R: model and measurement noise

    Neural mass model parameters have been modified from Jansen & Rit (1995)

    For further references see:

        [1] Freestone, D. R., Karoly, P. J., Neši?, D., Aram, P., Cook, M. J., & Grayden, D. B. (2014).
        Estimation of effective connectivity via data-driven neural modeling. Frontiers in neuroscience, 8, 383

        [2] Ahmadizadeh, S., Karoly, P. J., Neši?, D., Grayden, D. B., Cook, M. J., Soudry, D., & Freestone, D. R. (2018).
        Bifurcation analysis of two coupled Jansen-Rit neural mass models. PloS one, 13(3), e0192842.

        [3] Kuhlmann, L., Freestone, D. R., Manton, J. H., Heyse, B., Vereecke, H. E., Lipping, T., ... & Liley, D. T. (2016).
        Neural mass model-based tracking of anesthetic brain states. NeuroImage, 133, 438-456.
    """

    scale = 50 # This is to get states and derivatives on the same order of magnitude

    mV = True # Set units of membrane potentials (True for mV, False for V)
    V2mVfactor = 1e3

    # Time parameters
    dt = 1/Fs #Time step (s)
    N_samples = int(np.round(TimeOfSim/dt)) # No. of time poitns to simulate

    input_offset = np.array(input_offset) # Making sure it is a NumPy array
    if input_offset.size > 0:
        N_inputs = 2
    else:
        N_inputs = 1

    N_syn = 4 # No. of synapses
    N_states = 3 * N_syn + N_inputs # No. of states / dimension of the state-space. Plus one for input

    # Define the distrubance covariance matrix
    sigma_all = 5e-8
    sigma_input = 5e-4
    sigma_params = 5e-5
    sigma_offset = 5e-6


    Q = np.identity(N_states) * (scale * np.sqrt(dt) * sigma_all)**2 # add a tiny bit of noise to all states (added for numerical stability)
    Q[2*N_syn:, 2*N_syn:] = np.identity(N_syn + N_inputs) * (scale * np.sqrt(dt) * sigma_params)**2
    Q[2*N_syn, 2*N_syn] = (scale * np.sqrt(dt) * sigma_input)**2
    if N_inputs > 1:
        Q[2*N_syn+1, 2*N_syn+1] = (scale * np.sqrt(dt) * sigma_offset)**2

    # Measurement disturbance covariance
    R = sigma_R**2

    # General parameters from Jansen and Rit
    # Sigmoid bits
    f_max = 2.5 # maximum firing rate (spikes/s)
    r = 560
    varsigma = 1.699/r # spikes/(Vs)
    varsigma_sq = varsigma**2 # V

    # Synaptic gains
    alpha_e = 3.25e-3   # gain of excitatory synapses (V)
    alpha_i = -22e-3    # gain of inhibitory synapses (V)
    v0 = 0.006

    # Synaptic kernel time constants
    ex_tau = 0.010  # excitatory synaptic time constant (s)
    in_tau = 0.020  # inhibitory synaptic time constant (s)

    # input to py population
    #
    # SCALE 1 - this is to avoid large differences between states upsetting the filter
    # (magnitude of the membrane potentials and their derivatives)
    ext_input = ext_input*scale
    # SCALE 2 - this converts a constant input to its effect on the pyramidal
    # membrane potential by taking the steady state limit of the synaptic kernel
    # (assumption that the input varies much slower than the state variables).
    ext_input = ext_input * alpha_e/ex_tau * ex_tau**2
    #       ~~~~~   ~~~~~~~~~~~~~~   ~~~~~~~~~
    #       input   synaptic gain    integral of kernel

    # measurement DC offset
    input_offset = input_offset * scale
    input_offset = input_offset * alpha_e/ex_tau * ex_tau**2

    if mV:
        Q = (V2mVfactor**2) * Q
        R = (V2mVfactor**2) * R

        r = r/V2mVfactor
        varsigma = 1.699/r  # (spikes/(Vs))
        varsigma_sq = varsigma**2
        v0 = v0*V2mVfactor
        alpha_e = alpha_e*V2mVfactor    # gain of excitatory synapses (V)
        alpha_i = alpha_i*V2mVfactor    # gain of inhibitory synapses (V)

        ext_input= ext_input*V2mVfactor
        input_offset = input_offset*V2mVfactor

    # Connectivity constants to relate to Jansen and Rit 1995 model
    ConnectivityConst = 270     # Jansen and Rit connectivity parameters. Either 135, 270 or 675
    C1 = ConnectivityConst
    C2 = 0.8*ConnectivityConst
    C3 = 0.25*ConnectivityConst
    C4 = 0.25*ConnectivityConst

    #   Model structure
    #  ~~~~~~~~~~~~~~~~~
    #
    #           X
    #       __  |  __
    #      /  \ | /  \
    #     /  04 | 01  \
    #     |     P     |
    #  ^  |     | |   |  ^
    #  |  E     | v   I  |  direction of info
    #     03   /|\   02
    #     |   / | \   |
    #      \_/  |  \_/
    #           v
    # population types: E, P, I, X
    # synapses: 01 (IP), 02 (PI), 03 (PE), 04 (EP)

    # Initialize some variables
    tau = np.zeros([4,])
    alpha = np.zeros([4,])

    # This is the Observation function
    H = np.zeros([1,N_states]) # Initialize to zeros and later add 1s to states that contribute to EEG

    # Initialize adjancy matrix
    Gamma = np.zeros([2 * N_syn + N_inputs, 2 * N_syn + N_inputs]) # - plus 1 for input

    # Specify synapses
    syn_index = 0 # Python indexing starts at 0, hence syn_index 0 is equivalent to the syn_index 1 in Matlab

    # Syn1, connection from I to P
    tau[syn_index,] = in_tau
    alpha[syn_index,] = alpha_i * 2 * f_max * C4 * dt / tau[syn_index,] # note the time constant and time step are in the gains
    presyn_inputs = np.array([2]) # the presynaptic population is getting inputs from synapses 2
    # Check if presyn_inputs not empty
    if presyn_inputs.size > 0:
        """ Note! the folllowing line was brought from Matlab code. -1 was
        turned into -2 to comply with Python indexing"""
        Gamma[2 * (syn_index + 1) -1, 2 * presyn_inputs -2] = 1 # set the entries of Gamma corresponding to indices of presynaptic inputs to 1
    H[0, 2 * (syn_index + 1) -2] = 1

    # Syn2, connection from P to I
    syn_index = syn_index + 1
    tau[syn_index,] = ex_tau
    alpha[syn_index,] = alpha_e * 2 * f_max * C3 * dt / tau[syn_index,] # note the time constsnt and time step are in the gains
    presyn_inputs = np.array([1, 4, 5]) # the presynaptic population is getting inputs from synapses 1, 4, 5
    if presyn_inputs.size > 0:
        Gamma[2 * (syn_index + 1) -1, 2 * presyn_inputs -2] = 1
    H[0, 2 * (syn_index + 1) -2] =  0 # set to one if it contributes to the EEG (i.e. if the synapse is to Py cells)

    # Syn3, connection from P to E
    syn_index = syn_index + 1
    tau[syn_index,] = ex_tau
    alpha[syn_index,] = alpha_e * 2 * f_max * C1 * dt / tau[syn_index,] # note the time constsnt and time step are in the gains
    presyn_inputs = np.array([1, 4, 5]) # the presynaptic population is getting inputs from synapses 1, 4, 5
    if presyn_inputs.size > 0:
        Gamma[2 * (syn_index + 1) -1, 2 * presyn_inputs -2] = 1
    H[0, 2 * (syn_index + 1) -2] =  0

    # Syn4, connection from E to P
    syn_index = syn_index + 1
    tau[syn_index,] = ex_tau
    alpha[syn_index,] = alpha_e * 2 * f_max * C2 * dt / tau[syn_index,] # note the time constsnt and time step are in the gains
    presyn_inputs = np.array([3]) # the presynaptic population is getting inputs from synapses 3
    if presyn_inputs.size > 0:
        Gamma[2 * (syn_index + 1) -1, 2 * presyn_inputs -2] = 1
    H[0, 2 * (syn_index + 1) -2] =  1

    # For input
    syn_index = syn_index + 1
    H[0, 2 * (syn_index + 1) -2] =  1 # The input contributes to the observation function

    if N_inputs > 1:
        # Offset term
        H[0, 2 * (syn_index + 1) -1] =  1 # Offset contributes to the observation function. Notice the -1 instead of -2 in the indexing.

    # Rescale
    H = H/scale # Scale! This helps deal with our numerical issues.

    # Define A
    #
    # A is made up of the submatrices Psi in a block diagonal structure.
    # There is a Psi for each connection in the model. This is where all the
    # synaptic time constants enter the system. Further, the scale paramter
    # enters here (and with C (multiplicative factor) and with the H (divisor).
    Psi = np.zeros([2*N_syn, 2*N_syn]) # initialise Psi, the component of A for fast states
    for n in range(0,N_syn):
        index = 2*n
        Psi[index : index + 2, index : index + 2] = np.array([[0, scale],
                                                              [-1/(scale * tau[n]**2), -2/tau[n]]])

    # A = [1+dt*Psi, 0
    #          0   , 1]
    # Where 1 is the identity matrix of the appropriate size.
    a11 = np.identity(2*N_syn) + dt*Psi # [1]+dt*Psi
    a12 = np.zeros([2*N_syn, N_syn+N_inputs]) # [0]
    a21 = np.zeros([N_syn+N_inputs, 2*N_syn]) # [0]
    a22 = np.identity(N_syn+N_inputs) # [1]
    # Concatenate horizontally
    a1 = np.concatenate((a11,a12), axis=1)
    a2 = np.concatenate((a21,a22), axis=1)
    # Concantenate vertically
    A = np.concatenate((a1,a2))

    # Define B
    #
    # Theta = [0 0 ... 0
    #          1 0 ... 0
    #          0 0 ... 0
    #          0 1 ... 0
    #              ...
    #          0 0 ... 0
    #          0 0 ... 1]
    # B = [0 Theta ; 0 0]
    Theta = np.zeros([2*N_syn, N_syn]) # Theta is used twice!
    for n in range(0,N_syn):
        index = 2*n
        Theta[index : index + 2, n : n+1] = np.array([[0],[1]])
    b1 = np.concatenate((np.zeros([2*N_syn, 2*N_syn+N_inputs]), Theta), axis=1)
    b2 = np.zeros([N_syn+N_inputs, 3*N_syn+N_inputs])
    B = np.concatenate((b1,b2))

    # Define C (adjacency matrix)
    c1 = np.concatenate((Gamma/scale, np.zeros([2*N_syn + N_inputs, N_syn])), axis=1)
    c2 = np.zeros([N_syn, 3 * N_syn + N_inputs])
    C = np.concatenate((c1,c2))

    # Model structure has been defined. WE ARE NOW BORN TO RUN!
    # Set up forward simulation

    # Define initial conditions
    ve = 0; ze = 0; vp1 = 0; zp1 = 0; vp2 = 0; zp2 = 0; vi = 0; zi = 0; # vp3 = 0; zp3 = 0; # <- unused variables. Legacy?
    x = np.array([ve, ze, vp1, zp1, vp2, zp2, vi, zi])
    xi = np.zeros([N_states, N_samples])
    i_off = np.array([input_offset])
    xi[:,0] = np.concatenate((x, np.array([ext_input]), np.reshape(i_off, (i_off.size,)), alpha)) # ext_input and alpha are set in set_params

    np.random.seed(1)
    w = np.random.multivariate_normal(mean = np.zeros(N_states),
                                       cov = np.array(Q),
                                       size = N_samples)
    w = np.array(w)
    w = np.transpose(w)

    phi_p = np.zeros([1, N_samples])

    # Run the model forward
    for n in range(0, N_samples-1):
        v = np.matmul(C,xi[:,n])
        phi = g(v,v0,varsigma)
        phi_p[0,n:n+1] = phi[3:4,0]
        xi[:,n+1:n+2] = np.matmul(A, xi[:,n:n+1]) + np.multiply(np.matmul(B, xi[:,n:n+1]), phi) + w[:,n:n+1]

    np.random.seed(1)
    v = np.sqrt(R) * np.random.randn(1, N_samples)
    y = np.matmul(H,xi) + v

    # Define output arguments
    output = {
        'A': A,
        'B': B,
        'C': C,
        'H': H,
        'N_states': N_states,
        'N_syn': N_syn,
        'N_inputs': N_inputs,
        'N_samples': N_samples,
        'xi': xi,
        'y': y,
        'v0': v0,
        'varsigma': varsigma,
        'Q': Q,
        'R': R,
    }
    return output

# Error function sigmoid
def g(v, v0, varsigma):
    """ g(.) is the sigmoidal activation function """

    v_shape = np.array(v.shape)
    if v_shape.size > 0:
        g_ = np.zeros(v_shape)
        for i in range(0, v_shape[0]):
            g_[i] = 0.5 * math.erf((v[i] - v0) / (np.sqrt(2)*varsigma)) + 0.5
    else:
        raise ValueError("The size of input parameter 'v' must be greater than zero.")

    g_ = np.reshape(g_, (g_.size,1))

    return g_

# Propagate mean and covariance
def propagate_metrics(N_syn, N_states, N_inputs, A, B, C, P_0p, xi_0p, varsigma, v0, Q):
    """
    PROPAGATE MEAN AND COVARIANCE TRHOUGH NMM EQUATIONS

    Results for the posterior mean and covariance of the multivariate
    distribution over the membrane potentials (v), derivatives (z) and alpha
    gain terms (alpha) in a NMM.

    Authors:
        Dean Freestone, Philippa Karoly 2016

    This code is licensed under the MIT License 2018


    Parameters
    ----------
    N_syn :
        N synapses in the model.
    N_states :
        N estimation states.
    N_inputs :
        N external inputs.
    A :
        A - Matrix describing the neural mass model equations.
    B :
        B - Matrix describing the neural mass model equations.
    C :
        C - Matrix describing the neural mass model equations.
    P_0p :
        Initial covariance matrix.
    xi_0p :
        Initial mean matrix.
    varsigma :
        Sigmoid function parameter - variance.
    v0 :
        Sigmoid function parameter - initial value.
    Q :
        Model unvertainty matrix (same dimensions as P).

    Returns
    -------
    xi_1m :
        Posterior mean matrix.
    P_1m :
        Posterior covariance matrix.

    References
    ----------
    Further details on the estimation method can be found in
    the following references:

    [1] Freestone, D. R., Karoly, P. J., Neši?, D., Aram, P., Cook, M. J., & Grayden, D. B. (2014).
    Estimation of effective connectivity via data-driven neural modeling. Frontiers in neuroscience, 8, 383

    [2] Ahmadizadeh, S., Karoly, P. J., Neši?, D., Grayden, D. B., Cook, M. J., Soudry, D., & Freestone, D. R. (2018).
    Bifurcation analysis of two coupled Jansen-Rit neural mass models. PloS one, 13(3), e0192842.

    [3] Kuhlmann, L., Freestone, D. R., Manton, J. H., Heyse, B., Vereecke, H. E., Lipping, T., ... & Liley, D. T. (2016).
    Neural mass model-based tracking of anesthetic brain states. NeuroImage, 133, 438-456.

    """
    # Give

    v_indexes = np.array(range(0, 2*N_syn+1, 2))
    z_indexes = np.array(range(1, 2*N_syn+1, 2))
    alpha_indexes = np.array(range(2*N_syn+N_inputs, N_states))

    w_fast = np.reshape(np.array([0.1713244923791705, 0.3607615730481384,  0.4679139345726904, 0.1713244923791705, 0.3607615730481384, 0.4679139345726904]), [1,6])
    y_fast = np.array([0.033765242898424, 0.169395306766868, 0.380690406958402, 0.966234757101576, 0.830604693233132, 0.619309593041599])

    CPC = np.matmul(np.matmul(C[z_indexes, :], P_0p) , np.transpose(C[z_indexes, :])) # CPC' # Consider adding np.newaxis as in Bxi line
    dCPC = np.diag(CPC)
    dCPB = np.diag(np.matmul(np.matmul(C[z_indexes, :], P_0p) , np.transpose(B[z_indexes, :]))) #diagonal of CPB'
    Bxi = np.matmul(B[z_indexes[:, np.newaxis], alpha_indexes[np.newaxis,:]], xi_0p[alpha_indexes])
    AP = np.matmul(A, P_0p)

    # Analytic mean
    gamma = np.divide(1,np.sqrt(abs(2*(dCPC + varsigma**2))))# gamma = np.divide(1,np.sqrt(2*(dCPC + varsigma**2))) # gamma = np.divide(1,np.sqrt(abs(2*(dCPC + varsigma**2)))) #
    beta = (np.matmul(C[z_indexes[:, np.newaxis], v_indexes[np.newaxis, :]], xi_0p[v_indexes]) - v0) * gamma
    Xi = np.zeros(beta.shape)
    for i in range(0, Xi.size):
        Xi[i] = (math.erf(beta[i]) + 1)/2

    Upsilon = np.divide((np.exp( -(pow(np.array(beta), 2))) * gamma), np.sqrt(np.pi))
    psi = (Bxi * Xi) + (dCPB * Upsilon)

    xi_1m = np.matmul(A, xi_0p)
    xi_1m[np.array(range(1, 2*N_syn,2))] = xi_1m[np.array(range(1, 2*N_syn,2))] + psi

    # Exact covariance
    # cov part 1 (Phi)
    q2 = Upsilon * (Bxi - dCPB * beta * (gamma ** 2) * 2)
    Phi_ = np.matmul(np.ones([N_states, 1]), np.reshape(q2, [1, q2.size])) * np.matmul(AP, np.transpose(C[z_indexes,:]))
    Phi = Phi_ + ( np.matmul(np.ones([N_states, 1]), np.reshape(Xi, [1, Xi.size])) * np.matmul(AP, np.transpose(B[z_indexes,:])) )

    # cov part 2 (Omega)
    # NOTE: we can reduce the dimensionality (only compute upper triangle
    # part) of the matrices we turn to vectors to increase speed in larger
    # networks.
    with np.errstate(invalid='raise'):
        try:
            CPCgammagamma = np.arcsin(CPC * np.matmul(np.reshape(gamma, [gamma.size, 1]), np.reshape(gamma, [1, gamma.size]))*2)
        except:
            CPC_lower_1 = np.array(CPC)
            CPC_lower_1[abs(CPC) > 1] = CPC[abs(CPC) > 1]/abs(CPC[abs(CPC) > 1])
            CPCgammagamma = np.arcsin(CPC_lower_1 * np.matmul(np.reshape(gamma, [gamma.size, 1]), np.reshape(gamma, [1, gamma.size]))*2)

    CPCgammagamma = np.reshape(CPCgammagamma, [CPCgammagamma.size,1]) # Change matrix into a vector
    CPCgammagamma_y = np.matmul(CPCgammagamma, np.reshape(y_fast, [1, y_fast.size])) # a Matrix of row vectors

    # Warning! Check if x2 should go or not. It goes in MATLAB version
    betabeta = np.matmul(np.reshape(beta, [beta.size, 1]), np.reshape(beta, [1, beta.size]))*2
    betabeta_mat = np.matmul(np.reshape(betabeta, [betabeta.size, 1]), np.ones([1, w_fast.size]))

    beta2mat = np.matlib.repmat(np.reshape(beta**2, [beta.size, 1]), 1, beta.size) # square and repmat
    beta2matT = np.transpose(beta2mat) # Transpose allows for permutation when we sum below
    b2b2T =  np.reshape(beta2matT, [beta2mat.size, 1]) + np.reshape(beta2mat, [beta2mat.size, 1]) # Change to vector, add, repmat
    b2b2T = np.matlib.repmat(b2b2T, 1, w_fast.size)

    # Put it together
    # Psi = reshape(sum(CPCgammagamma*w_fast.*exp(-(beta2_plus_beta2T - betabeta_mat.*sin(CPCgammagamma_y))./cos(CPCgammagamma_y).^2),2),N_syn,N_syn)/(4*pi);
    Psi = np.matmul(CPCgammagamma, w_fast) * np.exp(-(b2b2T - betabeta_mat * np.sin(CPCgammagamma_y))/(np.cos(CPCgammagamma_y ** 2)))
    Psi = np.reshape(np.sum(Psi, axis=1), [N_syn, N_syn]) / (4*np.pi)

    Xi = np.reshape(Xi, [Xi.size, 1]) # Fix Xi dimmensions
    Bxi = np.reshape(Bxi, [Bxi.size, 1]) # Fix Bxi dimmensions
    Omega = (np.matmul(Xi,np.transpose(Xi)) + Psi) * (np.matmul(Bxi,np.transpose(Bxi)) + P_0p[alpha_indexes[np.newaxis, :], alpha_indexes[:, np.newaxis]])
    #                  ~~~~~~~~~~~~~~~~~~~~~~~~~~                ~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                       E[g^2(Cxi)]           (alpha^2 + sigma_alpha^2)

    # Here we construct the cov matrix:
    psi = np.reshape(psi, [psi.size, 1]) # Fix psi dimmesnions
    P_1m = np.matmul(AP, np.transpose(A)) + Q
    P_1m[z_indexes[np.newaxis,:], z_indexes[:, np.newaxis]] = P_1m[z_indexes[np.newaxis,:], z_indexes[:, np.newaxis]] + Omega - np.matmul(psi, np.transpose(psi))
    P_1m[:, z_indexes] = P_1m[:, z_indexes] + Phi
    P_1m[z_indexes, :] = P_1m[z_indexes, :] + np.transpose(Phi)

    # Define output arguments
    output = {
        'xi_1m': xi_1m,
        'P_1m': P_1m,
    }
    return output

# save estimates to parquet files
def save_estimate(estimate, save_file_name):
    arrays = [
        pa.array(col)  # Create one arrow array per column
        for col in estimate
    ]

    table = pa.Table.from_arrays(
        arrays,
        names=[str(i) for i in range(len(arrays))] # give names to each columns
    )
    pq.write_table(table,save_file_name)

# load estimate from parquet files
def load_estimate(save_file_name):
    table_from_parquet = pq.read_table(save_file_name)
    matrix_from_parquet = table_from_parquet.to_pandas().T.to_numpy()
    return matrix_from_parquet

if __name__ == "__main__":
    # This part is hardcoded for a demo case
    # Please modify it according to your needs
    # You can use your own data by changing the following lines [503-511]
    print('##################################################')
    print('This demonstration requires a lot of storage space \n due to the parameter estimates file is very big.')
    print('##################################################')

    if not os.path.exists('./output/21'):
        os.makedirs('./output/21')

    data_dict = mat73.loadmat('./data/data_21.mat')
    data = data_dict['virtualdata_timeseries']
    time = 50
    Fs = 150
    channels = 20 # only take the first 5 channels as a demonstration
    downsampled_size = 33375 # downsample the meg data

    aEP_collection = np.zeros([channels, downsampled_size])
    aIP_collection = np.zeros([channels, downsampled_size])
    aPE_collection = np.zeros([channels, downsampled_size])
    aPI_collection = np.zeros([channels, downsampled_size])
    mu_collection = np.zeros([channels, downsampled_size])
    v_e_collection = np.zeros([channels, downsampled_size])
    v_i_collection = np.zeros([channels, downsampled_size])
    v_p_collection = np.zeros([channels, downsampled_size])

    # Run the filter for demonstration
    for iCh in range(0,channels):
        print('Channel %02d ...' % (iCh)) # Print current channel
        y = data[iCh,:]
        y = 5.4724e+12 * y # scale the data to the same magnitude of the model output
        y = signal.resample(y, downsampled_size) # downsample the data

        # Initialize input
        ext_input = 300#300 # External input
        input_offset = np.empty(0)

        # Parameter initialization
        params = set_params(ext_input, input_offset, time,Fs)

        # Retrive parameters into single variables
        A = params['A']
        B = params['B']
        C = params['C']
        H = params['H']
        N_inputs = params['N_inputs']
        N_samples = params['N_samples']
        N_states = params['N_states']
        N_syn = params['N_syn']
        Q = params['Q']
        R = params['R']
        v0 = params['v0']
        varsigma = params['varsigma']
        xi = params['xi']

        xi_hat_init = np.mean( params['xi'][:, int(np.round(N_samples/2))-1:] , axis = 1)
        P_hat_init = 10 * np.cov(params['xi'][:, int(np.round(N_samples/2))-1:])
        P_hat_init[2*N_syn:, 2*N_syn:] = np.eye(N_syn + N_inputs) * 10e-2

        # Set initial conditions for the Kalman Filter
        xi_hat = np.zeros([N_states, N_samples])
        P_hat = np.zeros([N_states, N_states, N_samples])
        P_diag = np.zeros([N_states, N_samples])

        xi_hat[:,0] = xi_hat_init
        P_hat[:,:,0] = P_hat_init

        anneal_on = 1 # Nice!
        kappa_0 = 10000
        T_end_anneal = N_samples/20

        N_samples = y.size

        # Set initial conditions for the Kalman Filter
        xi_hat = np.zeros([N_states, N_samples])
        P_hat = np.zeros([N_states, N_states, N_samples])
        P_diag = np.zeros([N_states, N_samples])
        xi_hat[:,0] = xi_hat_init
        P_hat[:,:,0] = P_hat_init

        try:
            for t in range(1,N_samples):

                xi_0p = xi_hat[:, t-1].squeeze()
                P_0p = P_hat[:, :, t-1].squeeze()

                # Predict
                metrics = propagate_metrics(N_syn, N_states, N_inputs, A, B, C, P_0p, xi_0p, varsigma, v0, Q)
                xi_1m = metrics['xi_1m']
                P_1m = metrics['P_1m']

                if (t <= T_end_anneal) & (anneal_on):
                    kappa = pow(kappa_0, (T_end_anneal-t-1)/(T_end_anneal-1))
                else:
                    kappa = 1

                # K = P_1m*H'/(H*P_1m*H' + kappa*R);
                K = np.divide(np.matmul(P_1m, np.transpose(H)), np.matmul(H, np.matmul(P_1m, np.transpose(H))) + kappa*R)

                # Correct
                xi_1m = np.reshape(xi_1m, [xi_1m.size, 1])
                xi_hat[:, t:t+1] = xi_1m + K*(y[t] - np.matmul(H, xi_1m))

                P_hat[:,:,t] = np.matmul((np.identity(N_states) - np.matmul(K,H)), P_1m)
                try:
                    P_hat[:,:,t] = (P_hat[:,:,t] + np.transpose(P_hat[:,:,t]))/2
                    chol_matrix = np.linalg.cholesky(P_hat[:,:,t])
                except(np.linalg.LinAlgError):
                    P_hat[:,:,t] , k = nearest_spd(P_hat[:,:,t])
                    if k == -1:
                        print(file, 'cannot find PSD')
                P_diag[:,t] = np.diag(P_hat[:,:,t])

            aEP_collection[iCh,:] = xi_hat[12,:]
            aIP_collection[iCh,:] = xi_hat[9,:]
            aPE_collection[iCh,:] = xi_hat[11,:]
            aPI_collection[iCh,:] = xi_hat[10,:]
            mu_collection[iCh,:] = xi_hat[8,:]
            v_e_collection[iCh,:] = xi_hat[4,:]
            v_i_collection[iCh,:] = xi_hat[2,:]
            v_p_collection[iCh,:] = np.matmul(H,xi_hat)
        except:
            continue

# save estimates to parquet files to minimise the storage required
    save_estimate(aEP_collection, './output/21/aEP_estimate.pq')
    save_estimate(aIP_collection, './output/21/aIP_estimate.pq')
    save_estimate(aPE_collection, './output/21/aPE_estimate.pq')
    save_estimate(aPI_collection, './output/21/aPI_estimate.pq')
    save_estimate(mu_collection, './output/21/mu_estimate.pq')
    save_estimate(v_e_collection, './output/21/v_e_estimate.pq')
    save_estimate(v_i_collection, './output/21/v_i_estimate.pq')
    save_estimate(v_p_collection, './output/21/v_p_estimate.pq')
