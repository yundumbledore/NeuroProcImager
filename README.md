# NeuroProcImager
#### Yun Zhao, Levin Kuhlmann (Monash University, Australia), Email: yun.zhao1@monash.edu, levin.kuhlmann@monash.edu

**NeuroProcImager** is associated with the Neuroimage paper *Space-time resolved inference-based neurophysiological process imaging: Application to resting-state alpha rhythm*. https://www.sciencedirect.com/science/article/pii/S1053811922007078

In general, **NeuroProcImager** is an inferenced-based whole-brain imaging framework, aiming at visualising nuerophysiological variables and thereby providing a novel way to understand neural mechanisms underlying given brain state. In this repository, there are two critical components, 
1. Analytical Kalman Filter: estimates time-resolved neurophysiological variable in the neural mass model from the actual measurement
2. Whole-brain imaging module: visualises neurophysiological variable estimates or statistical observations derived from neurophysiological variable estimates in a whole brain manner.

In the below sections, you can find: 
- **Recommended system requirements**
- **Installation guide**
- **Demonstration with example data**
- **Common issues and solutions**
- **Work with your own data and visualisation packages**
- **Code reusability**

## Recommended system requirements
- Operating system: Windows 10
- MATLB version: R2020a or later
- Storage: 10 GB available space

**We recommend that you use a multi-core CPU, as the time efficiency depends on the number of CPU cores in your computer.**

## Installation guide
- To run this demo, you need to install "signal processing toolbox", "statistics and machine learning toolbox", "parallel computing toolbox" in the MATLAB.

- Fieldtrip version equals to or later than 20200828 is also required and available to download from https://www.fieldtriptoolbox.org. Another link is helpful for you to add fieldtrip path to your MATLAB https://www.fieldtriptoolbox.org/faq/should_i_add_fieldtrip_with_all_subdirectories_to_my_matlab_path/.

- After above steps, you download this package (NeuroProcImager), unzip and change the current working directory to this package in MATLAB Command Window and run 'defaults', which sets the defaults and configures up the minimal required path settings.

**Fieldtrip is optional to download which is exclusively for Case 2 visualisation purpose. You may use your own preferred visualisation packages, e.g., MNE-python (see 'Work with your own data and visualisation packages' below).**

## Demonstration with example data
The example data consists of 4714 MEG-derived source timeseries in the cerebral cortex from subject 21 and is available to download at https://drive.google.com/drive/folders/18EvFR4cr6YfhNUgZaijj9L3M1sG6cusL?usp=sharing. The data is to be put under the directory ‘/data’.

Before running showcases, you need to make sure that the current directory is ‘NeuroProcImager’. In this directory, you can find a file called ‘test.m’ and this is the script you want to use to run showcases. Note that there are three critical variables in 'test.m' to check before you run the script:
1.	‘pipeline’: put “AKF estimation” or “Contrast imaging” in the list to indicate which showcase to run.
2.	‘ncpu’:  number of CPU cores you want to use. By default, NeuroProcImager uses all cores in your CPU.
3.	‘multiple_comparison_correction’:  put 0 or 1 to indicate whether to perform corrections for multiple comparisons problem in Case 2.

### Case 1: Neurophysiological variables inference by AKF. 
This showcase demonstrates estimating neurophysiological variables in the neural mass model from 4714 MEG-derived source time-series by the AKF. You need to open the MATLAB script ‘test.m’ and put “AKF estimation” in the ‘pipeline’ list. You can run the script by entering the following command in the MATLAB Command Window: ‘test’.  The input to the framework is subject 21 MEG data containing 4714 virtual timeseries and the output is a (.mat) file containing neurophysiological variable estimates for every source point at ‘/output/variable_estimates_21.mat’. There are eight matrices in ‘/output/variable_estimates_21.mat’ and they are the neurophysiological variable estimates for α_ep, α_ip, α_pe, α_pi, μ, v_e, v_i, v_p. Each element contains 4714 rows representing 4714 MEG source points.

P.S. We have a python-version AKF called "run_akf_python.py". This is an alternative solution when Matlab is unavailable for some users. The program takes M/EEG as input and saves parameter estimates in the output folder.

### Case 2: Contrast imaging between strong and weak occipital alpha oscillation. 
This showcase is to image the contrast of neurophysiological variables between when the occipital alpha power is strong and when it is weak. One needs to open the MATLAB script ‘test.m’ and put “Contrast imaging” in the ‘pipeline’ list. The variable ‘multiple_comparison_correction’ toggles between performing corrections for multiple comparisons or not. Note that the correction implemented in this case study is exclusively for individual-level analysis and not the same as what we did in the group-level, but the principal mechanism is the same with group-level analysis. You can then run the script by entering the following command in the MATLAB Command Window: ‘test’. The input to the framework is the neurophysiological variable estimates derived in Case 1 and the output is imaging showing the contrast of neurophysiological variables during strong and weak occipital alpha power. The output images can be found under the directory ‘/demo_cases/contrast_imaging/visualisation_outputs’. Three dimensions to view the brain are provided: lateral, rear, and dorsal view.

**Note: Case 1 should be run prior to Case 2.**

#### Case 2 visulisation outputs (only shows lateral view with multiple comparison correction applied)
Figure legend    |         Figure snapshot                                                                                                                                                                                                                     
-------------------------------------- |--------------------------------------------------  
Individual-level (subject 21) contrast imaging of the alpha power during times of strong and weak occipital alpha rhythm. The colour bar refers to t-statistics indicating the mean difference between alpha power during strong and weak occipital alpha. Note this figure should not be interpreted in the context of hypothesis testing and only uses the t-statistic as a visualisation tool.           | ![](assets/sub_21_D1_v_pyr_corrected.jpg)
Individual-level (subject 21) contrast imaging of the excitatory to pyramidal connection strength, α_ep, during times of strong and weak occipital alpha rhythm, with corrections for multiple comparisons applied. The colour bar refers to t-statistics indicating the mean difference between connection strength estimates during strong and weak occipital alpha. | ![](assets/sub_21_D1_aEP_corrected.jpg)
Individual-level (subject 21) contrast imaging of the inhibitory to pyramidal connection strength, α_ip, during times of strong and weak occipital alpha rhythm, with corrections for multiple comparisons applied. The colour bar refers to t-statistics indicating the mean difference between connection strength estimates during strong and weak occipital alpha. | ![](assets/sub_21_D1_aIP_corrected.jpg)
Individual-level (subject 21) contrast imaging of the pyramidal to excitatory connection strength, α_pe, during times of strong and weak occipital alpha rhythm, with corrections for multiple comparisons applied. The colour bar refers to t-statistics indicating the mean difference between connection strength estimates during strong and weak occipital alpha. | ![](assets/sub_21_D1_aPE_corrected.jpg)
Individual-level (subject 21) contrast imaging of the pyramidal to inhibitory connection strength, α_pi, during times of strong and weak occipital alpha rhythm, with corrections for multiple comparisons applied. The colour bar refers to t-statistics indicating the mean difference between connection strength estimates during strong and weak occipital alpha. | ![](assets/sub_21_D1_aPI_corrected.jpg)
Individual-level (subject 21) contrast imaging of the external cortical input, μ, during times of strong and weak occipital alpha rhythm, with corrections for multiple comparisons applied. The colour bar refers to t-statistics indicating the mean difference between cortical input estimates under strong and weak posterior alpha. | ![](assets/sub_21_D1_input_corrected.jpg)
Individual-level (subject 21) contrast imaging of the membrane potential of excitatory population, V_e, during times of strong and weak occipital alpha rhythm, with corrections for multiple comparisons applied. The colour bar refers to t-statistics indicating the mean difference between membrane potential estimates under strong and weak posterior alpha. | ![](assets/sub_21_D1_v_es_corrected.jpg)
Individual-level (subject 21) contrast imaging of the membrane potential of inhibitory population, V_i, during times of strong and weak occipital alpha rhythm, with corrections for multiple comparisons applied. The colour bar refers to t-statistics indicating the mean difference between membrane potential estimates under strong and weak posterior alpha. | ![](assets/sub_21_D1_v_ii_corrected.jpg)
Individual-level (subject 21) contrast imaging of the membrane potential of pyramidal population, V_p, during times of strong and weak occipital alpha rhythm, with corrections for multiple comparisons applied. The colour bar refers to t-statistics indicating the mean difference between membrane potential estimates under strong and weak posterior alpha. | ![](assets/sub_21_D1_v_pyr_corrected.jpg)

#### Time efficiency
Here we show basic time statistics of running this framework using a MacBook Pro (Processor 2.9 GHz 6-Core Intel Core i9, Memory 32 GB 2400 MHz DDR4). The actual running time may vary depending on different computers/servers.

"AKF estimation": 3.16 hr on 6 cpus (14.5 seconds each channel)

"Contrast imaging": 1 hr on 6 cpus (with corrections for multiple comparisons); 0.14 hr on 6 cpus (without corrections for multiple comparisons)

## Common issues and solutions
Issue: *For Mac users, the error "invalid MEX-file" might happen at the first time of running.*

Solution: This link might help, https://in.mathworks.com/help/matlab/matlab_external/invalid-mex-file-error.html.

Issue: *MATLAB complains that mexmaci64 cannot be opened because the developer cannot be verified.*

Solution: This fieldtrip webpage might help, https://www.fieldtriptoolbox.org/faq/mexmaci64_cannot_be_opened_because_the_developer_cannot_be_verified/.

**Please send us an email at yun.zhao1@monash.edu for further help and support.**

## Work with your own data and visualisation packages
### Use your own data with AKF
NeuroProcImager works with EEG, ECoG, MEG source imaging data. You need to name your data file in this format 'data_[subject_index].mat', e.g., 'data_29.mat' means the data coming from subject 29, and put the data file in the directory '/data'. Then, you need to open the MATLAB script ‘test.m’ and put "AKF estimation" in the ‘pipeline’ list. You can run the script by entering the following command in the MATLAB Command Window: ‘test’. Finally, the output is a file containing neurophysiological variable estimates for every channel (source point) at ‘/output/variable_estimates_[subject_index].mat’, e.g., 'variable_estimates_29.mat' means the neurophysiological variable estimates for subject 29. There are eight matrices in the output file and they are the neurophysiological variable estimates for α_ep, α_ip, α_pe, α_pi, μ, v_e, v_i, v_p. Note that you can put more than 1 data file in the '/data' directory.

### Use your own visualisation packages
The output generated from "AKF estimation", the variable estimates file(s), can be further analysed by other statistical methods and visualised by other packages. In this repository, we provide one option, **Fieldtrip**, to visualise t-statistics derived from two-sample t-test which reveals the mean difference of variable estimates during strong and weak occipital alpha rhythm in a whole-brain manner.

To use other visualisation toolboxes or packages, one needs to note that in the variable estimates file there are eight matrices and they are the neurophysiological variable estimates for α_ep, α_ip, α_pe, α_pi, μ, v_e, v_i, v_p. Each matrix is of size (number of channels x timesteps of the recording).

## Code reusability
NeuroProcImager is an open-source software and all users can use it for research, teaching and learning purposes. We kindly ask you to email us if you want to use it in your publication or integrate it in your toolboxes.

Currently, we offer NeuroProcImager as a prototype containing two demo cases in terms of resting state alpha rhythm study: analytic Kalman Filter for neurophysiological variables estimation, whole-brain contrast imaging. In the upcoming updates, we will provide other estimation methods, other visualisation modules.
