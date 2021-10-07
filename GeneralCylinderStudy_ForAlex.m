%% General cylinders study

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was created to obtain R2* synthetic data for any cylinder
% configuration and be fitted afterwards by R2* signal models.
% The code operates as follows: A set of
% directions (defined by a file or created on fly) are created, 
% each one representing a cylinder. Each one contributes to the
% total signal by its intra-axonal and myelin (if not neglected)
% compartments, and its effect to the extra-axonal compartment.
% All variables are recommended to be defined at the beginning of
% the code. When the signal is created, it is fitted immediately 
% by each model (classic and quadratics).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Directions (and checking if they are normalised):
path_for_data = fileparts(pwd);

Data = load('Directions1500_Cartesian.mat'); % If defined by text file
DAxons1 = Data.DAxons1;
DAxons1_Norm = [DAxons1(:,1), DAxons1(:,2), DAxons1(:,3)]./sqrt(sum(DAxons1.^2,2));

% Angle sampling and time sampling.
angle_values = pi/90:pi/90:pi/2; % In radians. 

kappa_values = [0.1:0.1:6.0,7.5,10,20]; % kappa values undersampled
                              % for SNR - noise
                              % study                                                     

g_samples = 0.5:0.01:0.9; % Here it includes all the experimental ones used for the paper
FVF_samples = 0.5:0.1:0.9;

%pc_information = memory;

%if pc_information.MemUsedMATLAB > ((5000*63*241*45)*8) + 100
    % In case that memory is large enough, run this:
    time_range = (0:0.00025:0.06)*1000; % Oversampling time in msec. (time sampling
                                       % reduced from 0.1 ms to 0.25 ms and max
                                       % of 60 ms for SNR - noise study)
    samples = 5000;
%elseif pc_information.MemUsedMATLAB > ((2500*63*121*45)*8) + 100
    % Otherwise, execute instead these values:
%    time_range = (0:0.0005:0.06)*1000; % Oversampling time in msec. (time sampling
                                       % reduced from 0.1 ms to 0.25 ms and max
                                       % of 60 ms for SNR - noise study)
%    samples = 2500;
%else
%    error('No possible configuration available to run this code. Change accordingly');
%end

max_time = [18, 27, 36, 45, 54];
exp_time_range = 3.4:3.34:53.5; % Experimental undersampled time


% Variables:
DispNoMye_Sampling_points = [];
NoDispNoMye_Sampling_points = [];
DispMye_Sampling_points = [];
NoDispMye_Sampling_points = [];
datatypes = {'DispNoMye','NoDispNoMye','DispMye','NoDispMye'};

% Auxiliar variables:
% NOTE: For data synthesis (next section), it is advisible to use a large
% ranges of g-ratios and FVF values (from 1:5 in both cases). However, for
% data_fitting (the next-next section), please perform it by using only ONE
% VALUE PER G-RATIO/FVF! Otherwise it would produce memory issues. Sorry
% for the inconvinience.
%g_indx = 3;
%FVF_indx = 1;
datatype_indexes = [1,3]; % Change HERE the indexes used for data analysis.
% The order, once more, goes as follows:
% datatype_indexes = 1 -> No myelin with dispersion
% datatype_indexes = 2 -> No myelin with no dispersion
% datatype_indexes = 3 -> Myelin with dispersion
% datatype_indexes = 4 -> Myelin with no dispersion

% Extra parameters useful for renaming the file name accordingly
%suffix_string = {'_T2ExVivo_GratioExVivo';'_T2ExVivo_GratioWB';'_T2WB_GratioExVivo';'_T2ExVivo_GratioMiddle'};
%suffix_string = {'_T2WB_GratioMiddle'};
suffix_string = {'_T2ExVivo_GratioExVivo';'_T2ExVivo_GratioMiddle';...
                 '_T2ExVivo_GratioWB';'_T2WB_GratioExVivo';...
                 '_T2WB_GratioMiddle';'_T2WB_GratioWB'};
g_samples_for_study = [0.66,0.73,0.8];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Creating synthetic R2* signal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FJv21(28.01): Ex vivo setup (Note: find the references for using these
% numbers. 
% FJv21(11.02): Ok, the "ex vivo setup" was performed with the following
% parameters: T2m ~ 13.26 ms and T2a = T2e ~ 53.96 ms from (https://onlinelibrary.
%wiley.com/doi/epdf/10.1002/mrm.22267) which, doing a monoexponential study
%on it, the "fitted" mean T2 is similar of what Evgeniya reported (T2_mean ~
%45 ms, from https://advances.sciencemag.org/content/6/41/eaaz9281), FVF and G-
%ratio were from Laurin's presentation, resulting in FVF ~ (0.213 + 0.276) = 0.49
% and G-ratio ~ sqrt(1 - 0.276/0.49) ~ 0.66. For this study, just with and without
% exchange studies were performed.
%params.E = 0; % Do for exchange and not exchange.
%params.FVF = 0.7;
% FJv21(16.08): Due to space and consistency, only the R2* signal decay at
% SNR = infty will be saved. Noise addition and such is "added on fly" in
% data fitting.

R2a_test = [1/53.96,1/36];
R2e_test = [1/53.96,1/36];
R2m_test = [1/13.26,1/8];

%This is overloaded just to avoid hard code modification, and mostly for
%".mat file name" and "loop".
FVF_indx = 1;
for R2_indx = 1:numel(R2a_test)
    params.R2m = R2m_test(R2_indx);
    params.R2a = R2a_test(R2_indx);
    params.R2e = R2e_test(R2_indx);
    h = waitbar(0,'Initializing waitbar...');
    
    for g_indx = 1:numel(g_samples_for_study)     
        g_sample = g_samples(g_samples == g_samples_for_study(g_indx));
        for k = 1:numel(kappa_values)  
            for j = 1:numel(angle_values)
                params.theta = angle_values(j);
                params.kappa = kappa_values(k);
                params.FiberLimit = 0.5;
                % Getting equations for datapoints
                Signal_time = SignalModelR2(DAxons1_Norm,...
                                            FVF_samples(FVF_indx),...
                                            g_sample,...
                                            params);

                % Signal at SNR = infty is in complex, so complex noise is
                % added afterwards
                DispNoMye_Sampling_points(k,j,1:numel(time_range)) = Signal_time{1}(time_range);
                NoDispNoMye_Sampling_points(k,j,1:numel(time_range)) = Signal_time{3}(time_range);
                DispMye_Sampling_points(k,j,1:numel(time_range)) = Signal_time{2}(time_range);
                NoDispMye_Sampling_points(k,j,1:numel(time_range)) = Signal_time{4}(time_range);
                
                waitbar((j + numel(angle_values)*(k-1) + numel(angle_values)*numel(kappa_values)*(g_indx-1))...
                    /(numel(g_samples_for_study)*numel(kappa_values)*numel(angle_values)),h,...
                    ['FVF-Gratio: ' num2str(FVF_samples(FVF_indx)) '-' num2str(g_sample) ', with \kappa = ' num2str(kappa_values(k)) ', \theta_\mu ' num2str(angle_values(j)*180/pi)]);
            end
        end
        
        if ~exist(fullfile(path_for_data,'InSilico_SignalDecay'),'dir')
            mkdir(path_for_data,'InSilico_SignalDecay');
        end
        
        % data saving per signal equation and SNR:
        save(fullfile(path_for_data,'InSilico_SignalDecay',['DispNoMye_data_Gratio' num2str(g_sample*100) '_FVF' num2str(FVF_samples(FVF_indx)*10) suffix_string{g_indx + (R2_indx-1)*numel(g_samples_for_study)} '.mat']), 'DispNoMye_Sampling_points','-v7.3');
        save(fullfile(path_for_data,'InSilico_SignalDecay',['NoDispNoMye_data_Gratio' num2str(g_sample*100) '_FVF' num2str(FVF_samples(FVF_indx)*10) suffix_string{g_indx + (R2_indx-1)*numel(g_samples_for_study)} '.mat']), 'NoDispNoMye_Sampling_points','-v7.3');   
        save(fullfile(path_for_data,'InSilico_SignalDecay',['DispMye_data_Gratio' num2str(g_sample*100) '_FVF' num2str(FVF_samples(FVF_indx)*10) suffix_string{g_indx + (R2_indx-1)*numel(g_samples_for_study)} '.mat']), 'DispMye_Sampling_points','-v7.3');   
        save(fullfile(path_for_data,'InSilico_SignalDecay',['NoDispMye_data_Gratio' num2str(g_sample*100) '_FVF' num2str(FVF_samples(FVF_indx)*10) suffix_string{g_indx + (R2_indx-1)*numel(g_samples_for_study)} '.mat']), 'NoDispMye_Sampling_points','-v7.3');    
        
        clear 'DispNoMye_Sampling_points' 'NoDispNoMye_Sampling_points' 'DispMye_Sampling_points' 'NoDispMye_Sampling_points'               
    end
close(h);
end
%
