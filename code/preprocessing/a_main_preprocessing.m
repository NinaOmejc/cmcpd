clear all
close all
clc

%% 
% This is the main script to set the settings and from which the
% preprocessing script (eemg_preprocessing.m) is run. All settings are
% saved in config file.

% Needed plugins:
% eeglab, firfilt, dipfit, zapline-plus, ICLabel
% 


%%
% general settings
config.path_root = ['D:' filesep 'Experiments' filesep 'corticomuscular_analysis' filesep 'cmcpd'];
config.path_data_raw = [config.path_root filesep 'data' filesep 'raw'];
config.fullfname_eemg_channels =  [config.path_root filesep 'electrodes' filesep 'ChanLocs209.loc'];

config.n_eeg_channels = 128;
config.n_emg_channels = 64;
config.verbose = 1;             

% EEMG preprocessing settings
lets_do_eemg_preprocessing = 1;
config.force_recompute_preprocessing = 0;
 
config.subjects = ["PDP06", "PDH12"];
config.tasks = ["IZO"];  
config.plot_raw_data = 1;
config.resample_freq = 500;
config.do_detrend = 1;
config.remove_line_noise_with_zapline = 1;
config.low_freq_cutoff = 1;
config.high_freq_cutoff = 150;
config.band_freq = [];
config.do_channel_interpolation = 1;
config.do_auto_removal_of_break_segments = 1;
config.do_manual_removal_of_data_segments = 1;
config.rectify_emg = 0;
config.do_csd = 0;
config.do_reref = 1;
config.prepare_for_hpc_amica = 0;
config.load_amica_results = 1;
config.do_data_task_splits = 1;
config.make_equal_segments = 0;
config.manual_check_of_equal_segments = 1;
config.additional_resampling_before_task_split = 300;

config.do_manual_inspection_of_channel_interpolation = 1;
config.flatLineCriterion = 10;
config.channelCriterion = 0.85;
config.channelCriterionMaxBadTime = 0.5; 

config.amica_max_iter = 2000;                             % number amica iterations
config.reject_ica_components = 1;
config.threshold_for_eye_component_removal = 0.85;                           % threshold, if above -> reject eye related IC
config.threshold_for_other_than_eye_component_removal = 0.95;                  % threshold, if above -> reject muscle, line, heart or channel related ICs
config.save_ic_label_plots = 1;
config.do_manual_inspection_of_IC_removal = 1;
config.reject_components_during_preprocessing = 1;
config.save_brain_ICs = 1;


%%
if lets_do_eemg_preprocessing
    eemg_preprocessing(config)
    if config.verbose, disp('Preprocessing successfully ended.'); end
end
