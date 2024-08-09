% create config.mat for connectivity part (after preprocessing steps)

function create_config()

config.version = 'test';

config.do_cwt = 1;  % continous wavelet transform
config.do_pc = 1;   % phase coherence
config.do_dbi = 1;  % dynamic bayesian inference

%% general settings
% config.path_root = '/d/hpc/home/nomejc/corticomuscular_analysis';
% config.file_to_ic_clusters = 'D:\Experiments\corticomuscular_analysis\cmcpd\data\STUDY\cluster_info.mat';
config.path_root = ['D:' filesep 'Experiments' filesep 'corticomuscular_analysis' filesep 'cmcpd'];
config.n_eeg_channels = 128;
config.n_emg_channels = 64;
config.make_plot_visible = 'off';
config.verbose = 1;
config.save_to_dbx = 0;
config.dbx_access_token = '';
config.save_figures = 0;

%% EEMG preprocessing settings
%config.subjects = {"PDH04", "PDH05", "PDH06", "PDH07", "PDH09", "PDH10", "PDH12", "PDH14", "PDH18", "PDH19", "PDH20", "PDH22", "PDH23", "PDH24", "PDH25", "PDH26", ...
%                   "PDP02", "PDP03", "PDP04", "PDP06", "PDP09", "PDP10", "PDP11", "PDP13", "PDP15", "PDP17", "PDP18", "PDP19", "PDP20", "PDP21", "PDP22"}; %, ...
config.subjects = {"PDH12", "PDP06"};
config.subjects_removed = {"PDH08"};
config.tasks = ["IZO"];
config.with_resample_freq = 300;
config.with_removed_line_noise_with_zapline = 1;
config.with_manual_removal_of_data_segments = 1;
config.with_auto_removal_of_break_segments = 1;
config.with_low_freq_cutoffs = [1];
config.with_high_freq_cutoffs = [150];
config.with_band_freq = [];
config.with_detrend = 1;
config.with_chan_interp = 1;
config.with_rerefs = [1];
config.with_rectify_emgs = [0];
config.with_csds = [0];
config.with_amica_clean = 1;
config.with_brain_icss = [0];
config.task_splits = {'SL', 'DL', 'SR', 'DR', 'C'};
config.with_equal_segments = 0;

%% Time freq decomposition settings
config.cwt_freq_min = 4;
config.cwt_freq_max = 90;
config.cwt_f0 = 1;
config.cwt_nv = 30;
config.cwt_chans_eeg = {'all'}; % can also go 'all', 'all_central' or 'avgbrain' or 'ics', or particular channels {'Cz', 'C1'} (see function get_chans in utils folder)
config.cwt_chans_emg = {'average'}; % can also go 'every_second' or 'every_fourth'. Average is calculation of average for L and R (emgl and emgr)
config.cwt_save_power_insteadof_complex = 1; % so its appropriate to downsample cause its real number
config.cwt_reduce_size_before_save = 1;
config.cwt_ds_factor = 6; % 300/50 (initial/new fs) downsampling factor

%% Phase coherence settings
config.pc_ns = [1]; % integer ratios for phases 
config.pc_ms = [1]; %2 3 4 5]; % for eeg
config.pc_save_single_emg_results = 0;
config.pc_add_surrogates = 0;
config.pc_create_surrogates_on_spot = 0; % to create surrogates of different sizes
config.pc_n_surrogates = 3;           % 30 time series -> at the end ( (n*n-1) / 2 = 435 coherences)
config.pc_npnts_surr = 20000;         % config.n_pnts_surr = 35341 for srate = 500 Hz. Just for the fixed surrogates (first versions, not important now)
config.pc_chans_eeg = {'all'};        % can also go 'all', 'all_central' or 'avgbrain' or particular channels {'Cz', 'C1'} (see function get_chans in utils folder)
config.pc_chans_emg = {'average'};    % can also go 'every_second' or 'every_fourth'
config.pc_reduce_size_before_save = 1;
config.pc_ds_factor = 6; % 300/50 (initial/new fs) downsampling factor
config.pc_allow_to_load_cwt = 0;

%% Dynamic Bayesian inference settings
config.dbi_freq_bands = struct('alpha', [7 12], 'beta_low', [12.5 20], 'beta_high', [20.5 30], 'beta', [12.5 30], 'gamma_low', [30.5 48], 'gamma_high', [52 90]); % config.dbi_freq_bands = struct('mu', [11 15]);         % 
config.dbi_win_sizes = [3];      % in seconds
config.dbi_winoverlap = 0.5;     % in procentages
config.dbi_chans_eeg = {'Cz'}; % can also go 'all', 'all_central' or 'avgbrain' or particular channels {'Cz', 'C1'} (see function get_chans in utils folder)
config.dbi_chans_emg = {'average'}; % can also go 'all', 'every_second' or 'every_fourth'
config.dbi_compute_n_surrogates = 0;
config.dbi_save_func_over_time_plots = 0;

%% exceptions
config.exceptions_to_remove = struct('PDH05', 'emgr', ...
                                     'PDP09', 'emgr');

config.force_recompute_cwt = 0;
config.force_recompute_pc = 0;
config.force_recompute_dbi = 0;

% save config file 
if config.do_cwt; analysis_txt = '_cwt'; else; analysis_txt = ''; end
if config.do_pc; analysis_txt = [analysis_txt '_pc']; else;  analysis_txt = [analysis_txt '']; end
if config.do_dbi; analysis_txt = [analysis_txt '_dbi']; else;  analysis_txt = [analysis_txt '']; end
if contains(config.path_root, '\'); environment = 'loc'; else; environment = 'hpc'; end

if ~exist([config.path_root filesep 'configs'], 'dir'); mkdir([config.path_root filesep 'configs']); end
fullfname_config = [config.path_root filesep 'configs' filesep 'config_' environment '_' config.version analysis_txt '.mat'];
save(fullfname_config, "config");
disp(fullfname_config)

end