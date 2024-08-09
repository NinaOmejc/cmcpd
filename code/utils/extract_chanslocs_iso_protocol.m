%% extract protocol time series and make a figure
load('D:\Experiments\corticomuscular_analysis\data\real\raw\PDH12_IZO_EMG_EEG_realigned.mat');
force = EEG_realigned.data(209, :);
protocol = EEG_realigned.data(208, :);

fs_new = 500;
ds_factor = 2000/fs_new;
force_ds = downsample(force, ds_factor);
protocol_ds = downsample(protocol, ds_factor);

x1 = 93000; x2 = 110000;
x3 = 127000; x4 = 143500;
protocol_to_plot = protocol_ds([x1:x2 x3:x4]);
force_to_plot = force_ds([x1:x2 x3:x4]);

t = 1/fs_new : 1/fs_new : length(protocol_to_plot)/fs_new;
figure, plot(t, protocol_to_plot, 'k-', 'LineWidth', 2); hold on; plot(t, force_to_plot, 'r-');

save(['D:\Experiments\corticomuscular_analysis\data\protocol_ISO_timeseries.mat'], "protocol_to_plot", "force_to_plot", "t");


%% extract chanlocs
load('D:\Experiments\corticomuscular_analysis\data\real\preproc\IZO\PDH05\PDH05_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_splitC_mc_icacl_icat.mat')
eeglab nogui; EEMG = pop_select(EEMG, 'chantype', 'EEG');
chanlocs = EEMG.chanlocs;
save(['D:\Experiments\corticomuscular_analysis\data\chanlocs.mat'], "chanlocs");
