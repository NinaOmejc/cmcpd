config_version = 'v1c_avgbrain';
load(['D:\Experiments\corticomuscular_analysis\configs\config_loc_' config_version '_dbi.mat'])
config.config_version = config_version;
eqsplit_txt = '';
config.fs = 50;
subject = 'PDH09';
config.task_split = 'SL';
config.remove_first_nwins = 5;
config.path_dbi_base = [config.path_root filesep 'data' filesep 'real' filesep 'dbi' filesep ...
    'IZO_' config.version '_' eqsplit_txt 'split' config.task_split '_dbiwin3'];

freq_band_eeg = 'eeg_beta';
freq_band_emg = 'emg_beta';
eeg_channel = 'avgcen';
emg_channel = 'emgl';

name_base = [freq_band_eeg '_' eeg_channel '_' freq_band_emg '_' emg_channel];
fname_dbi_cc = [subject '_IZO_eemg_dbi_dbiwin3_' name_base '.mat'];
fname_dbi_strength = [subject '_IZO_eemg_dbi_dbiwin3_' name_base '_cstrength_withsurr.mat'];

dbi_cc = load([config.path_dbi_base filesep emg_channel filesep fname_dbi_cc]);
dbi_strength = load([config.path_dbi_base filesep emg_channel filesep fname_dbi_strength]);

% plot coupling strength over time and average
tm = dbi_strength.tm;
tm_avgstart = tm(config.remove_first_nwins);
cs_time = dbi_strength.cpl_strength2;
surr_cs_time_m = mean(dbi_strength.surrs_cpl_strength2, 1);
surr_cs_time_std = std(dbi_strength.surrs_cpl_strength2);
surr_cs_time_m2std = surr_cs_time_m + 2*surr_cs_time_std;
cs_avg = mean(cs_time);
surr_cs_m2std = mean(surr_cs_time_m2std);
fig_title = [config.task_split ' |' subject '| ' strrep(name_base, '_', ' | ')];

figure('Visible', 'on'); hold on;
plot(tm, cs_time, 'k-', 'LineWidth',2);
plot(tm, surr_cs_time_m2std, 'b-', 'LineWidth',2)
yline(cs_avg, 'k--');
yline(surr_cs_m2std, 'c--');
xline(tm_avgstart, 'k-');
xlabel('Time'); ylabel('Coupling strength')
legend({'true', 'surr', 'true avg', 'surr avg'})
title(['CPL Strength | ' fig_title])





