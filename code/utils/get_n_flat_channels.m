path_main = 'D:\Experiments\corticomuscular_analysis\data\real\preproc\IZO\';

subjects = {"PDH04", "PDH05", "PDH06", "PDH07", "PDH09", "PDH10", "PDH12", "PDH14", "PDH18", "PDH19", "PDH20", "PDH22", "PDH23", "PDH24", "PDH25", "PDH26", ...
            "PDP02", "PDP03", "PDP04", "PDP06", "PDP09", "PDP10", "PDP11", "PDP13", "PDP15", "PDP17", "PDP18", "PDP19", "PDP20", "PDP21", "PDP22"}; %, ...

flat_emgs = {};
flat_eegs = {};
n_flat_emgs = zeros(31, 1);

for isub = 1:length(subjects)
    sub = char(subjects{isub});
    path_sub = [path_main sub];
    load([path_sub filesep sub '_IZO_eemg_onlytask_fs500_dtrnd_lf1_hf150_interp.mat'])

    flat_emg = EEMG.etc.flat_emg_chans_removed;
    flat_eeg = EEMG.etc.flat_eeg_chans_interp;

    flat_emgs{isub} = flat_emg;
    flat_eegs{isub} = flat_eeg;
    n_flat_emgs(isub) = length(flat_emg);
end

flat_ch_ratio_nonzeros = n_flat_emgs(n_flat_emgs ~= 0);
median(flat_ch_ratio_nonzeros)
max(flat_ch_ratio_nonzeros)