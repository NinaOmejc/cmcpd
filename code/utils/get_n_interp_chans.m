path_main = 'D:\Experiments\corticomuscular_analysis\data\real\preproc\IZO\';

subjects = {"PDH04", "PDH05", "PDH06", "PDH07", "PDH09", "PDH10", "PDH12", "PDH14", "PDH18", "PDH19", "PDH20", "PDH22", "PDH23", "PDH24", "PDH25", "PDH26", ...
            "PDP02", "PDP03", "PDP04", "PDP06", "PDP09", "PDP10", "PDP11", "PDP13", "PDP15", "PDP17", "PDP18", "PDP19", "PDP20", "PDP21", "PDP22"}; %, ...

for isub = 1:length(subjects)
        sub = char(subjects{isub});
        path_sub = [path_main sub];
        load([path_sub filesep sub '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp.mat'])
    
        n_interp_chans(isub) = length(EEMG.etc.chans_interp_eeg);

    end
end

flat_ch_ratio_nonzeros = n_interp_chans(n_interp_chans ~= 0);
median(flat_ch_ratio_nonzeros)
max(flat_ch_ratio_nonzeros)