path_main = 'D:\Experiments\corticomuscular_analysis\data\real\preproc\IZO\';

subjects = {"PDH04", "PDH05", "PDH06", "PDH07", "PDH09", "PDH10", "PDH12", "PDH14", "PDH18", "PDH19", "PDH20", "PDH22", "PDH23", "PDH24", "PDH25", "PDH26", ...
            "PDP02", "PDP03", "PDP04", "PDP06", "PDP09", "PDP10", "PDP11", "PDP13", "PDP15", "PDP17", "PDP18", "PDP19", "PDP20", "PDP21", "PDP22"}; %, ...
tasks = {'SL', 'SR', 'DL', 'DR', 'C'};
i = 1;
for isub = 1:length(subjects)
    for itask = 1:length(tasks)

        sub = char(subjects{isub});
        task = char(tasks{itask});
        path_sub = [path_main sub];
        try
            load([path_sub filesep sub '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_split' task '_mc_icacl_icatdic.mat'])
        
            n_comps_rej(i) = size(EEMG.etc.spatial_filter.rejected_amica_comps, 1);
            n_comps_brain(i) = length(EEMG.etc.spatial_filter.brain_comps);
            i = i+1;
        catch
            disp(['Doesnt exist: ' sub ' ' task])
        end
    end
end

mean(n_comps_brain)
std(n_comps_brain)