
subjects = {"PDH04", "PDH05", "PDH06", "PDH07", "PDH09", "PDH10", "PDH12", ...
           "PDH14", "PDH18", "PDH19", "PDH20", "PDH22", "PDH23", "PDH24", "PDH25", "PDH26",  ...
           "PDP02", "PDP03", "PDP04", "PDP06",  "PDP09", "PDP10", "PDP11", "PDP13", "PDP15", ...
           "PDP17", "PDP18", "PDP19", "PDP20", "PDP21", "PDP22"};
task = 'IZO';
task_split_names = {'SL', 'DL', 'SR', 'DR', 'C'};
sampling_freq = 300;

for isub = 1:length(subjects)
    subject = char(subjects{isub});
    for itasksplit = 1:length(task_split_names)
        
        task_split = task_split_names{itasksplit};
        path_preproc_data = ['D:\Experiments\corticomuscular_analysis\data\real\preproc\IZO\' subject];
        fname_preproc_data = [subject '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_split' task_split '_mc_icacl_icatdic.mat'];

        path_preproc_data_hpc = strrep(path_preproc_data, 'preproc', 'preproc_hpc');
        % path_preproc_data_hpc = ['D:\Experiments\corticomuscular_analysis\data\real\STUDY\' subject];
        if ~exist(path_preproc_data_hpc, 'dir'); mkdir(path_preproc_data_hpc); end
        try
            copyfile([path_preproc_data filesep fname_preproc_data], [path_preproc_data_hpc filesep fname_preproc_data])
            % load([path_preproc_data filesep fname_preproc_data])
            % EEG = EEMG; EEG = eeg_checkset(EEG); 
            % EEG = pop_saveset(EEG, [path_preproc_data_hpc filesep subject '_' task_split '.set']);
        catch
            disp(['Failed for ' fname_preproc_data])
        end
    end
end

