% Load each subjects data (of all 5 tasks) and run ICA on all of them. At
% the end, save each file (for each task) separately.

clear all;
path_data = 'D:\Experiments\corticomuscular_analysis\cmcpd\data\preproc\IZO';
config.subjects = ["PDH12", "PDP06"];
config.task_splits = {'SL', 'DL', 'SR', 'DR', 'C'};

for isubject = 1:length(config.subjects)
    
    subject = config.subjects{isubject};
    disp(['Subject ' subject])
    
    % Skip if already finished.
    if exist([path_data filesep subject filesep subject '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_splitDL_mc_ica.mat'], 'file')
        disp(['Skipping ' subject])
        continue    
    end
    
    % Try to load all files
    eeglab nogui
    successful_task_splits = [];
    isuccessful_task = 0;
    for itask = 1:length(config.task_splits)
        try
            task_split = config.task_splits{itask};
            load([path_data '\' subject '\' subject '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_mc_reref_norect_nocsd_split' task_split '_mc.mat']);
            EEMG = eeg_checkset(EEMG);
            isuccessful_task = isuccessful_task + 1;
            if isempty(ALLEEG); ALLEEG = EEMG; else; ALLEEG(isuccessful_task) = EEMG;end
            clear EEMG;
            successful_task_splits = [successful_task_splits; {task_split}];
        catch
            disp('Couldnt load file. Skipping.')
        end
    end
    n_succesful_files = length(successful_task_splits);    
    
    % 
    ALLEEG2 = pop_runica(ALLEEG, ...
                        'icatype', 'runica', ...
                        'extended',1, ...
                        'rndreset','yes', ...
                        'interrupt','on', ...
                        'pca',100, ...
                        'dataset', [1:n_succesful_files], ...
                        'concatenate', 'on', ...
                        'chanind', {'EEG'});

    for ieeg = 1:length(ALLEEG2)
        EEMG = ALLEEG2(ieeg);
        save([path_data '\' subject '\' subject '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_split'...
                                        successful_task_splits{ieeg} '_mc_ica.mat'], "EEMG")
        clear EEMG;
    end
    clear ALLEEG
end
