% run by calling the function. Change the itask variable to complete all 5
% tasks.

function d_local_ica_single_task()
addpath 'D:\Experiments\corticomuscular_analysis\cmcpd\code\preprocessing'
eeglab nogui

path_data = 'D:\Experiments\corticomuscular_analysis\cmcpd\data\preproc\IZO';
config.subjects = ["PDH12", "PDP06"];
config.task_splits = {'SL', 'DL', 'SR', 'DR', 'C'};
itask = 1;

task_split = config.task_splits{itask};
successful_subs = [];
ALLEEG = [];
isuccessful_sub = 0;

for isubject = 1:length(config.subjects)
    subject = config.subjects{isubject};
    disp(['Subject ' subject])
    try
        load([path_data '\' subject '\' subject '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_split' task_split '_mc_icacl.mat']);
        EEMG = eeg_checkset(EEMG);
        isuccessful_sub = isuccessful_sub + 1;
        % average EMGs l and r
        config.cwt_chans_eeg = {'all'};
        config.cwt_chans_emg = {'average'};
        [chans_eeg, chans_emg] = get_chans(EEMG, config, 'cwt');
        [data, datanames] = get_data_of_channels(EEMG, chans_eeg, chans_emg, config);
        EEMG.data = [data.eeg; data.emg];
        cwt_chans = [datanames.eeg, datanames.emg];

        for icluster = 1: size(EEMG.data, 1)
            EEMG.chanlocs(icluster).labels = cwt_chans{icluster};
        end
        EEMG.chanlocs(size(EEMG.data, 1)+1:end) = [];
        EEMG = eeg_checkset(EEMG);

        if isempty(ALLEEG); ALLEEG = EEMG; else; ALLEEG(isuccessful_sub) = EEMG;end
        clear EEMG;
        successful_subs = [successful_subs; {subject}];
    catch
        disp('Couldnt load file. Skipping.')
    end
end
n_succesful_files = length(successful_subs);


ALLEEG2 = pop_runica(ALLEEG, ...
    'icatype', 'runica', ...
    'extended',1, ...
    'rndreset','yes', ...
    'interrupt','on', ...
    'pca',90, ...
    'dataset', [1:n_succesful_files], ...
    'concatenate', 'on', ...
    'chanind', {'EEG'});

for ieeg = 1:length(ALLEEG2)
    EEMG = ALLEEG2(ieeg);
    save([path_data '\' successful_subs{ieeg} '\' successful_subs{ieeg} '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_split'...
        task_split '_mc_icacl_icat.mat'], "EEMG")
    clear EEMG;
end
clear ALLEEG ALLEEG2
end


%%

function [chans_eeg, chans_emg] = get_chans(EEMG, config, analysis_type)

    n_eeg_chans = sum(strcmp({EEMG.chanlocs(:).type}, 'EEG'));

    if strcmp(analysis_type, 'cwt')
        chosen_chans_eeg = config.cwt_chans_eeg;
        chosen_chans_emg = config.cwt_chans_emg;
    elseif strcmp(analysis_type, 'pc')
        chosen_chans_eeg = config.pc_chans_eeg;
        chosen_chans_emg = config.pc_chans_emg;
    elseif strcmp(analysis_type, 'dbi')
        chosen_chans_eeg = config.dbi_chans_eeg;
        chosen_chans_emg = config.dbi_chans_emg;
    end

    % EEG
    if strcmp(chosen_chans_eeg{1}, 'all')
        chans_eeg = {EEMG.chanlocs(1:n_eeg_chans).labels}; 
    elseif strcmp(chosen_chans_eeg{1}, 'all_central_large')
        all_channels = {EEMG.chanlocs(1:n_eeg_chans).labels};
        central_channels_mask = contains(all_channels, 'C');
        chans_eeg = all_channels(central_channels_mask);
    elseif strcmp(chosen_chans_eeg{1}, 'central')
        chans_eeg = {'C1', 'C2', 'Cz', 'FCC1h', 'FCC2h', 'CCP1h', 'CCP2h'};
    elseif strcmp(config.dbi_chans_eeg{1}, 'avgcen')
        chans_eeg.avgcen = {'C1', 'C2', 'Cz', 'FCC1h', 'FCC2h', 'CCP1h', 'CCP2h'};
    elseif strcmp(chosen_chans_eeg{1}, 'avgbrain')
        chans_eeg.avgcen = {'C1', 'C2', 'Cz', 'FCC1h', 'FCC2h', 'CCP1h', 'CCP2h'};
        chans_eeg.avgpar = {'P1', 'P2', 'Pz', 'CPP1h', 'CPP2h', 'PPO1h', 'PPO2h'};
        chans_eeg.avgteml = {'FC5', 'C5', 'CP5', 'CCP5h', 'FCC5h'};
        chans_eeg.avgtemr = {'FC6', 'C6', 'CP6', 'CCP6h', 'FCC6h'};
        chans_eeg.avgocc = {'O1', 'O2', 'Oz', 'POO1', 'POO2', 'OI1h', 'OI2h'};
        chans_eeg.avgfro = {'F1', 'F2', 'Fz', 'AFF1h', 'AFF2h', 'FFC1h', 'FFC2h'};
    else
        chans_eeg = chosen_chans_eeg;
    end
    
    % EMG
    if strcmp(chosen_chans_emg{1}, 'all')
        chans_emg = {EEMG.chanlocs(n_eeg_chans+1:end).labels}; 
    elseif strcmp(chosen_chans_emg{1}, 'every_second')
        chans_emg = {EEMG.chanlocs(n_eeg_chans+1:2:end).labels}; 
    elseif strcmp(chosen_chans_emg{1}, 'every_third')
        chans_emg = {EEMG.chanlocs(n_eeg_chans+1:3:end).labels}; 
    elseif strcmp(chosen_chans_emg{1}, 'every_sixst')
        chans_emg = {EEMG.chanlocs(n_eeg_chans+1:6:end).labels}; 
    elseif strcmp(chosen_chans_emg{1}, 'average')
        emgr_indices = contains({EEMG.chanlocs(:).labels}, 'EMuoviI-');
        emgl_indices = contains({EEMG.chanlocs(:).labels}, 'EMuoviII-');
        chans_emg.emgr = {EEMG.chanlocs(emgr_indices).labels}; 
        chans_emg.emgl = {EEMG.chanlocs(emgl_indices).labels}; 
    else
        chans_emg = chosen_chans_emg;
    end

end


%%
function [data, datanames] = get_data_of_channels(EEMG, chans_eeg, chans_emg, config)
    
    % -----
    % EEG
    % -----

    % get EEG indices of the channels
    if isstruct(chans_eeg)
        names = fieldnames(chans_eeg);
        for ifield = 1:length(names)
            [~, eeg_channel_idx] = ismember(chans_eeg.(names{ifield}), {EEMG.chanlocs.labels});
            datatemp.(names{ifield}) = EEMG.data(eeg_channel_idx, :);
            data.eeg(ifield, :) = mean(datatemp.(names{ifield}), 1);
            datanames.eeg = names';
        end
    else
        [~, eeg_channel_idx] = ismember(chans_eeg, {EEMG.chanlocs.labels});
        data.eeg = EEMG.data(eeg_channel_idx, :);
        datanames.eeg = chans_eeg;
    end

    % -----
    % EMG
    % -----
    
    % get EMG indices of the channels and get data
    if isstruct(chans_emg)
        names = fieldnames(chans_emg);
        [~, emgr_channel_idx] = ismember(chans_emg.emgr, {EEMG.chanlocs.labels});
        [~, emgl_channel_idx] = ismember(chans_emg.emgl, {EEMG.chanlocs.labels});
        data.emg = [mean(EEMG.data(emgr_channel_idx, :), 1); ...
                    mean(EEMG.data(emgl_channel_idx, :), 1)];        
        datanames.emg = names';
    else
        [~, emg_channel_idx] = ismember(chans_emg, {EEMG.chanlocs.labels});
        data.emg = EEMG.data(emg_channel_idx, :);
        datanames.emg = chans_emg;
    end
end

