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
    elseif strcmp(analysis_type, 'ac')
        chosen_chans_eeg = config.ac_chans_eeg;
        chosen_chans_emg = config.ac_chans_emg;
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
    elseif strcmp(chosen_chans_eeg{1}, 'ic')
        % eeg_brain_ics = eeg_getdatact(EEMG, 'component', 1:size(EEMG.icaweights, 1));
        if ~isfield(config, 'file_to_ic_clusters') || isempty(config.file_to_ic_clusters)
            chans_eeg = cellstr(strcat('IC', arrayfun(@num2str, 1:size(EEMG.icaweights, 1), 'UniformOutput', false)));
        else
            load(config.file_to_ic_clusters)
            chans_eeg = fieldnames(cluster_info)';
        end
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
        if ismember('emgr', {EEMG.chanlocs(:).labels})
            emgr_indices = contains({EEMG.chanlocs(:).labels}, 'emgr');
            emgl_indices = contains({EEMG.chanlocs(:).labels}, 'emgl');
            chans_emg.emgr = {EEMG.chanlocs(emgr_indices).labels}; 
            chans_emg.emgl = {EEMG.chanlocs(emgl_indices).labels};            
        else
            emgr_indices = contains({EEMG.chanlocs(:).labels}, 'EMuoviI-');
            emgl_indices = contains({EEMG.chanlocs(:).labels}, 'EMuoviII-');
            chans_emg.emgr = {EEMG.chanlocs(emgr_indices).labels}; 
            chans_emg.emgl = {EEMG.chanlocs(emgl_indices).labels}; 
        end
    else
        chans_emg = chosen_chans_emg;
    end

end