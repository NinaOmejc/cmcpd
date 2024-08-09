function inner_eemg_preprocessing(config)
    
    config = get_combinations(config);

    for idx_job = 1:size(config.combinations, 1)
        
        [config, to_skip_all, start_after_interp, start_after_amica] = set_up_job(config, idx_job);
        if to_skip_all, continue, end
        eeglab nogui
        
        if ~start_after_amica
            if ~start_after_interp
                EEMG = get_data(config);            
                [EEMG, config] = get_channels_info(EEMG, config); % Global_1.
                EEMG = resample_data(EEMG, config); 
                EEMG = check_events(EEMG, config);
                EEMG = remove_flat_line_chans(EEMG, config);
                EEMG = auto_remove_break_segments(EEMG, config);
                EEMG = detrend_data(EEMG, config);
                EEMG = filter_data(EEMG, config);  
                EEMG = channel_interpolation_wrapper(EEMG, config);
            else
                disp('Loading file with interpolation finished.')
                load([config.path_data_preproc filesep config.fname_data_preproc_after_interp]);
                EEMG = eeg_checkset(EEMG);
            end
            EEMG = manually_remove_data_segments(EEMG, 10, 150, 32, config); 
            EEMG = reref(EEMG, config);
            EEMG = rectify_emg(EEMG, config);
            EEMG = do_csd(EEMG, config);
            prepare_for_amica(EEMG, config)
            split_data_in_tasks(EEMG, config);
        else
            [EEMG, EEG] = load_data_and_amica_results(config);
            [EEG, eeg_brain_ics, brain_comps] = decide_on_components(EEG, config);
            EEMG.data(1:config.n_eeg_channels, :) = EEG.data;
            save_data(EEMG, config, 'after_amica_clean')
            split_data_in_tasks(EEMG, config)
            save_EEMG_with_brain_comps(EEMG, eeg_brain_ics, brain_comps, config)
            split_data_in_tasks(EEMG, config)
        end
    end
end


%%

function config = get_combinations(config)
    
    % Create grids for task_types and subjects and Reshape the grids to get all combinations
    [taskGrid, subjectGrid] = ndgrid(config.tasks, config.subjects);
    config.combinations = [taskGrid(:), subjectGrid(:)];

end



%%

function [config, to_skip_all, start_after_interp, start_after_amica] = set_up_job(config, idx_job)
    
    config.task = char(config.combinations{idx_job, 1});
    config.subject = char(config.combinations{idx_job, 2});
    if config.verbose, fprintf('\nWorking on Job idx %d: | Task Type: %s | Subject: %s\n', idx_job, config.task, config.subject); end

    config = set_paths(config);
    [to_skip_all, start_after_interp, start_after_amica] = check_if_already_done(config);

end


%% 

function config = set_paths(config)
    
    path_data_split = strsplit(config.path_data_raw, filesep);
    config.path_data_base = strjoin(path_data_split(1:end-1), filesep);
    config.path_data_preproc = [config.path_data_base filesep 'preproc' filesep char(config.task) filesep char(config.subject)];
    config.path_local_amica = [config.path_data_preproc filesep 'amica_results_' config.task '_' config.subject filesep];
    config.path_local_for_hpc_amica = [config.path_data_base filesep 'preproc' filesep 'amica_on_hpc'];
    config.path_amica_on_hpc = 'amica_on_hpc/';
    config.fname_data_raw = sprintf('%s_%s_EMG_EEG_realigned.mat', config.subject, config.task);    
    
    if config.do_detrend; detrend_text = 'dtrnd'; else; detrend_text = 'nodtrnd'; end
    if config.low_freq_cutoff ~= 0; lowfiltstext = ['lf' num2str(config.low_freq_cutoff)]; else lowfiltstext = 'lf0'; end
    if config.high_freq_cutoff ~= 0; highfiltstext = ['hf' num2str(config.high_freq_cutoff)]; else highfiltstext = 'hf0'; end
    if config.do_manual_removal_of_data_segments; mclean_text = 'mc'; else; mclean_text = 'nomc'; end
    if config.do_reref ~= 0; rereftext = 'reref'; else; rereftext = 'noreref'; end
    if config.rectify_emg; rect_text = 'rect'; else; rect_text = 'norect'; end
    if config.do_csd ~= 0; csdtext = 'csd'; else; csdtext = 'nocsd';end
    if config.load_amica_results ~= 0; amica_txt= 'amica'; else; amica_txt = 'noamica';end
    if config.reject_ica_components ~=0; icaclean_txt = 'compclean'; else; icaclean_txt = 'nocompclean'; end 
    if config.do_auto_removal_of_break_segments; data_type = 'onlytask'; else; data_type = 'allduration'; end
    if config.do_data_task_splits; split_txt = '_split'; else; split_txt = ''; end
    if config.save_brain_ICs; brainicstxt = '_amica_brain_ICs'; else; brainicstxt = ''; end
    if config.resample_freq ~= 0; fstext = ['fs' num2str(config.resample_freq)]; else; fstext = ['fs2000']; end
    
    config.fname_data_preproc_after_interp =  [config.subject '_' config.task '_eemg_' data_type '_' fstext '_' detrend_text '_' lowfiltstext '_' highfiltstext '_interp'];
    config.fname_data_preproc_after_csd = [config.fname_data_preproc_after_interp '_' mclean_text '_' rereftext '_' rect_text '_' csdtext];
    config.fname_data_preproc_after_amica_clean = [config.fname_data_preproc_after_csd '_' amica_txt '_' icaclean_txt];
    config.fname_data_preproc_after_brain_ics = [config.fname_data_preproc_after_csd brainicstxt];
    config.fname_data_preproc_after_split = strrep(config.fname_data_preproc_after_interp, fstext, ['fs' num2str(config.additional_resampling_before_task_split)]);

    config.path_local_for_hpc_amica_outdir = [config.path_local_for_hpc_amica filesep config.fname_data_preproc_after_csd '_amicaout'];

    if ~exist(config.path_data_preproc, 'dir')
        mkdir(config.path_data_preproc)
    end

    if ~exist(config.path_local_for_hpc_amica, 'dir') && config.prepare_for_hpc_amica
        mkdir(config.path_local_for_hpc_amica)
    end

    % correction
    if exist([config.path_data_preproc filesep strrep(config.fname_data_preproc_after_interp, ['_fs' num2str(config.resample_freq)], ''), '.mat'], 'file')
        config.fname_data_preproc_after_interp = strrep(config.fname_data_preproc_after_interp, ['_fs' num2str(config.resample_freq)], '');
    end

end

%%

function [to_skip, start_after_interp, start_after_amica, do_final_check] = check_if_already_done(config)

    % if there is a file with finished preprocessing until amica
    if exist([config.path_data_preproc filesep config.fname_data_preproc_after_split '_splitC_mc.mat'], 'file')
        if config.force_recompute_preprocessing 
            to_skip = 0;
            if config.verbose, disp(['It seems csd has already been done for this job configuration. BUT FORCED TO RECOMPUTE, continuing...' newline]); end
        else
            to_skip = 1;
            if config.verbose, disp(['It seems csd has already been done for this job configuration. Skipping...']); end
        end
    else
        to_skip = 0;
    end

    % if there is a file with finished interpolation
    if exist([config.path_data_preproc filesep config.fname_data_preproc_after_interp '.mat'], 'file') || exist([config.path_data_preproc filesep strrep(config.fname_data_preproc_after_interp, ['_fs' num2str(config.resample_freq)], ''), '.mat'], 'file')
        if config.force_recompute_preprocessing 
            start_after_interp = 0;
            if config.verbose, disp(['It seems interp has already been done for this job configuration. BUT FORCED TO RECOMPUTE, continuing...']); end
        else
            start_after_interp = 1;
            if config.verbose, disp(['It seems interpolation has already been done for this job configuration. Skipping...']); end
        end
    else
        start_after_interp = 0;
    end
    
    % check if amica is finished
    if exist([config.path_local_amica], 'dir') || exist(config.path_local_for_hpc_amica_outdir, 'dir')
        if config.force_recompute_preprocessing
            start_after_amica = 0;
            if config.verbose, disp([newline 'It seems amica has already been done for this job configuration. BUT FORCED TO RECOMPUTE, continuing...' newline]); end
        else
            start_after_amica = 1;
            if config.verbose, disp([newline 'It seems amica has already been done for this job configuration. Skipping...' newline]); end
        end
    else
        start_after_amica = 0;
    end

end

%%

function EEMG = get_data(config)
    
    % load data
    if config.verbose,  disp('Loading raw files...');   end
    loaded_variable = whos('-file', [config.path_data_raw filesep config.fname_data_raw]);
    load([config.path_data_raw filesep config.fname_data_raw])
    EEMG = eval(loaded_variable.name);

    EEMG.setname = [config.subject '_' config.task '_EMG_EEG_realigned'];
    EEMG.filename = [config.subject '_' config.task '_EMG_EEG_realigned.set'];
    EEMG.filepath = config.path_data_raw;
    EEMG.subject = config.subject;
    EEMG.group = config.subject(3);
    EEMG.condition = config.task;
    EEMG = eeg_checkset(EEMG);

    % check for NaNs
    [~, col] = find(isnan(EEMG.data));
    col = unique(col);
    if ~isempty(col)
        if config.verbose, warning('NaNs values were detected in data. Removing that part.');  end
        EEMG = pop_select(EEMG, 'rmpoint', [min(col),max(col)]);
        EEMG = eeg_checkset(EEMG);
    end

    if config.plot_raw_data
        spacing = 800; dispchans = 64; wait_to_close = 1;
        title = ['subject ' config.subject ' | ' config.task '. Raw data.'];
      
        plot_data(EEMG, spacing, EEMG.xmax, dispchans, title, wait_to_close)
    end
end


%%

function [EEMG, config] = get_channels_info(EEMG, config)
    
    EEMG = pop_chanedit(EEMG, 'load', config.fullfname_eemg_channels, 'filetype', 'autodetect');
    
    for idx_channel = 1:EEMG.nbchan

        if idx_channel <= config.n_eeg_channels
            EEMG.chanlocs(idx_channel).type = 'EEG';
            EEMG.chanlocs(idx_channel).chantype = 'EEG';
        elseif regexp(EEMG.chanlocs(idx_channel).labels, 'EMuovi.*\d$')
            EEMG.chanlocs(idx_channel).type = 'EMG';
            EEMG.chanlocs(idx_channel).chantype = 'EMG';
        elseif strcmp(EEMG.chanlocs(idx_channel).labels, 'TRIGGER')
            EEMG.chanlocs(idx_channel).type = 'TRIG_EEG';
            EEMG.chanlocs(idx_channel).chantype = 'TRIG_EEG';         
        elseif strcmp(EEMG.chanlocs(idx_channel).labels,'EStation-aux-2')
            EEMG.chanlocs(idx_channel).type = 'TRIG_EMG';
            EEMG.chanlocs(idx_channel).chantype = 'TRIG_EMG';        
        else
            EEMG.chanlocs(idx_channel).type = 'OTH';
            EEMG.chanlocs(idx_channel).chantype = 'OTH';  
        end
    end
    
    EEMG = pop_select(EEMG, 'chantype', {'EEG', 'EMG'});
    EEMG = eeg_checkset(EEMG);

    if config.verbose disp(['EEG channel information successfully found ... ' ...
            'Auxillary channels were removed. Only type EEG and EMG were kept.']); end
    
end

%%
function EEMG = resample_data(EEMG, config)
    if config.resample_freq ~= 0 || config.resample_freq ~= EEMG.srate
        EEMG = pop_resample(EEMG, config.resample_freq);
        disp('Resampling')
    else
        disp('Skipping resampling.')
    end
end

%% 

function EEMG = detrend_data(EEMG, config)
    if config.do_detrend
        disp('Detrending...')

        for ichan = 1:EEMG.nbchan
            [tsd,~,~] = inner_detrend_flow(EEMG.data(ichan, :), EEMG.srate, 10);
            EEMG.data(ichan, :) = tsd;
        end
    else
        disp('skipping detrending')
    end
end

%%

function EEMG = remove_flat_line_chans(EEMG, config)
    flat_eeg_chans_idx = [];
    flat_emg_chans_idx = [];
    for ichan = 1:EEMG.nbchan
        if ~any(EEMG.data(ichan, :))
            if ichan > config.n_eeg_channels
                flat_emg_chans_idx = [flat_emg_chans_idx; ichan];
            else
                flat_eeg_chans_idx = [flat_eeg_chans_idx; ichan];
            end
        end
    end
    
    disp(['flat line EMG chans: ' char(strjoin(string(flat_emg_chans_idx), ' ')) '. Removing.']);
    EEMG = pop_select(EEMG, 'rmchannel', flat_emg_chans_idx);
    EEMG.etc.flat_emg_chans_removed = flat_emg_chans_idx;

    disp(['flat line EEG chans: ' char(strjoin(string(flat_eeg_chans_idx), ' ')) '. Interpolating.']);
    if ~isempty(flat_eeg_chans_idx)
        EEMG = pop_interp(EEMG, flat_eeg_chans_idx);
    end
    EEMG.etc.flat_eeg_chans_interp = flat_eeg_chans_idx;
    EEMG = eeg_checkset(EEMG);
end

%% 
function [EEMG] = check_events(EEMG, config)
        
    % check if types are correcty named
    if any(contains({EEMG.event(:).type}, 'condition'))
        for ievent = 1:length({EEMG.event(:).type})
            if contains(EEMG.event(ievent).type, 'condition')
                old_type = EEMG.event(ievent).type;
                new_type = strrep(old_type, 'condition ', '');
                EEMG = pop_editeventvals(EEMG, 'changefield', {ievent 'type' new_type});
            end
        end
    end

    % rename 0 and 32768 event codes to 9 and 10
    events_without_boundaries = EEMG.event(~strcmp({EEMG.event.type}, 'boundary'));
    event_STsub_latencies_start = [];
    event_STsub_latencies_end = [];
    for ie = 8:length(events_without_boundaries)
        if strcmp(events_without_boundaries(ie-7).type, '32768') && ...
                strcmp(events_without_boundaries(ie-6).type, '0') && ...
                strcmp(events_without_boundaries(ie-5).type, '32768') && ...
                strcmp(events_without_boundaries(ie-4).type, '0') && ...
                strcmp(events_without_boundaries(ie-3).type, '32768') && ...
                strcmp(events_without_boundaries(ie-2).type, '0') && ...
                strcmp(events_without_boundaries(ie-1).type, '32768') && ...
                strcmp(events_without_boundaries(ie).type, '0') && ...
                events_without_boundaries(ie-4).latency - events_without_boundaries(ie-7).latency < 2*EEMG.srate && ...
                events_without_boundaries(ie-4).latency - events_without_boundaries(ie-6).latency < 2*EEMG.srate && ...
                events_without_boundaries(ie-4).latency - events_without_boundaries(ie-5).latency < 2*EEMG.srate && ...
                events_without_boundaries(ie).latency - events_without_boundaries(ie-3).latency < 2*EEMG.srate && ...
                events_without_boundaries(ie).latency - events_without_boundaries(ie-2).latency < 2*EEMG.srate && ...
                events_without_boundaries(ie).latency - events_without_boundaries(ie-1).latency < 2*EEMG.srate

            event_STsub_latencies_start = [event_STsub_latencies_start; {events_without_boundaries(ie-4).latency}];
            event_STsub_latencies_end = [event_STsub_latencies_end; {events_without_boundaries(ie-2).latency}];

        end
    end

    % add events to EEMG:
    if ~isempty(event_STsub_latencies_start)
        for ie = [event_STsub_latencies_start{:}]
            event_STsub_idx = find([EEMG.event.latency]' == ie);
            EEMG = pop_editeventvals(EEMG, 'changefield', {event_STsub_idx 'type' '9'});
        end

        for ie = [event_STsub_latencies_end{:}]
            event_STsub_idx = find([EEMG.event.latency]' == ie);
            EEMG = pop_editeventvals(EEMG, 'changefield', {event_STsub_idx 'type' '10'});
        end

    end
    % check if last event has a pair
    if rem(str2num(events_without_boundaries(end).type), 2) ~= 0
        event_type = str2num(events_without_boundaries(end).type) + 1;
        if event_type < 11
            event_latency = EEMG.event(end).latency - 10;
            event_code = 'Stimulus';
            n_events = length(EEMG.event);
            EEMG.event(n_events+1).latency = event_latency;
            EEMG.event(n_events+1).code = event_code;
            EEMG.event(n_events+1).type = num2str(event_type);
            EEMG.event(n_events+1).duration = 0;
            EEMG = eeg_checkset(EEMG);
        end
    end

    % do some corrections if needed.
    aa = 1;
    % plot_data(EEMG, 100, 20, 64, 'after event edit', 1)
    % EEMG = pop_editeventvals(EEMG, 'changefield', {2 'type' '9'});
    % EEMG = pop_editeventvals(EEMG, 'delete', [7 8]);
    % n_events = length(EEMG.event);
    % EEMG.event(n_events+1).type = '10';
    % EEMG.event(n_events+1).latency = 16201;
    % EEMG.event(n_events+1).code = 'Stimulus';
    % EEMG.event(n_events+1).duration = 0;
    % EEMG = eeg_checkset(EEMG);

    % remove zeros and 32768
    event_indices = strcmp({EEMG.event.type}, '0') | ...
                    strcmp({EEMG.event.type}, '32768');
    EEMG = pop_editeventvals(EEMG, 'delete', event_indices);
    
    % check if ok
    plot_data(EEMG, 800, EEMG.xmax, 64, 'after event edit', 1)

end

%%

function [event_tbl] = create_event_tbl(EEMG, config)

    event_codes_new = {'ST-Left-Start', 'ST-Left-End', 'DT-Left-Start', 'DT-Left-End', ...
        'ST-Right-Start', 'ST-Right-End', 'DT-Right-Start', 'DT-Right-End', 'C-Start', 'C-End'};
    event_types_new = {'SL1', 'SL2', 'DL1', 'DL2', 'SR1', 'SR2', 'DR1', 'DR2', 'C1', 'C2'};
    
    i = 1;
    for event_type = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10'}
        latencies = {EEMG.event(strcmp({EEMG.event.type}, event_type{:})).latency}';
        types = {EEMG.event(strcmp({EEMG.event.type}, event_type{:})).type}';
        types_new = repmat(event_types_new(i), length(latencies), 1);
        code = repmat(event_codes_new(i), length(latencies), 1);
        subject_entry = repmat({config.subject}, length(latencies), 1);
        task_entry = repmat({config.task}, length(latencies), 1);
    
        if i == 1
            event_tbl = table(subject_entry, task_entry, latencies, types, types_new, code, 'VariableNames', {'subject', 'task', 'latency', 'type', 'type_new', 'code'});
        else
            event_tbl_addition = table(subject_entry, task_entry, latencies, types, types_new, code, 'VariableNames', {'subject', 'task', 'latency', 'type', 'type_new','code'});
            event_tbl = [event_tbl; event_tbl_addition];
        end
        i = i+1;
    end
    event_tbl = sortrows(event_tbl, 'latency');
end

%%

function EEMG = auto_remove_break_segments(EEMG, config)    
    
    if config.do_auto_removal_of_break_segments

        addt = 10; % in points
        event_tbl = create_event_tbl(EEMG, config);

        points_to_keep = [];
        for ievent = 1:size(event_tbl, 1)-1
            current_event_type = event_tbl.type_new{ievent};
            next_event_type = event_tbl.type_new{ievent+1};
            current_lat = event_tbl.latency{ievent};
            next_lat = event_tbl.latency{ievent+1};
    
            if strcmp(current_event_type(1:end-1), next_event_type(1:end-1)) && ...
               str2num(current_event_type(end)) == str2num(next_event_type(end)) - 1     
                points_to_keep = [points_to_keep; current_lat-addt, next_lat+addt];
            end
        end
        % points_to_remove = keep2remove(points_to_keep, EEMG.pnts);
        EEMG = pop_select(EEMG, 'point', points_to_keep);
        % EEMG = eeg_eegrej(EEMG, points_to_remove);
        EEMG = eeg_checkset(EEMG);

    end  
    
    plot_data(EEMG, 800, EEMG.xmax, 64, 'after removal of breaks.', 1)
end


%%

%%

function EEMG = filter_data(EEMG, config)

    if config.remove_line_noise_with_zapline
        try
            [EEMG.data, EEMG.etc.zapline.config, EEMG.etc.zapline.analyticsResults, plothandles] = clean_data_with_zapline_plus(EEMG.data, EEMG.srate, ...
                            struct('noisefreqs', 'line', 'plotResults', 1)); % 'winSizeCompleteSpectrum', floor(EEMG.pnts/8/EEMG.srate)));
            zapline_succesfully_finished = 1;
        catch
            disp('Zapline failed first time, trying later.')
            zapline_succesfully_finished = 0;
        end
    end

    if config.low_freq_cutoff ~= 0
        EEMG = pop_eegfiltnew(EEMG, 'locutoff', config.low_freq_cutoff); 
    end

    if config.high_freq_cutoff ~= 0
        EEMG = pop_eegfiltnew(EEMG, 'hicutoff', config.high_freq_cutoff); 
    end

    if config.remove_line_noise_with_zapline && ~zapline_succesfully_finished
        try
            [EEMG.data, EEMG.etc.zapline.config, EEMG.etc.zapline.analyticsResults, plothandles] = clean_data_with_zapline_plus(EEMG.data, EEMG.srate, ...
                            struct('noisefreqs', 'line', 'plotResults', 1)); % 'winSizeCompleteSpectrum', floor(EEMG.pnts/8/EEMG.srate)));
        catch
            disp('Zapline failed for the second time as well. Skipping.')
        end
    end
    close gcf

    if ~isempty(config.band_freq)
        EEMG = pop_eegfiltnew(EEMG, 'locutoff', config.band_freq(1), ...
                                    'hicutoff', config.band_freq(2), 'revfilt', 1); 
    end
   
    EEMG.etc.filtering_info.lowf = config.high_freq_cutoff;
    EEMG.etc.filtering_info.highf = config.low_freq_cutoff;
    EEMG.etc.filtering_info.bandf = config.band_freq;
    EEMG.etc.filtering_info.zapline_performed = config.remove_line_noise_with_zapline;

end

%% 

function [EEG, EMG] = separate_eemg_data(EEMG)
    EEG = pop_select(EEMG, 'chantype', 'EEG');
    EMG = pop_select(EEMG, 'chantype', 'EMG');
end


%% 

function EEMG = join_eemg_data(EEMG, EEG, EMG, config)
    EEMG.data(1:config.n_eeg_channels, :) = EEG.data;
    EEMG.data(config.n_eeg_channels+1:EEMG.nbchan, :) = EMG.data;
   
    %add etc
    EEMG.etc = EEG.etc;
    EEMG = eeg_checkset(EEMG);

end



%% 

function EEMG = channel_interpolation_wrapper(EEMG, config)
    
    if config.do_channel_interpolation
        keep_cleaning = 1;
    
        while keep_cleaning
            [EEG, EMG] = separate_eemg_data(EEMG);       
            EEG = inner_channel_interpolation(EEG, config, 'EEG');
            
            % put channels to zero
            channels_to_zero = input('Channels to zero: ');
            if ~isempty(channels_to_zero)
                channel_indices_to_zero = find(ismember({EEG.chanlocs.labels}, channels_to_zero));
                EEG.data(channel_indices_to_zero, :) = zeros(length(channel_indices_to_zero), EEG.pnts);
                EEG.etc.chans_to_zero = channels_to_zero;
            end

            EEMG = join_eemg_data(EEMG, EEG, EMG, config);
            EEMG.etc.chans_interp_eeg = EEG.etc.interpolated_channels;
            
            % remove emg channels 
            emg_chan_idx_to_remove = emg_channel_removal(EMG, config);
            if ~isempty(emg_chan_idx_to_remove)
                EEMG = pop_select(EEMG, 'rmchannel', emg_chan_idx_to_remove);
            end
            EEMG.etc.emg_bad_chans_removed = emg_chan_idx_to_remove;

            plot_data(EEMG, 100, 50, 64, 'After cleaning', 1);
            keep_cleaning = input("Would you like to do another round of data inspection and cleaning? [1-yes, 0-no]");
            close all
        end

        EEMG = save_data(EEMG, config, 'after_interp'); 

    else
        disp('Channel interpolation skipped.')
    end

end


%% 

function data = manually_remove_data_segments(data, win_length, spacing, dispchan, config)

   if config.do_manual_removal_of_data_segments
       screen_size = get(0, 'ScreenSize');
       width = screen_size(3);
       height = screen_size(4);
    
       eegplot(data.data, ...
                'eloc_file', data.chanlocs, ...
                'srate', data.srate, ...
                'events', data.event, ...
                'winlength', win_length, ...
                'submean', 'on', ...
                'spacing', spacing, ...
                'dispchans', dispchan, ...
                'command', '1', ...
                'position', [round(0.05*width), round(0.075*height), width-round(0.10*width), height-round(0.15*height)], ...
                'title', ['subject ' config.subject ' | ' config.task '. Here you can remove bad time intervals. Click REJECT at the end.']);

        waitfor( findobj('parent', gcf, 'string', 'REJECT'), 'userdata');
        try
            tmprej = evalin('base', 'TMPREJ');
            data = eeg_eegrej(data, eegplot2event(tmprej));
            data = eeg_checkset(data);
            if ~isfield(data.etc, 'manual_segment_rejection')
                data.etc.manual_segment_rejection = tmprej;
            else
                data.etc.manual_segment_rejection = [data.etc.manual_segment_rejection; tmprej];
            end
        catch
            disp([newline 'Seems like no data was removed' newline])
        end
   end

end

%% 
function emg_chans_to_remove = emg_channel_removal(EMG, config)

    plot_data(EMG, 800, 50, 64, 'After cleaning', 0);
    emg_manual_removal = input('Do you want to do additional EMG channel manual removal  [1-yes, 0-no]? ');
    if emg_manual_removal
        
        emg_chans_to_remove = input(['Give indices of emg channels you want to remove as a mat array (e.g. [1, 2]): ']);        
        % increase to account for eeg channels
        if emg_chans_to_remove == 0
            emg_chans_to_remove = [];
        end
        if ~isempty(emg_chans_to_remove) 
            emg_chans_to_remove = emg_chans_to_remove + config.n_eeg_channels;
        end
    else
        emg_chans_to_remove = [];
    end
    close gcf

end

%%

function EEMG = reref(EEMG, config)
    % separate, reref and join back
    if config.do_reref
        disp('rereferencing EEG data to average.')
        EEG = pop_select(EEMG, 'chantype', 'EEG');
        EEG = pop_reref(EEG, []);
        EEMG.data(1:config.n_eeg_channels, :) = EEG.data; 
    end
end

%%

function EEMG = rectify_emg(EEMG,config)
    if config.rectify_emg
        disp('rectifying emg signal')
        emg_chan_idx = find(strcmp({EEMG.chanlocs.type}, 'EMG'));
        EEMG.data(emg_chan_idx, :) = abs(EEMG.data(emg_chan_idx, :));
    else
        disp('skipping rectification of emg.')
    end
end

%%

function EEMG = do_csd(EEMG, config)
    
    if config.do_csd
        disp('doing csd')
        EEG = pop_select(EEMG, 'chantype', 'EEG');
        EEG = pop_currentsourcedensity(EEG); 
        EEMG.data(1:config.n_eeg_channels, :) = EEG.data;
    else
        disp('skippind csd calculation.')
    end

    EEMG = save_data(EEMG, config, 'after_csd');

end

%%

function [EEMG, event_tbl] = get_event_times(EEMG, config, event_tbl_type)
    
    if strcmp(event_tbl_type, 'during_preproc')
        disp('Getting events...')
        event_codes_new = {'ST-Left-Start', 'ST-Left-End', 'DT-Left-Start', 'DT-Left-End', ...
                           'ST-Right-Start', 'ST-Right-End', 'DT-Right-Start', 'DT-Right-End'};
        event_types_new = {'SL1', 'SL2', 'DL1', 'DL2', 'SR1', 'SR2', 'DR1', 'DR2'};
        i = 1;
        for event_type = {'1', '2', '3', '4', '5', '6', '7', '8'}
            latencies = {EEMG.event(strcmp({EEMG.event.type}, event_type{:})).latency}';
            types = {EEMG.event(strcmp({EEMG.event.type}, event_type{:})).type}';
            types_new = repmat(event_types_new(i), length(latencies), 1);
            code = repmat(event_codes_new(i), length(latencies), 1);
            subject_entry = repmat({config.subject}, length(latencies), 1);
            task_entry = repmat({config.task}, length(latencies), 1);
    
            if i == 1
                event_tbl = table(subject_entry, task_entry, latencies, types, types_new, code, 'VariableNames', {'subject', 'task', 'latency', 'type', 'type_new', 'code'});
            else
                event_tbl_addition = table(subject_entry, task_entry, latencies, types, types_new, code, 'VariableNames', {'subject', 'task', 'latency', 'type', 'type_new','code'});
                event_tbl = [event_tbl; event_tbl_addition];
            end
            i = i+1;
        end
        
        events_without_boundaries = EEMG.event(~strcmp({EEMG.event.type}, 'boundary'));
        events_without_bound_and_zeros = events_without_boundaries(~strcmp({events_without_boundaries.type}, '32768'));
        event_STsub_latencies_start = [];
        event_STsub_latencies_end = [];
        for ie = 4:length(events_without_bound_and_zeros)
            if strcmp(events_without_bound_and_zeros(ie-3).type, '0') && ...
                    strcmp(events_without_bound_and_zeros(ie-2).type, '0') && ...
                    strcmp(events_without_bound_and_zeros(ie-1).type, '0') && ...
                    strcmp(events_without_bound_and_zeros(ie).type, '0') && ...
                    events_without_bound_and_zeros(ie-2).latency - events_without_bound_and_zeros(ie-3).latency < 3*EEMG.srate && ...
                    events_without_bound_and_zeros(ie).latency - events_without_bound_and_zeros(ie-1).latency < 3*EEMG.srate

                event_STsub_latencies_start = [event_STsub_latencies_start; {events_without_bound_and_zeros(ie-2).latency}];
                event_STsub_latencies_end = [event_STsub_latencies_end; {events_without_bound_and_zeros(ie-1).latency}];

            end
        end

        % add events to EEMG:
        if ~isempty(event_STsub_latencies_start)
            for ie = [event_STsub_latencies_start{:}]
                event_STsub_idx = find([EEMG.event.latency]' == ie);
                EEMG = pop_editeventvals(EEMG, 'changefield', {event_STsub_idx 'type' '9'});
                EEMG = pop_editeventvals(EEMG, 'changefield', {event_STsub_idx 'code' 'STsub-None-Start'});
            end
    
            for ie = [event_STsub_latencies_end{:}]
                event_STsub_idx = find([EEMG.event.latency]' == ie);
                EEMG = pop_editeventvals(EEMG, 'changefield', {event_STsub_idx 'type' '10'});
                EEMG = pop_editeventvals(EEMG, 'changefield', {event_STsub_idx 'code' 'STsub-None-End'});
            end
    
            subject_entry = repmat({config.subject}, length(event_STsub_latencies_start), 1);
            task_entry = repmat({config.task}, length(event_STsub_latencies_start), 1);
            type_entry_start = repmat({'9'}, length(event_STsub_latencies_start), 1);
            type_entry_end = repmat({'10'}, length(event_STsub_latencies_end), 1);
            code_entry_start = repmat({'STsub-None-Start'}, length(event_STsub_latencies_start), 1);
            code_entry_end = repmat({'STsub-None-End'}, length(event_STsub_latencies_start), 1);
            type_new_entry_start = repmat({'C1'}, length(event_STsub_latencies_start), 1);
            type_new_entry_end = repmat({'C2'}, length(event_STsub_latencies_end), 1);
        
            event_tbl_addition_start = table(subject_entry, task_entry, event_STsub_latencies_start, ...
                                             type_entry_start, type_new_entry_start, code_entry_start, 'VariableNames', {'subject', 'task', 'latency', 'type', 'type_new', 'code'});
            event_tbl_addition_end = table(subject_entry, task_entry, event_STsub_latencies_end, ...
                                             type_entry_end, type_new_entry_end, code_entry_end, 'VariableNames', {'subject', 'task', 'latency', 'type', 'type_new', 'code'});
            
            event_tbl = [event_tbl; event_tbl_addition_start];
            event_tbl = [event_tbl; event_tbl_addition_end];
        end

        event_tbl = sortrows(event_tbl, 'latency');

        if rem(str2double(event_tbl.type{end}), 2) ~= 0
            event_tbl_addition_row = event_tbl(end, :);
            event_tbl_addition_row.latency = {[EEMG.pnts]};
            event_tbl_addition_row.type = {num2str(str2double(event_tbl_addition_row.type{1}) +1)};
            event_tbl_addition_row.type_new = {[event_tbl_addition_row.type_new{1}(1:2) '2']};
            event_tbl_addition_row.code = strrep(event_tbl_addition_row.code, 'Start', 'End');
            event_tbl = [event_tbl; event_tbl_addition_row];
        end

        % add events to EEMG:
        EEMG.etc.event_tbl = event_tbl;
    else
        event_codes_new = {'ST-Left-Start', 'ST-Left-End', 'DT-Left-Start', 'DT-Left-End', ...
                           'ST-Right-Start', 'ST-Right-End', 'DT-Right-Start', 'DT-Right-End', 'C-Start', 'C-End'};
        event_types_new = {'SL1', 'SL2', 'DL1', 'DL2', 'SR1', 'SR2', 'DR1', 'DR2', 'C1', 'C2'};
        i = 1;
        for event_type = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10'}
            latencies = {EEMG.event(strcmp({EEMG.event.type}, event_type{:})).latency}';
            types = {EEMG.event(strcmp({EEMG.event.type}, event_type{:})).type}';
            types_new = repmat(event_types_new(i), length(latencies), 1);
            code = repmat(event_codes_new(i), length(latencies), 1);
            subject_entry = repmat({config.subject}, length(latencies), 1);
            task_entry = repmat({config.task}, length(latencies), 1);
    
            if i == 1
                event_tbl = table(subject_entry, task_entry, latencies, types, types_new, code, 'VariableNames', {'subject', 'task', 'latency', 'type', 'type_new', 'code'});
            else
                event_tbl_addition = table(subject_entry, task_entry, latencies, types, types_new, code, 'VariableNames', {'subject', 'task', 'latency', 'type', 'type_new','code'});
                event_tbl = [event_tbl; event_tbl_addition];
            end
            i = i+1;
        end
        event_tbl = sortrows(event_tbl, 'latency');
    end
end



%%

function EEG = amica_calculation(EEG, config)
    
    if config.do_amica
        n_interpolated_chans = length(EEG.etc.chans_interp_eeg);
        data_rank = EEG.nbchan - n_interpolated_chans - 1 ;
        
        [w, s, mods] = runamica15(EEG.data,...
            'num_models', 1,...
            'max_threads', 4,...
            'max_iter', config.amica_max_iter,...
            'outdir', config.path_local_amica,...
            'num_chans', EEG.nbchan,...
            'writestep', 2000,...
            'pcakeep', data_rank,...
            'do_reject',0,...
            'write_nd',1,...
            'do_history',0,...
            'histstep',2,...
            'min_dll',1e-9,...
            'min_grad_norm',0.0000005);
    
        EEG.etc.spatial_filter.algorithm = 'AMICA';
        EEG.etc.spatial_filter.AMICAmods = mods;
        EEG.icaweights = w;
        EEG.icasphere = s;
        EEG = eeg_checkset(EEG);
    
    elseif config.prepare_for_hpc_amica

        if rank(EEG.data') < EEG.nbchan
            data_rank = rank(EEG.data');
        else
            interpolated_chans = unique({EEG.chaninfo.removedchans(strcmp({EEG.chaninfo.removedchans.type}, 'EEG')).labels});
            n_interpolated_chans = length(interpolated_chans);
            data_rank = rank(EEG.data') - n_interpolated_chans - 1 ;
        end
 
        set_amica_params([config.path_amica_on_hpc EEG.filename(1:end-4) '.fdt'], ...
                          'outdir', config.path_amicaout_on_hpc, ...
                          'path_param_file', [config.path_amica_for_hpc EEG.filename(1:end-4) '_amica_params'], ...
                          'num_chans', EEG.nbchan, ...
                          'num_frames', EEG.pnts, ...
                          'num_models', 1, ...
                          'max_iter', config.amica_max_iter, ...
                          'num_mix_comps', 3, ...
                          'share_comps', 0, ...
                          'pcakeep', data_rank)

        copyfile([config.path_preproc filesep EEG.filename(1:end-4) '.fdt'], [config.path_amica_for_hpc filesep EEG.filename(1:end-4) '.fdt'])
    end
end

%%

function EEG_to = transfer_components(EEG_from, EEG_to)
    EEG_to.icawinv = EEG_from.icawinv;
    EEG_to.icaweights = EEG_from.icaweights;
    EEG_to.icasphere = EEG_from.icasphere;
    EEG_to = eeg_checkset(EEG_to);
end

%%

%% 

function [EEG, eeg_brain_ics, brain_comps] = decide_on_components(EEG, config)

    if ~config.load_amica_results && ~config.do_amica
        return
    end
    
    if isfield(EEG.etc.spatial_filter, 'rejected_amica_comps')
        % print if there were some components already rejected.
        disp('Currently rejected AMICA comps:')
        disp(EEG.etc.spatial_filter.rejected_amica_comps)

    else
        [EEG, fig_handle] = ic_labeling(EEG, config);
    
        % Reject components (all but brain category above 85% certainty)
        % for label names check "EEG.etc.ic_classification.ICLabel.classes"
        % {'1: Brain', '2: Muscle', '3: Eye', '4: Heart','5: Line Noise', '6: Channel Noise', '7: Other'}
        
        rejComps = [];
        rejCompsLabs = [];
        for comp = 1:length(EEG.etc.ic_classification.ICLabel.classifications)
            [maxVal, maxLab] = max(EEG.etc.ic_classification.ICLabel.classifications(comp, :));
            if  maxLab == 3 && maxVal > config.threshold_for_eye_component_removal
                rejComps(end+1) = comp;
                rejCompsLabs(end+1) = maxLab;
            elseif any(maxLab == [2, 4, 5, 6]) && maxVal > config.threshold_for_other_than_eye_component_removal
                rejComps(end+1) = comp;
                rejCompsLabs(end+1) = maxLab;
            end
        end
        
        % manually check the automatically rejected components
        disp([newline 'These components would be automatically rejected: ' num2str(rejComps) newline])
        doKeep = 1;
        while doKeep && config.do_manual_inspection_of_IC_removal
              doKeep = input('Do you want to KEEP any marked components? yes-1/ no-0.  ');
              if isempty(doKeep); doKeep = 0; end
              if doKeep
                    components_to_keep = input('Enter the components you want to KEEP as a list of numbers seperated by a comma (e.g. [1, 3, 4]).');
                    for ICtoKeep = components_to_keep
                        try
                            rejComps = rejComps(rejComps ~= ICtoKeep);
                        catch
                            disp(['There is no rejected IC labeled' num2str(ICtoKeep) '.'])                        
                        end
                    end
                    disp(['These are the updated components that will be rejected: ' num2str(rejComps)])
              end
        end
    end

    keep_cleaning = 1;
    n_cleans = 0;

    while keep_cleaning
        
        if isfield(EEG.etc.spatial_filter, 'rejected_amica_comps')
            [EEG, fig_handle] = ic_labeling(EEG, config);
            rejComps = [];
            rejCompsLabs = [];
        end

        % plot components time series
        pop_eegplot( EEG, 0, 1, 1);
    
        % additionally manually reject components 
        nextManualRej = 1;
        while config.do_manual_inspection_of_IC_removal && nextManualRej
            nextManualRej = input('Do you want to REJECT more components? yes-1/ no-0.   ');
    
            if isempty(nextManualRej) || nextManualRej == 0
                nextManualRej = 0;
            else
                components_to_remove = input('Enter the components you want to REMOVE as a list of numbers seperated by a comma (e.g. [1, 3, 4]).');
                
                for ICtoRemove = components_to_remove
                    if any(rejComps == ICtoRemove)
                        disp(['The component ' num2str(ICtoRemove) ' is already included for rejection, skipping.'])
                    else
                        try
                            rejComps(end+1) = ICtoRemove;
                            [~, maxLab] = max(EEG.etc.ic_classification.ICLabel.classifications(ICtoRemove, :));
                            rejCompsLabs(end+1) = maxLab;
                            
                        catch
                            disp(['There was an error for the rejection of IC label' num2str(ICtoRemove) ', skipping.'])
                        end
                    end
                end
                disp(['These are the updated components that will be rejected: ' num2str(rejComps)])
    
            end
        end
        
        % add info about which ICs were removed to the EEG struct    
        disp(['Final rejected components: ' num2str(rejComps)])
        if isfield(EEG.etc.spatial_filter, 'rejected_amica_comps')
            n_already_removed_comps = length(EEG.etc.spatial_filter.rejected_amica_comps);
            for irc = 1: length(rejComps)
                EEG.etc.spatial_filter.rejected_amica_comps(n_already_removed_comps+irc, 1) = num2cell(rejComps(irc))';
                EEG.etc.spatial_filter.rejected_amica_comps(n_already_removed_comps+irc, 2) = EEG.etc.ic_classification.ICLabel.classes(rejCompsLabs(irc));
            end  
        else
            EEG.etc.spatial_filter.rejected_amica_comps = num2cell(rejComps)';
            for irc = 1: length(rejCompsLabs)
                EEG.etc.spatial_filter.rejected_amica_comps(irc, 2) = EEG.etc.ic_classification.ICLabel.classes(rejCompsLabs(irc));
            end
        end

        % brain components to save
        if config.save_brain_ICs
            
            brain_comps = [];
            for comp = 1:length(EEG.etc.ic_classification.ICLabel.classifications)
                [maxVal, maxLab] = max(EEG.etc.ic_classification.ICLabel.classifications(comp, :));
                if  maxLab == 1
                    brain_comps(end+1) = comp;
                end
            end
            disp(['Brain components: ' string(brain_comps)])

            additional_brain_comps = input('Add more brain comps, e.g. [1 2]: ');
            brain_comps = [brain_comps, additional_brain_comps];
            EEG.etc.spatial_filter.brain_comps = brain_comps;
            
            % get component activations
            eeg_brain_ics = eeg_getdatact(EEG, 'component', brain_comps); 

            save([config.path_data_preproc filesep config.fname_data_preproc_after_csd '_brain_ICs.mat'], 'eeg_brain_ics')
        end
   
        
        % reject components (they will also be removed from the data)
        if config.reject_components_during_preprocessing
            EEG = pop_subcomp( EEG, rejComps, 0, 0);
            EEG.etc.spatial_filter.rejection_completed = 1;
        else
            EEG.etc.spatial_filter.rejection_completed = 0;
        end
    
        for ifig = 1:size(fig_handle, 2)
            if isvalid(fig_handle(ifig)), close(fig_handle(ifig)); end
        end

        pop_eegplot(EEG, 1, 1, 0);
        keep_cleaning = input("Would you like to do another round of data inspection and cleaning? [1-yes, 0-no]");
        n_cleans = n_cleans + 1;
    end
end

%%

function [EEG, fig_handle] = ic_labeling(EEG, config)
    
    disp('Estimating components'' type using ICLabel package ...');

    EEG = iclabel(EEG, 'lite');

    %  Plot ICA components and optionally save figures
    igroupStart = 1;
    numComps = size(EEG.icawinv,2);
    for igroup = 1 :2 % ceil(numComps/35)
        if igroup == ceil(numComps/35)
            pop_viewprops(EEG, 0, igroupStart:numComps);
        else
            pop_viewprops(EEG, 0, igroupStart:igroup*35);
        end
        igroupStart = igroupStart+35;
        fig_handle(igroup) = gcf;

        if config.save_ic_label_plots
            print( fig_handle(igroup), fullfile( config.path_local_for_hpc_amica, [config.subject '_' config.task '_amica_ICs' num2str(igroup) '.png']), '-dpng')
            % savefig( fig_handle(igroup), fullfile(config.path_local_for_hpc_amica, [config.subject '_' config.task '_amica_ICs' num2str(igroup) '.fig']))
            % close(fig_handle(igroup));
        end
    end

end

%%

function prepare_for_MODA(EEMG, config)
    channels_to_extract = [config.eeg_channels_to_extract, config.emg_channels_to_extract];
    EEMG_final = pop_select(EEMG, 'channel', channels_to_extract);
    writematrix(EEMG_final.data, [EEMG.filepath filesep EEMG.setname '_moda.csv'])
end

%%

function EEMG = save_data(EEMG, config, at_step)
    
    EEMG.etc.config = config;

    if strcmp(at_step, 'after_interp')
        EEMG.setname = [config.fname_data_preproc_after_interp];
        EEMG.filename = [config.fname_data_preproc_after_interp '.set'];
        EEMG.filepath = config.path_data_preproc;
        save([EEMG.filepath filesep EEMG.setname '.mat'], "EEMG")
        % EEG = pop_saveset(EEG, 'filename', EEG.filename, 'filepath', EEG.filepath);
        % EEG = pop_loadset('filename', EEG.filename, 'filepath', EEG.filepath);
    
    elseif strcmp(at_step, 'after_mclean')
        EEMG.setname = [config.fname_data_preproc_after_mclean];
        EEMG.filename = [config.fname_data_preproc_after_mclean '.set'];
        EEMG.filepath = config.path_data_preproc;
        save([EEMG.filepath filesep EEMG.setname '.mat'], "EEMG")
        % EEMG = pop_saveset(EEMG, 'filename', EEMG.filename, 'filepath', EEMG.filepath);
        % EEG = pop_loadset('filename', EEG.filename, 'filepath', EEG.filepath);

    elseif strcmp(at_step, 'after_csd')
        EEMG.setname = [config.fname_data_preproc_after_csd];
        EEMG.filename = [config.fname_data_preproc_after_csd '.set'];
        EEMG.filepath = config.path_data_preproc;
        save([EEMG.filepath filesep EEMG.setname '.mat'], "EEMG")
        % EEMG = pop_saveset(EEMG, 'filename', EEMG.filename, 'filepath', EEMG.filepath);
        % EEG = pop_loadset('filename', EEG.filename, 'filepath', EEG.filepath);
        disp([newline newline 'Preprocessing until amica has successfully ended for: ' newline ...
            'Subject: ' EEMG.subject newline 'Condition: '  EEMG.session ' & ' EEMG.condition newline newline]); 

    elseif strcmp(at_step, 'after_amica_clean')
        EEMG.setname = [config.fname_data_preproc_after_amica_clean];
        EEMG.filename = [config.fname_data_preproc_after_amica_clean '.set'];
        EEMG.filepath = config.path_data_preproc;
        save([EEMG.filepath filesep EEMG.setname '.mat'], "EEMG")
        % EEMG = pop_saveset(EEMG, 'filename', EEMG.filename, 'filepath', EEMG.filepath);
        % EEG = pop_loadset('filename', EEG.filename, 'filepath', EEG.filepath);
    
    elseif strcmp(at_step, 'after_breakcut')
        EEMG.setname = [config.fname_data_preproc_after_breakcut];
        EEMG.filename = [config.fname_data_preproc_after_breakcut '.set'];
        EEMG.filepath = config.path_data_preproc;
        save([EEMG.filepath filesep EEMG.filename], "EEMG")
        % EEMG = pop_saveset(EEMG, 'filename', EEMG.filename, 'filepath', EEMG.filepath);
    
    elseif strcmp(at_step, 'after_taskcut')
        EEMG.setname = [config.fname_data_preproc_after_taskcut];
        EEMG.filename = [config.fname_data_preproc_after_taskcut '.set'];
        EEMG.filepath = config.path_data_preproc;
        save([EEMG.filepath filesep EEMG.filename], "EEMG")
        % EEMG = pop_saveset(EEMG, 'filename', EEMG.filename, 'filepath', EEMG.filepath);
    elseif strcmp(at_step, 'brain_ics')
        EEMG.setname = [config.fname_data_preproc_after_brain_ics];
        EEMG.filename = [config.fname_data_preproc_after_brain_ics '.set'];
        EEMG.filepath = config.path_data_preproc;
        save([EEMG.filepath filesep EEMG.setname '.mat'], "EEMG")
    end
end

%%

function save_EEMG_with_brain_comps(EEMG, eeg_brain_ics, brain_comps, config)
    n_brain_comps = size(eeg_brain_ics, 1);
    EEMG = pop_select(EEMG, 'rmchannel', [1:(config.n_eeg_channels - n_brain_comps)]);
    EEMG.data(1:n_brain_comps, :) = single(eeg_brain_ics);
    for ic = 1:n_brain_comps
        EEMG.chanlocs(ic).labels = ['IC' num2str(brain_comps(ic))];
    end
    save_data(EEMG, config, 'brain_ics');
end

%%

function [EEMG, keep_data] = decide_to_keep_data_overall(EEMG)
    keep_data = input('Do you want to keep data as a whole [1-yes, 0-no]? ');
    if ~keep_data
        EEMG = pop_select(EEMG, 'point', [0 1]);
        EEMG.comments = "This data was to bad. It was removed.";
    end
end


%% 

function prepare_for_amica(EEMG, config)
    
    if config.prepare_for_hpc_amica

        % do additional removal of segments and interpolation just for amica. 
        EEMG = channel_interpolation_wrapper(EEMG, config);
        EEMG = manually_remove_data_segments(EEMG, 10, 150, config); 
        save([config.path_data_preproc filesep EEMG.setname '_cleaned_for_amica.mat'], "EEMG");
    
        % save data
        [EEG, ~] = separate_eemg_data(EEMG);
          
        pop_saveset(EEG, 'filename', EEG.filename, 'filepath', config.path_local_for_hpc_amica);
        delete([config.path_local_for_hpc_amica filesep EEG.filename])
    
        % rank
        n_interpolated_chans = length(EEG.etc.chans_interp_eeg);
        data_rank = EEG.nbchan - n_interpolated_chans - 1 ;
        
        % save AMICA settings in the file input.param
        set_amica_params(['amica_on_hpc/' EEG.filename(1:end-4) '.fdt'], ...
                          'outdir', ['amica_on_hpc/' EEG.filename(1:end-4) '_amicaout'], ...
                          'path_param_file', [config.path_local_for_hpc_amica filesep EEG.filename(1:end-4) '_amica_params'], ...
                          'num_chans', EEG.nbchan, ...
                          'num_frames', EEG.pnts, ...
                          'num_models', 1, ...
                          'max_iter', config.amica_max_iter, ...
                          'num_mix_comps', 3, ...
                          'share_comps', 0, ...
                          'pcakeep', data_rank)
    else
        disp('Skipping amica preparation')
    end
end

%%

function [EEMG, EEG] = load_data_and_amica_results(config)
    
    % load data
    load([config.path_data_preproc filesep config.fname_data_preproc_after_csd]);
    EEMG = eeg_checkset(EEMG);
    EEG = pop_select(EEMG, 'chantype', 'EEG');

    % load amica results
    modout = loadmodout15(config.path_local_for_hpc_amica_outdir);
    EEG.etc.spatial_filter.algorithm = 'AMICA';
    EEG.etc.spatial_filter.AMICAmods = modout;

    EEG.etc.amica = modout;
    EEG.icaweights = modout.W;
    EEG.icasphere = modout.S(1:modout.num_pcs, :);
    EEG.icawinv = modout.A;
    EEG.icachansind = 1:EEG.nbchan;
    EEG = eeg_checkset(EEG);
end


%%

function split_data_in_tasks(EEMG_orig, config)
    
    if config.additional_resampling_before_task_split ~= 0
        setname = EEMG_orig.setname;
        previous_fs_txt = ['fs' num2str(EEMG_orig.srate)];
        new_fs_txt = ['fs' num2str(config.additional_resampling_before_task_split)];
        EEMG_orig = pop_resample(EEMG_orig, config.additional_resampling_before_task_split);
        if contains(setname, '_fs')
            fname_preproc_data_new = strrep(setname, previous_fs_txt, new_fs_txt);
        else
            fname_preproc_data_split = strsplit(setname, '_');
            fname_first_part = strjoin(fname_preproc_data_split(1:4), '_');
            fname_second_part = strjoin(fname_preproc_data_split(5:end), '_');
            fname_preproc_data_new = [fname_first_part '_' new_fs_txt '_' fname_second_part];
        end
        EEMG_orig.setname = fname_preproc_data_new;
    end

    if config.do_data_task_splits
        event_codes_new = {'ST-Left-Start', 'ST-Left-End', 'DT-Left-Start', 'DT-Left-End', ...
            'ST-Right-Start', 'ST-Right-End', 'DT-Right-Start', 'DT-Right-End', 'C-Start', 'C-End'};
        event_types_new = {'SL1', 'SL2', 'DL1', 'DL2', 'SR1', 'SR2', 'DR1', 'DR2', 'C1', 'C2'};
        task_names = {'SL', 'DL', 'SR', 'DR', 'C'};
        start_task_types = {'1', '3', '5', '7', '9'};
    
        for i = 1:length(start_task_types)
            start_task_type = start_task_types{i};
            start_events = EEMG_orig.event(strcmp({EEMG_orig.event.type}, start_task_type));
            end_events = EEMG_orig.event(strcmp({EEMG_orig.event.type}, num2str(str2double(start_task_type)+1)));
            task_periods = [];
            for period = 1:length(start_events)
                task_periods = [task_periods; start_events(period).latency-2 end_events(period).latency+2];
            end
            if ~any(task_periods)
                continue
            end
            if size(task_periods, 1) ~= 1
                task_periods = task_periods';
            end
            if config.make_equal_segments
                if ~strcmp(start_task_type, '9'); segment_duration = 33; else;  segment_duration = 27; end %seconds
                segment_n_points = segment_duration*EEMG_orig.srate;
                duration_of_task_periods = diff(task_periods);
                if duration_of_task_periods < segment_n_points
                    warning('Duration of task segment is smaller than requested! ')
                    too_short = duration_of_task_periods < segment_n_points*0.9;
                    if any(too_short)
                        for i = find(too_short)
                            disp(['Duration of task ' num2str(i) ' too short. Removing.'])
                            task_periods(:, i) = [];
                            duration_of_task_periods(i) = [];
                        end
                    end
                end
                if isempty(task_periods); continue; end
                middle_of_task_periods = task_periods(1, :) + floor(duration_of_task_periods/2);
                equal_start = floor(middle_of_task_periods-(segment_n_points/2));
                equal_end = ceil(middle_of_task_periods+(segment_n_points/2));
                task_periods = [];
                for period = 1:length(equal_start)
                    task_periods = [task_periods; equal_start(period) equal_end(period)];
                end
                if length(equal_start) == 2; figure, hold on, plot(EEMG_orig.data(140, :)), xline(equal_start(1)), xline(equal_start(2)), xline(equal_end(1)), xline(equal_end(2)); end
            end
            EEMG_orig = eeg_checkset(EEMG_orig);
            EEMG = pop_select(EEMG_orig, 'point', task_periods);
            EEMG.session = task_names{i};
            
            if config.make_equal_segments
                EEMG.setname = [EEMG_orig.setname '_eqsplit' task_names{i}];
            else
                EEMG.setname = [EEMG_orig.setname '_split' task_names{i}];
            end
            
            save([config.path_data_preproc filesep EEMG.setname '_nomc.mat'], 'EEMG');
            EEMG = manually_remove_data_segments(EEMG, 80, 800, 192, config); % win_length, spacing, dispchan
       %     EEMG = manual_channel_interpolation(EEMG, config, 'EEG');
            EEMG = manually_remove_data_segments(EEMG, 20, 100, 64, config); 
            
            if config.make_equal_segments
                EEMG.setname = [EEMG_orig.setname '_eqsplit' task_names{i}];
            else
                EEMG.setname = [EEMG_orig.setname '_split' task_names{i}];
            end

            save([config.path_data_preproc filesep EEMG.setname '_mc.mat'], 'EEMG');
            disp(['Saved: ' task_names{i}])
            clear EEMG;
        end
    end
end


%% 

function plot_data(EEMG, spacing, winlength, dispchans, fig_title, wait_to_close)
   
   screen_size = get(0, 'ScreenSize');
   width = screen_size(3);
   height = screen_size(4);

   eegplot(EEMG.data, ...
            'srate', EEMG.srate, ...
            'events', EEMG.event, ...
            'winlength', winlength, ...
            'submean', 'on', ...
            'spacing', spacing, ...
            'dispchans', dispchans, ...
            'position', [round(0.05*width), round(0.075*height), width-round(0.10*width), height-round(0.15*height)], ...
            'title', fig_title);

   if wait_to_close; uiwait(gcf); end
end


%%


function data = manual_channel_interpolation(data, config, data_type)
    
    n_manual_intepolation_loops = 0;
    plot_data_eeglab(data, 20, 64, 100, ['subject ' config.subject ' | ' config.task '. Plotting interpolated data'], 0)

    repeat_manual_interpolation = input('Do you want to do additional manual interpolation [1-yes, 0-no]? ');
    if iscell(repeat_manual_interpolation)
        repeat_manual_interpolation = input('Do you want to do additional manual interpolation [1-yes, 0-no]? ');
    end

    while repeat_manual_interpolation
      
        chans_idx_to_interp_manually = [];
        chans_to_interp_manually = input(['Give names of channels you want to manually interpolate as a cell array ' ...
                                          '(e.g. {''T8'', ''F1''}, not indices): ']);
        if ischar(chans_to_interp_manually)
            chans_to_interp_manually = cellstr(chans_to_interp_manually);
        end
        
        if isstring(chans_to_interp_manually{1})
          for ich = 1:length(chans_to_interp_manually)
              chans_to_interp_manually{ich} = char(chans_to_interp_manually{ich});
          end
        end

        if ~isempty(chans_to_interp_manually)
            for channel = chans_to_interp_manually
                try
                    % chans_idx_to_interp_manually(end+1) = data.chanlocs(find(ismember({data.chanlocs.labels}, channel))).urchan;
                    chans_idx_to_interp_manually(end+1) = find(ismember({data.chanlocs.labels}, channel));
                    disp(['Electrode with the name ' char(channel) ' will be interpolated.'])
                catch
                    disp(['There is no electrode with the name: ' char(channel) '. It was not included. Try again.'])
                end
            end

            % Interpolate missing channels
            if ~isempty(chans_idx_to_interp_manually)
                if n_manual_intepolation_loops > 1
                    data = pop_interp(data, chans_idx_to_interp_manually, 'invdist');
                end
                data = pop_interp(data, chans_idx_to_interp_manually);
                data = eeg_checkset(data);
                % save info on which channels were interpolated
                data = add_interpolation_info(data, chans_idx_to_interp_manually);
            end
    
            close gcf
            disp('Plotting interpolated data...')
            if strcmp(data_type, 'EEG'); spacing = 100; else; spacing = 500; end
            plot_data_eeglab(data, 50, 64, spacing, ['subject ' config.subject ' | ' config.task '. Plotting interpolated data'], 0)
     
        end
        n_manual_intepolation_loops = n_manual_intepolation_loops + 1;
        repeat_manual_interpolation = input('Do you want to do additional manual interpolation [1-yes, 0-no]? ');
    end

    if n_manual_intepolation_loops == 0; data = add_interpolation_info(data, []); end
    close gcf
end

%%

function data = add_interpolation_info(data, chans_idx_to_interp)
    channel_labels = {data.chanlocs(chans_idx_to_interp).labels};
    if ~isfield(data.etc, 'interpolated_channels')
        if isempty(channel_labels)
            data.etc.interpolated_channels = [];
        else
            data.etc.interpolated_channels = channel_labels;
        end
    else
        data.etc.interpolated_channels = unique([data.etc.interpolated_channels, channel_labels]);
    end
end

%%


function plot_data_eeglab(EEG, win_length, disp_chans, spacing, fig_title, wait_to_close)
    
    screen_size = get(0, 'ScreenSize');
    width = screen_size(3);
    height = screen_size(4);

    eegplot(EEG.data, ...
        'eloc_file', EEG.chanlocs, ...
        'srate', EEG.srate, ...
        'events', EEG.event, ...
        'winlength', win_length, ...
        'dispchans', disp_chans, ...
        'submean', 'on', ...
        'spacing', spacing, ...
        'position', [round(0.05*width), round(0.075*height), width-round(0.10*width), height-round(0.15*height)], ...
        'title', fig_title);
    
    if wait_to_close, uiwait(gcf); end

end
