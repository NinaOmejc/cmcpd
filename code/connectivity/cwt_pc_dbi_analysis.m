function cwt_pc_dbi_analysis(config, idx_job)
    
    % set paths and load data
    [config] = set_up_job(config, idx_job);
    load([config.path_preproc filesep config.fname_data_preproc '.mat'], 'EEMG')

    if config.do_cwt
        do_cwt(EEMG, config);
    end
    
    if config.do_pc
        do_pc(EEMG, config);
    end

    if config.do_dbi
        [chans_eeg, chans_emg] = get_chans(EEMG, config, 'dbi');
        [data, datanames] = get_data_of_channels(EEMG, chans_eeg, chans_emg, config);
        data_filt = filter_data_in_bands(data, EEMG.srate, config);
        phases = extract_phases(data_filt, config);
        do_dbi(phases, datanames, EEMG.srate, config);
    end

end


%%


function [config] = set_up_job(config, idx_job)
    
    % set the combination settings in the config structure
    config.idx_job = idx_job;
    config.comb = config.combinations(idx_job, :);

    for ic = 1:size(config.combinations, 2)
        
        if isnan(str2double(config.combinations{idx_job, ic}))
            config.(config.combination_names{ic}) = char(config.combinations{idx_job, ic});
        else
            config.(config.combination_names{ic}) = str2double(config.combinations{idx_job, ic});
        end
    end
    
    if ~isfield(config, 'pc_n'), config.pc_n = 0; end; if ~isfield(config, 'pc_m'), config.pc_m = 0; end
    if ~isfield(config, 'dbi_win_size'), config.dbi_win_size = 0; end

    % set paths based on the combination
    config = set_paths(config);
     
    if config.verbose; disp(newline); fprintf('\nWorking on Job idx: %d | Task: %s | Subject: %s \n\t\t\t\t\t  | Recify: %d | CSD: %d | LF: %d | HF: %d | ics: %d \n\t\t\t\t\t  | split: %s | pc_n: %d | pc_m: %d | dbi_win_size: %d \n', idx_job, config.task, ...
            config.subject, config.with_rectify_emg, config.with_csd, config.with_low_freq_cutoff, config.with_high_freq_cutoff, config.with_brain_ics, config.task_split, config.pc_n, config.pc_m, config.dbi_win_size); end
 end


%%


function config = set_paths(config)
   
    if config.with_rectify_emg; rect_text = 'rect'; else rect_text = 'norect'; end
    if config.with_auto_removal_of_break_segments; data_type = 'onlytask'; else; data_type = 'allduration'; end
    if config.with_resample_freq ~= 0; fstxt = ['fs' num2str(config.with_resample_freq)]; else; fstxt = 'fs2000'; end
    if config.with_detrend; detrend_text = 'dtrnd'; else; detrend_text = 'nodtrnd'; end
    if config.with_low_freq_cutoff ~= 0; lowfiltstext = ['lf' num2str(config.with_low_freq_cutoff)]; else lowfiltstext = 'lf0'; end
    if config.with_high_freq_cutoff ~= 0; highfiltstext = ['hf' num2str(config.with_high_freq_cutoff)]; else highfiltstext = 'hf0'; end
    if config.with_csd ~= 0; csdtext = 'csd'; else; csdtext = 'nocsd';end
    if config.with_chan_interp ~= 0; interptxt = 'interp'; else interptxt = 'nointerp'; end
    if config.with_manual_removal_of_data_segments ~= 0; mctxt = 'mc'; else mctxt = 'nomc'; end
    if config.with_rerefs ~=0; rereftxt = 'reref'; else; rereftxt = 'noreref'; end
    if config.with_amica_clean ~= 0; amicacleantxt = '_icacl'; else; amicacleantxt = ''; end
    if config.with_brain_ics ~= 0; brainicstxt = '_icatdic'; else brainicstxt = ''; end
    if config.with_equal_segments; eqtxt = 'eq'; else; eqtxt = ''; end
    if strcmp(config.task_split, 'none'), split_txt = ''; else; split_txt = [eqtxt 'split' config.task_split]; end
    if config.dbi_win_size; dbi_winsize_txt = ['dbiwin' num2str(config.dbi_win_size)]; else; dbi_winsize_txt = ''; end
    preproc_name_base = [data_type '_' fstxt '_' detrend_text '_' lowfiltstext '_' highfiltstext '_' interptxt '_' split_txt '_' mctxt amicacleantxt brainicstxt];

    config.path_data_base = [config.path_root filesep 'data' filesep];
    config.path_preproc = [config.path_data_base filesep 'preproc' filesep char(config.task) filesep char(config.subject)];
    config.path_cwt = [config.path_data_base filesep 'cwt' filesep char(config.task) '_' config.version '_' split_txt filesep char(config.subject)];
    config.path_pc = [config.path_data_base filesep 'pc' filesep char(config.task) '_' config.version '_' split_txt filesep char(config.subject)];
    config.path_dbi = [config.path_data_base filesep 'dbi' filesep char(config.task) '_' config.version '_' split_txt filesep char(config.subject)];      
    config.path_surr_pc = [config.path_data_base filesep 'pc_surr'];
    config.fullfname_surr_pc = [config.path_data_base filesep 'pc_surr' filesep 'pc_surr_size' num2str(config.pc_npnts_surr) '_fs' num2str(config.with_resample_freq) '_n' num2str(config.pc_n) '_m' num2str(config.pc_m) '_all.mat'];
    
    config.path_pc_single = [config.path_pc filesep 'pc_single_emg'];
    config.fullfname_freqs = [config.path_data_base filesep 'cwt' filesep 'time_freqs.mat'];
    config.fname_data_preproc = [config.subject '_' config.task '_eemg_' preproc_name_base];
    config.fname_data_cwt = [config.subject '_' config.task '_' config.version '_cwt'];
    config.fname_data_pc = [config.subject '_' config.task '_' config.version '_pc'];
    config.fname_data_dbi = [config.subject '_' config.task '_' config.version '_dbi'];

    if ~exist(config.path_cwt, 'dir') && config.do_cwt;      mkdir(config.path_cwt); end
    if ~exist(config.path_pc, 'dir') && config.do_pc;        mkdir(config.path_pc); end
    if ~exist(config.path_pc_single, 'dir') && config.pc_save_single_emg_results ; mkdir(config.path_pc_single); end
    if ~exist(config.path_dbi, 'dir') && config.do_dbi;      mkdir(config.path_dbi); end
    if ~exist(config.path_surr_pc, 'dir');                   mkdir(config.path_surr_pc); end
   
    if contains(config.path_data_base, '\') && ~exist([config.path_data_base filesep 'preproc_to_analysis_version_control_' config.version '.txt'], 'file')
        fileID = fopen([config.path_data_base filesep 'preproc_to_analysis_version_control_' config.version '.txt'], 'w');
        fprintf(fileID, ['Analysis version ' config.version ' corresponds to preprocessing version: \n']);
        fprintf(fileID, [ preproc_name_base '\n']);
        fclose(fileID);
    end
end


%%


function do_cwt(EEMG, config)

    if config.verbose; disp('Started cwt analysis.'); end
    if isfield(EEMG.etc, 'event_tbl'); event_tbl = EEMG.etc.event_tbl; else; event_tbl = []; end
    
    % get appropriate channels 
    [chans_eeg, chans_emg] = get_chans(EEMG, config, 'cwt');
    
    % get appropriate data
    if isstruct(chans_eeg) || isstruct(chans_emg)
        [data, datanames] = get_data_of_channels(EEMG, chans_eeg, chans_emg, config);
        EEMG.data = [data.eeg; data.emg];
        cwt_chans = [datanames.eeg, datanames.emg];
        
        for icluster = 1: size(EEMG.data, 1)
            EEMG.chanlocs(icluster).labels = cwt_chans{icluster};
        end
    else
        cwt_chans = [chans_eeg; chans_emg];
    end
    
    % do cwt for each channel
    for idx = 1:length(cwt_chans)
        
        % find channel idx
        channel_label = cwt_chans{idx};
        channel_idx = find(strcmp({EEMG.chanlocs.labels}, channel_label), 1, 'first');
        if isempty(channel_idx); continue; end

        % check if already exist
        if exist([config.path_cwt filesep config.fname_data_cwt '_' channel_label '.mat'], 'file') && ~config.force_recompute_cwt
            disp(['WT of ' channel_label ' already exists. Skipping...'])
            continue            
        end
        
        % do cwt
        if config.save_figures; make_plot = 'pow+tm'; else; make_plot = false; end
        fig_title = ['CWT | ' config.task ' | ' config.task_split ' | ' config.subject ' | ' channel_label ' | lf: ' num2str(config.with_low_freq_cutoff) ' | hf: ' num2str(config.with_high_freq_cutoff) ' | brain ics: ' num2str(config.with_brain_ics), ' | '];
        disp(fig_title)

        [WT, f] = wt(EEMG.data(channel_idx, :), EEMG.srate, ...
                        'fmin', config.cwt_freq_min, 'fmax', config.cwt_freq_max, ...
                        'Preprocess', 'off', 'Wavelet', 'Morlet', 'f0', config.cwt_f0, ...
                        'Plot', make_plot, 'CutEdges', 'on', 'nv', config.cwt_nv, ...
                        'event_tbl', event_tbl, 'Padding', 'none', 'fig_title', fig_title, 'Display', 'notify');
                    
        % save results and plots
        if config.save_to_dbx
            upload_to_dbx([config.path_cwt filesep config.fname_data_cwt '_' channel_label '.mat'], WT, config)
        else
            if config.cwt_save_power_insteadof_complex
                if config.cwt_reduce_size_before_save 
                    WT = abs(WT).^2;
                    WT = single(downsample(WT', config.cwt_ds_factor)');
                end
            end
            save([config.path_cwt filesep config.fname_data_cwt '_' channel_label '.mat'], "WT", "f")
            if config.save_figures; saveas(gcf,[config.path_cwt filesep config.fname_data_cwt '_' channel_label '.png']); end
        end
    end
end



%% 


function do_pc(EEMG, config)
    
    if config.verbose; disp('Started pc analysis.'); end
    if isfield(EEMG.etc, 'event_tbl'); event_tbl = EEMG.etc.event_tbl; else; event_tbl = []; end
    
    % get appropriate channels
    [chans_eeg, chans_emg] = get_chans(EEMG, config, 'pc');
    
    % get appropriate data
    if isstruct(chans_eeg) || any(contains(chans_eeg, 'IC'))
        [data, datanames] = get_data_of_channels(EEMG, chans_eeg, chans_emg, config);
        EEMG.data(1:size(data.eeg, 1), :) = data.eeg;
        chans_eeg = datanames.eeg;
        for icluster = 1: size(data.eeg, 1)
            EEMG.chanlocs(icluster).labels = chans_eeg{icluster};
        end
    end
    
    % for through each EEG channel
    for eeg_channel_label = chans_eeg
        
        % get eeg channel index
        eeg_channel_label = eeg_channel_label{:};
        eeg_channel_idx = find(strcmp({EEMG.chanlocs.labels}, eeg_channel_label));
        if isempty(eeg_channel_idx); continue; end

        % check if EEG CWT already exist
        if config.pc_allow_to_load_cwt && exist([config.path_cwt filesep config.fname_data_cwt '_' eeg_channel_label '.mat'], 'file')
            disp(['CWT of ' eeg_channel_label ' already exists. Loading...'])
            load([config.path_cwt filesep config.fname_data_cwt '_' eeg_channel_label '.mat'])
        
        % else first calculate EEG CWT
        else
            fig_title = ['CWT | ' config.task ' | ' config.task_split ' | ' config.subject ' | ' eeg_channel_label ' | lf: ' num2str(config.with_low_freq_cutoff) ' | hf: ' num2str(config.with_high_freq_cutoff) ' | brain ics: ' num2str(config.with_brain_ics), ' | '];
            [WT, f] = wt(EEMG.data(eeg_channel_idx, :), EEMG.srate, ...
                            'fmin', config.cwt_freq_min, 'fmax', config.cwt_freq_max, 'Preprocess', 'on', ...
                            'Wavelet', 'Morlet', 'f0', config.cwt_f0, 'Plot','off', 'CutEdges', 'on', 'nv', config.cwt_nv, ...
                            'event_tbl', event_tbl, 'Padding', 'none', 'fig_title', fig_title, 'Display', 'notify');
            
            % save([config.path_cwt filesep config.fname_data_cwt '_' eeg_channel_label '.mat'], "WT", "f")
            % saveas(gcf,[config.path_cwt filesep config.fname_data_cwt '_' eeg_channel_label '.png']); close gcf
        end
        WT_eeg = WT;  clear WT;
        
        % else, go through all emg channels and create variables to store
        % all emg pc results (to do at the end average over all EMG)
        pc_allemg = {{}, {}};                           % 1 - right leg, 2 - left leg
        tpc_avgemg_right = [];
        tpc_avgemg_left = [];
        
        if isstruct(chans_emg); chans_emg = [chans_emg.emgr, chans_emg.emgl]; end
        emg_exists = 0;
        for emg_channel_label = chans_emg
            emg_channel_label = emg_channel_label{:};
            if ~isempty(regexp(emg_channel_label, 'II-', 'once')) || strcmp(emg_channel_label, 'emgl'); muscle_location = 'left'; else; muscle_location = 'right'; end
            
            % check if EEG to AVG EMG TPC already exist and skip if so
            if exist([config.path_pc filesep config.fname_data_pc '_' eeg_channel_label '_allemg_' muscle_location '.mat'], 'file') && ~config.force_recompute_pc
                disp(['AVG TPC ' config.task_split ' | ' config.subject ' | between ' eeg_channel_label ' and all emgs already exists. Skipping.'])
                emg_exists = 1;
                break
            end

            % check if already exist EMG CWT analysis
            if config.pc_allow_to_load_cwt && exist([config.path_cwt filesep config.fname_data_cwt '_' emg_channel_label '.mat'], 'file')
                disp(['CWT of ' emg_channel_label ' already exists. Loading...'])
                load([config.path_cwt filesep config.fname_data_cwt '_' emg_channel_label '.mat'])
            
            % else calculate EMG CWT 
            else
                emg_channel_idx = find(strcmp({EEMG.chanlocs.labels}, emg_channel_label));
                if isempty(emg_channel_idx); continue; end

                fig_title = ['CWT | ' config.task ' | ' config.task_split ' | ' config.subject ' | ' emg_channel_label ' ' muscle_location ' | lf: ' num2str(config.with_low_freq_cutoff) ' | hf: ' num2str(config.with_high_freq_cutoff) ' | brain ics: ' num2str(config.with_brain_ics), ' | '];
                disp(fig_title)

                [WT, f] = wt(EEMG.data(emg_channel_idx, :), EEMG.srate, ...
                    'fmin', config.cwt_freq_min, 'fmax', config.cwt_freq_max, 'Preprocess', 'on', ...
                    'Wavelet', 'Morlet', 'f0', config.cwt_f0, 'Plot','off', 'CutEdges', 'on', 'nv', config.cwt_nv, ...
                    'event_tbl', event_tbl, 'Padding', 'none', 'fig_title', fig_title, 'Display', 'notify');
                
                % save
                % save([config.path_cwt filesep config.fname_data_cwt '_' emg_channel_label '.mat'], "WT")
                % saveas(gcf,[config.path_cwt filesep config.fname_data_cwt '_' emg_channel_label '.png']);  close gcf
            end

            WT_emg = WT;  clear WT;

            % -----------------
            % phase coherence
            % -----------------

            % check if already exist EEG to single EMG TPC
            if config.pc_save_single_emg_results && config.save_figures; pc_plot = 1; else; pc_plot = 0; end
            
            % either load or calculate surrogates, if needed
            if config.pc_add_surrogates
                n_pnts_data = size(WT_emg, 2);
                if config.pc_create_surrogates_on_spot
                    phcoh_surrogates = create_surrogates_phcoh(config, EEMG.srate, n_pnts_data, length(f));
                    points_diff = 0;
                else
                    load(config.fullfname_surr_pc)
                    points_diff = n_pnts_data - config.pc_npnts_surr;
                end
                if points_diff >= 0
                    points_diff_txt = ['Data longer than surr for ' num2str(points_diff) ' pnts.']; 
                else 
                    points_diff_txt = ['Data shorter than surr for ' num2str(abs(points_diff)) ' pnts.'];
                end                
            else
                phcoh_surrogates = [];
                points_diff = 0;
                points_diff_txt = '';
            end
            
            fig_title = ['PC ' config.task ' | ' config.subject ' | ' config.task_split ' | rect: ' num2str(config.with_rectify_emg) 'csd: ' num2str(config.with_csd) ...
                        '| lf: ' num2str(config.with_low_freq_cutoff) ' | hf: ' num2str(config.with_high_freq_cutoff) ' | ' eeg_channel_label ' | ' ...
                        emg_channel_label ' ' muscle_location ' | n:' num2str(config.pc_n) ' | m:' num2str(config.pc_m) newline points_diff_txt];
            disp(fig_title)
            
            % do TPC for single EMG channel
            [TPC, PC] = tlphcoh(WT_emg, WT_eeg, f, ...
                                        config.pc_n, config.pc_m, ...
                                        EEMG.srate, pc_plot, event_tbl, ...
                                        phcoh_surrogates, fig_title); 
            
            if config.pc_save_single_emg_results
                save([config.path_pc_single filesep config.fname_data_pc '_' eeg_channel_label '_' emg_channel_label '.mat'], "TPC", "PC", "f")
                saveas(gcf,[config.path_pc_single filesep config.fname_data_pc '_' eeg_channel_label '_' emg_channel_label '.png'])
                %saveas(gcf,[config.path_pc filesep config.fname_data_pc '_' eeg_channel_label '_' emg_channel_label '.fig'])
                close gcf;
            end

            % add to average 
            if strcmp(muscle_location, 'right')
                pc_allemg{1} = [pc_allemg{1}; PC];
                if isempty(tpc_avgemg_right)
                    tpc_avgemg_right = TPC; 
                else
                    tpc_avgemg_right = tpc_avgemg_right + TPC; 
                end                        
            else
                pc_allemg{2} = [pc_allemg{2}; PC];
                if isempty(tpc_avgemg_left)
                    tpc_avgemg_left = TPC; 
                else
                    tpc_avgemg_left = tpc_avgemg_left + TPC; 
                end    
            end

            clear TPC WT_emg WT PC
        end
        if emg_exists; continue; end
            
        % finish calculating average over TPC and phoch (divide by N)
        tpc_avgemg_right = tpc_avgemg_right ./ length(pc_allemg{1});
        tpc_avgemg_left = tpc_avgemg_left ./ length(pc_allemg{2});
        pc_avgemg_right = nanmean(cell2mat(pc_allemg{1}), 1);
        pc_avgemg_left = nanmean(cell2mat(pc_allemg{2}), 1);

        % plot average
        for imuscle = 1:2
            if imuscle == 1
                muscle_location = 'right';
                tpc_avgemg = tpc_avgemg_right;
                pc_avgemg = pc_avgemg_right;
                pcs = pc_allemg{1};
            else
                muscle_location = 'left'; 
                tpc_avgemg = tpc_avgemg_left;
                pc_avgemg = pc_avgemg_left;
                pcs = pc_allemg{2};
            end
            
            % plot and save figures 
            if config.save_figures
                fig_title_avg = ['AVG TPC ' config.task ' | ' config.subject ' | ' config.task_split ' | rect: ' num2str(config.with_rectify_emg) ' | csd: ' num2str(config.with_csd) ...
                            '| lf: ' num2str(config.with_low_freq_cutoff) ' | hf: ' num2str(config.with_high_freq_cutoff) ' | n:' num2str(config.pc_n) ' | m:' num2str(config.pc_m) ' | ' newline ...
                            ' | ' eeg_channel_label ' | all EMGs on ' muscle_location ' | ' points_diff_txt];
                
                plot_average_phcoh(tpc_avgemg, pc_avgemg, pcs, f, ...
                    EEMG.srate, event_tbl, phcoh_surrogates, fig_title_avg, config)
                saveas(gcf, [config.path_pc filesep config.fname_data_pc '_' eeg_channel_label '_allemg_' muscle_location '.png'])
            end
            
            % save results
            if config.save_to_dbx
                upload_to_dbx([config.path_pc filesep config.fname_data_pc '_' eeg_channel_label '_allemg_' muscle_location '_tpc_avgemg.mat'], tpc_avgemg, config)
                upload_to_dbx([config.path_pc filesep config.fname_data_pc '_' eeg_channel_label '_allemg_' muscle_location '_pc_avgemg.mat'], pc_avgemg, config)
            else
                % downsample to 50 Hz and transform from double to single (to reduce size)
                if config.pc_reduce_size_before_save 
                    tpc_avgemg = single(downsample(tpc_avgemg', config.pc_ds_factor)');
                end
                save([config.path_pc filesep config.fname_data_pc '_' eeg_channel_label '_allemg_' muscle_location '.mat'], "tpc_avgemg", "pc_avgemg", "f")
            end
            clear pc_avgemg tpc_avgemg
        end
        clear WT_eeg
    end
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
    elseif strcmp(analysis_type, 'ac')
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
    elseif any(contains(chans_eeg, 'IC'))
        if ~isfield(config, 'file_to_ic_clusters') || isempty(config.file_to_ic_clusters)
            data.eeg = eeg_getdatact(EEMG, 'component', 1:size(chans_eeg, 2));
            datanames.eeg = chans_eeg;
        else
            load(config.file_to_ic_clusters)
            comp_indices = [];
            for i = 1:length(chans_eeg)
                ic = chans_eeg{i};
                idx = find(contains(cluster_info.(ic).cluster_setname, EEMG.setname));
                if ~isempty(idx)
                    comp_indices(i) = cluster_info.(ic).cluster_compidx(idx);
                else
                    comp_indices(i) = nan;                
                end
            end
            data.eeg = eeg_getdatact(EEMG, 'component', comp_indices);
            datanames.eeg = chans_eeg;
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

%%

function phcoh_surrogates = create_surrogates_phcoh(config, fs, data_len, n_freqs)
    
    % check if already exists, if so, load them
    fullfname_pc_surrs = [config.path_surr_pc filesep  'pc_surr_' config.subject '_' config.task_split '.mat']; 
    if exist(fullfname_pc_surrs, 'file')
        disp('Surr data already exist. Loading.')
        load(fullfname_pc_surrs)
        return
    end
    
    % if not, calculate and save
    rng(1); 
    i = 1;
    surr = randn(config.pc_n_surrogates, data_len);
    n_pc_surrs = (config.pc_n_surrogates * (config.pc_n_surrogates-1))/2;
    phcoh_surrogates = zeros(n_pc_surrs, n_freqs);
    
    for idx_S1 = 1:config.pc_n_surrogates

        [WT_S1, f] = wt(surr(idx_S1, :), fs, ...
            'fmin', config.cwt_freq_min, 'fmax', config.cwt_freq_max, 'Preprocess', 'off', ...
            'Wavelet', 'Morlet', 'f0', config.cwt_f0, 'Plot', 'off', 'CutEdges', 'on', 'nv', config.cwt_nv, ...
            'Padding', 'none', 'fig_title', ['surr ' num2str(idx_S1)], 'event_tbl', [], 'Display', 'notify');

        % CWT2 for S2
        for idx_S2 = idx_S1+1:config.pc_n_surrogates

            [WT_S2, f] = wt(surr(idx_S2, :), fs, ...
                'fmin', config.cwt_freq_min, 'fmax', config.cwt_freq_max, 'Preprocess', 'on', ...
                'Wavelet', 'Morlet', 'f0', config.cwt_f0, 'Plot', 'off', 'CutEdges', 'on', 'nv', config.cwt_nv, ...
                'Padding', 'none', 'fig_title', ['surr ' num2str(idx_S2)], 'event_tbl', [], 'Display', 'notify');

            % phase coherence
            disp(['Doing surrogate analysis for S1: ' num2str(idx_S1), ' | S2: ' num2str(idx_S2)])
            [~, phcoh_surrogate] = tlphcoh(WT_S1, WT_S2, f, config.pc_n, config.pc_m, fs, 0, [], '');
            phcoh_surrogates(i, :) = phcoh_surrogate;
            i = i+1;
            clear WT_S2
        end
        clear WT_S1
    end
    % save TPC
    save(fullfname_pc_surrs, "phcoh_surrogates")
end



%%

function plot_average_phcoh(tpc_avgemg, pc_avgemg, pcs, f, ...
            fs, event_tbl, phcoh_surrogates, fig_title, config)
    
    L = size(tpc_avgemg, 2);
    YY=f; XX=(0:(L-1))/fs; ZZ=abs(tpc_avgemg);  ZZ=ZZ.^2;

    scrsz=get(0,'ScreenSize'); 
    figure('Position',[scrsz(3)/4,scrsz(4)/8,scrsz(3)/1.5,6*scrsz(4)/8], 'Visible', config.make_plot_visible);

    axes('Position',[0.1,0.15,0.6,0.75],'Layer','top','YScale','log','Box','on','FontSize',18);
    TL=length(XX); FL=length(YY);
    pc=pcolor(XX,YY,ZZ); set(pc,'EdgeColor','none'); %title(ZZname);
    set(gca,'yscale','log')
    ylabel('Frequency (Hz)'); xlabel('Time (s)'); sgtitle(fig_title);
    xlim([0,(L-1)/fs]);
    ylims = [4 150];
    hold on;
    ax = gca;
    % clim([0 ax.CLim(2)])
    %mycolormap = customcolormap(linspace(0,1,6), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd'});
    % mycolormap = customcolormap(linspace(0,1,6), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3'});
    mycolormap = jet;
    cb = colorbar('horiz');
    colormap(mycolormap);
    cb.Position = [cb.Position(1), cb.Position(2)-0.1, cb.Position(3), cb.Position(4)];

    axis tight
    ax = axis;
    sirnina_crt=0.5;
    frekvence = [7 12 30 50];
    y_tick_labels = {'4', '7', '12', '30', '50', '100', '150'};
    y_ticks = [4 7 12 30 50 100 150];

    for ifreq = 1:size(frekvence, 2)
        l = yline(frekvence(ifreq));
        set(l,'color','black','linewidth',sirnina_crt,'linestyle','--');
    end
    yticklabels(y_tick_labels)
    yticks(y_ticks)
    ylim(ylims)
    
    % added by nina
    if ~isempty(event_tbl)
        for i = 1:size(event_tbl, 1)
            latency = cell2mat(event_tbl.latency(i)) / fs;
            l = line([latency latency], gca().YLim);
            set(l,'color','k','linewidth',sirnina_crt,'linestyle','--');
            if i == 1
                mp = latency;
            elseif rem(i, 2) == 0
                mp = latency-1000/fs;
            else
                mp = latency+1000/fs;
            end
            t = text(mp, max(ylim), event_tbl.type_new{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
            set(t,'color', 'k')
        end
    end

    axes('Position',[0.75,0.15,0.2,0.75],'Layer','top','YScale','log','Box','on','FontSize',18);
    % mline=plot(phcoh',YY','-k','LineWidth',2); xlabel({'Time-averaged WPC'});
    set(gca,'yscale','log','YTickLabel', y_tick_labels, 'YTick', y_ticks)
    ylim(ylims);
    max_coh = max(pc_avgemg(f > 4));
    if max_coh > 0.2; xlim([0, max_coh]); else; xlim([0 0.2]); end

    %axis tight
    ax = axis;
    for ifreq = 1:size(frekvence, 2)
        l = yline(frekvence(ifreq));
        set(l,'color','black','linewidth',sirnina_crt,'linestyle','--');
    end
    yticklabels(y_tick_labels)
    yticks(y_ticks)
    xlabel('Phase Coherence')

    if ~isempty(phcoh_surrogates)

        true_signal = pc_avgemg';
        phcoh_surr_mean = mean(phcoh_surrogates, 1);
        phcoh_surr_2std = phcoh_surr_mean + 2*std(phcoh_surrogates, 1);
        n_phcohs = size(phcoh_surrogates, 1);
        
        hold on;
        for i_emg = 1:length(pcs)
            plot(pcs{i_emg, :}, f, 'color', '#A9A9A9', 'LineStyle', '-')
        end
        plot(phcoh_surr_mean, YY', 'color', 'b', 'LineStyle', '-', 'LineWidth', 2);
        plot(phcoh_surr_2std, YY', 'color', 'b', 'LineStyle', '--', 'LineWidth', 2)
        plot(true_signal, f, 'r-', 'LineWidth',2)
        hold on;
        shade_yx(gca, YY', true_signal, YY', phcoh_surr_2std,'FillType',[1 2], 'FillColor', {'red'}, 'FillAlpha', 0.5)
    end

end

%%

function sigs = filter_data_in_bands(data, fs, config) 
    
    sigs = struct();
    datatype_names = fieldnames(data);
    fband_names = fieldnames(config.dbi_freq_bands);
    % take different orders of the filters for different freq bands
    fband_ors = struct('alpha', 3, 'beta', 5, 'gamma', 6, 'mu', 3, 'beta_low', 4, 'beta_high', 5, 'gamma_low', 6, 'gamma_high', 6);
    
    for idatatype = 1:numel(datatype_names)
        datatype_name = datatype_names{idatatype};
        for ifrband = 1:numel(fband_names)
            fband_name = fband_names{ifrband};
            fband_or = fband_ors.(fband_name);
            fband_lf = config.dbi_freq_bands.(fband_name)(1);
            fband_hf = config.dbi_freq_bands.(fband_name)(2);
            
            % go through every time series and filter it 
            sigs.([datatype_name '_' fband_name]) = zeros(size(data.(datatype_name)));
            for its = 1:size(data.(datatype_name), 1)
                ts = data.(datatype_name)(its, :);
                ts_filt = dbi_bandpass(ts, fband_or, fband_lf, fband_hf, fs, 0);

                % if it fails, adjust order of the filter (lower it) and do it again
                if ((sum(isnan(ts_filt))>0)||(max(abs(ts_filt))>10000))
                    fband_or = fband_or - 1;
                    ts_filt = dbi_bandpass(ts, fband_or, fband_lf, fband_hf, fs, 0);
                    disp(['Reducing the order for ' fband_name ' to ' num2str(fband_or)])
                % else
                %     fband_or = fband_or + 1;
                %     ts_filt = dbi_bandpass(ts, fband_or, fband_lf, fband_hf, fs, 0);
                %     disp(['Increasing the order for ' fband_name ' to ' fband_or])
                end
                % disp(['Final filter order used for ' fband_name ' is: ' num2str(fband_or)]);

                % save filtered data
                if (sum(isnan(ts_filt))==0) && (max(abs(ts_filt)<10000))
                    sigs.([datatype_name '_' fband_name])(its, :) = ts_filt;
                end
            end
        end
    end
end

%%

function phases = extract_phases(data_filt, config)
    datanames = fieldnames(data_filt);
    for ifield = 1:numel(datanames)
        dataname = datanames{ifield};
        % phases.(dataname) = zeros(size(data_filt.(dataname)));
        for its = 1:size(data_filt.(dataname), 1)
            iphase = co_hilbproto(data_filt.(dataname)(its, :), 0, 0, 0, 0);
            [iphi, ~, ~] = co_fbtransf1(iphase', 256);
            phases.(dataname)(its, :) = iphi; % unwrap(iphi);
            % plot phases
            % figure,plot(iphase)
        end
    end
end


%%

function do_dbi(phases, channel_names, fs, config)
    
    % set names
    phases_fieldnames = fieldnames(phases);
    eeg_namefields = phases_fieldnames(contains(phases_fieldnames, 'eeg'));
    emg_namefields = phases_fieldnames(contains(phases_fieldnames, 'emg'));
    eeg_channames = channel_names.eeg; % config.dbi_chans_eeg;
    emg_channames = channel_names.emg; % {'emgl', 'emgr'}
      
    % -----------
    % DO DBI FOR EMG - EEG 
    % -----------

    for ieeg = 1:length(eeg_namefields)
        ieeg_name = eeg_namefields{ieeg};
        for iemg = 1:length(emg_namefields)
            iemg_name = emg_namefields{iemg};
            % if contains(iemg_name, 'alpha');  continue;    end      % skipping EMG alpha combination
            % if contains(ieeg_name, 'gamma_high');  continue;    end      % skipping EMG alpha combination
            phi1_eeg = phases.(eeg_namefields{ieeg});
            phi2_emg = phases.(emg_namefields{iemg});
            for ieeg_chans = 1:size(phi1_eeg, 1)
                for iemg_chans = 1:size(phi2_emg, 1)
                    
                    % set namings and paths
                    ieeg_chan_name = eeg_channames{ieeg_chans};
                    iemg_chan_name = emg_channames{iemg_chans};

                    dbi_path_base = [config.path_dbi filesep iemg_chan_name];
                    dbi_fname_base = [ieeg_name '_' ieeg_chan_name '_' iemg_name '_' iemg_chan_name];
                    dbi_fig_title = [config.task_split ' | ' config.subject ' | ' strrep([ieeg_name ': ' ieeg_chan_name], '_', '-') ' | ' strrep([iemg_name ': ' iemg_chan_name], '_', '-')];
                    if ~exist(dbi_path_base, 'dir'); mkdir(dbi_path_base); end
                    
                    % check if already exist and skip in that case
                    if exist([dbi_path_base filesep config.fname_data_dbi '_' dbi_fname_base '.mat'], 'file') && exist([dbi_path_base filesep config.fname_data_dbi '_' dbi_fname_base '_cstrength_withsurr.mat'], "file")
                        disp(['Skipping DBI for: ' dbi_fig_title '.'])
                        continue
                    end

                    % -------------------
                    % DBI on orig data
                    % -------------------

                    disp(['Doing DBI for ' dbi_fig_title '.'])
                    phi1 = phi1_eeg(ieeg_chans, :);
                    phi2 = phi2_emg(iemg_chans, :);
                    
                    [tm, cc, e] = bayes_main(unwrap(phi1), unwrap(phi2), ...
                                             config.dbi_win_size, ...
                                             1/fs, config.dbi_winoverlap, 0.2, 0, 2);
                    
                    % save results
                    if config.save_to_dbx
                        upload_to_dbx([dbi_path_base filesep config.fname_data_dbi '_' dbi_fname_base '.mat'], cc, config)
                    else
                        save([dbi_path_base filesep config.fname_data_dbi '_' dbi_fname_base '.mat'], "tm", "cc", "e")
                    end 

                    % calculate the coupling strength on orig data
                    nwins = size(cc, 1);
                    for iwin=1:nwins 
                        [cpl_strength1(iwin), cpl_strength2(iwin), drc] = dirc(cc(iwin,:), 2); 
                    end

                    % check if store coupling functions for every time point                   
                    if config.save_figures
                        if config.dbi_save_func_over_time_plots
                            if ~exist([dbi_path_base filesep 'cfplots_' dbi_fname_base ], 'dir'); mkdir([dbi_path_base filesep 'cfplots_' dbi_fname_base], 'dir'); end
                        end
                        
                        % calulate coupling func (and maybe save) for every time point
                        for iwin = 1:nwins
                            [q1, q2]= CFprint(cc(iwin, :), 2, config.dbi_save_func_over_time_plots, config.make_plot_visible);
                            if ~exist('q1_mean', 'var'); q1_mean=q1; else; q1_mean = q1_mean + q1; end
                            if ~exist('q2_mean', 'var'); q2_mean=q2; else; q2_mean = q2_mean + q2; end
                            
                            if config.dbi_save_func_over_time_plots
                                sgtitle([dbi_fig_title ' | win ' num2str(iwin)])
                                saveas(gcf, [dbi_path_base filesep 'cfplots_' dbi_fname_base filesep dbi_fname_base '_win' num2str(iwin) '.png'])
                                close gcf
                            end
                        end
    
                        % calculate and potentially plot and save coupling func that is averaged over time 
                        path_to_save = [dbi_path_base filesep config.fname_data_dbi '_' dbi_fname_base '_cfplotavg_q2.png'];
                        q2_mean = q2_mean ./ nwins;
                        plot_qq(q2_mean, dbi_fig_title, path_to_save, 'q_2(\phi_1,\phi_2)', config)
    
                        path_to_save = [dbi_path_base filesep config.fname_data_dbi '_' dbi_fname_base '_cfplotavg_q1.png'];
                        q1_mean = q1_mean ./ nwins;
                        plot_qq(q1_mean, dbi_fig_title, path_to_save, 'q_1(\phi_1,\phi_2)', config)
                    end
                    
                    % ----------------------
                    % do DBI on surrogates 
                    % ----------------------

                    if config.dbi_compute_n_surrogates ~= 0
                        
                        % calculate surragates
                        surrs_phi1 = surrogate(phi1, config.dbi_compute_n_surrogates, 'MCPP');
                        surrs_phi2 = surrogate(phi2, config.dbi_compute_n_surrogates, 'MCPP');
                        
                        % compute DBI for each surrogate 
                        surrs_cpl_strength1 = zeros(config.dbi_compute_n_surrogates, nwins);
                        surrs_cpl_strength2 = zeros(config.dbi_compute_n_surrogates, nwins);

                        for idx_surr = 1:config.dbi_compute_n_surrogates
                            
                            % disp(['Doing DBI Surr ' num2str(idx_surr) '/' num2str(config.dbi_compute_n_surrogates) ' between ' dbi_fig_title '.'])

                            [~, surr_cc, ~] = bayes_main(surrs_phi1(idx_surr, :), ...
                                                         surrs_phi2(idx_surr, :), ...
                                                         config.dbi_win_size, 1/fs, config.dbi_winoverlap, 0.2, 0, 2);

                            % calculate the coupling strength on surr data
                            for iwin=1:nwins 
                                [surr_cpl_strength1(iwin), surr_cpl_strength2(iwin), drc] = dirc(surr_cc(iwin,:), 2); 
                            end
                            surrs_cpl_strength2(idx_surr, :) = surr_cpl_strength2;
                            surrs_cpl_strength1(idx_surr, :) = surr_cpl_strength1;
                        end
                        
                        % --------
                        % compare surrogates to orig data
                        % --------
                        if config.save_figures
                            % calculate for direction 1 -> 2
                            path_to_save = [dbi_path_base filesep config.fname_data_dbi '_' dbi_fname_base '_cstrength2_withsurr.png'];
                            plot_coupling_strength(tm, cpl_strength2, surrs_cpl_strength2, config, dbi_fig_title, path_to_save)
                            % calculate for direction 2 -> 1
                            path_to_save = [dbi_path_base filesep config.fname_data_dbi '_' dbi_fname_base '_cstrength1_withsurr.png'];
                            plot_coupling_strength(tm, cpl_strength1, surrs_cpl_strength1, config, dbi_fig_title, path_to_save)
                        end

                        % save coupling strength and surrogates
                        if config.save_to_dbx
                            cpl_strength2_withsurrs = [tm; cpl_strength2; mean(surrs_cpl_strength2, 1); std(surrs_cpl_strength2, 1)];
                            upload_to_dbx([dbi_path_base filesep config.fname_data_dbi '_' dbi_fname_base '_cstrength2_withsurr.mat'], cpl_strength2_withsurrs, config)
                            cpl_strength1_withsurrs = [cpl_strength1, mean(surrs_cpl_strength1, 1), std(surrs_cpl_strength1, 1)];
                            upload_to_dbx([dbi_path_base filesep config.fname_data_dbi '_' dbi_fname_base '_cstrength1_withsurr.mat'], cpl_strength1_withsurrs, config)
                        else
                            save([dbi_path_base filesep config.fname_data_dbi '_' dbi_fname_base '_cstrength_withsurr.mat'], "tm", "cpl_strength2", "surrs_cpl_strength2", "cpl_strength1", "surrs_cpl_strength1")
                        end

                    else
                        % if no surrogates, just save results of coupling strength for original data
                        if config.save_figures
                            % calculate for direction 1->2
                            path_to_save = [dbi_path_base filesep config.fname_data_dbi '_' dbi_fname_base '_cstrength2_nosurr.png'];
                            plot_coupling_strength(tm, cpl_strength2, [], config, dbi_fig_title, path_to_save)
                            % calculate for direction 2->1
                            path_to_save = [dbi_path_base filesep config.fname_data_dbi '_' dbi_fname_base '_cstrength1_nosurr.png'];
                            plot_coupling_strength(tm, cpl_strength1, [], config, dbi_fig_title, path_to_save)
                        end
                        if config.save_to_dbx
                            upload_to_dbx([dbi_path_base filesep config.fname_data_dbi '_' dbi_fname_base '_cstrength2_nosurr.mat'], cpl_strength2, config)
                            upload_to_dbx([dbi_path_base filesep config.fname_data_dbi '_' dbi_fname_base '_cstrength1_nosurr.mat'], cpl_strength1, config)
                        else
                            save([dbi_path_base filesep config.fname_data_dbi '_' dbi_fname_base '_cstrength_nosurr.mat'], "tm", "cpl_strength2", "cpl_strength1")
                        end
                    end
                end
            end
        end
    end     
end


%% 

function plot_qq(qq, dbi_fig_title, path_to_save, zlabeltxt, config)    
    if ~config.save_to_dbx
        t1=0:0.13:2*pi; t2=0:0.13:2*pi;  
        figure('Visible', config.make_plot_visible), 
        surf(t1,t2, qq,'FaceColor','interp');                                               
        view([-40 50])
        set(gca,'fontname','Helvetica','fontsize',12,'Xgrid','off','Ygrid','off')            
        xlabel('\phi_1');ylabel('\phi_2');zlabel(zlabeltxt);axis tight
        colormap(hot)
        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
        sgtitle(['AVG CouplFunc | ' dbi_fig_title])
        saveas(gcf, path_to_save)    
        close gcf
    end
end

%%

function plot_coupling_strength(tm, cpl_strength, surrs_cpl_strength2, config, dbi_fig_title, path_to_save)
    
    if ~config.save_to_dbx
        figure('Visible', config.make_plot_visible); hold on;
        plot(tm, cpl_strength, 'k-', 'LineWidth',2);
        
        if ~isempty(cpl_strength)
            surrs_cpl_strength2_mean = mean(surrs_cpl_strength2);
            surrs_cpl_strength2_2std = surrs_cpl_strength2_mean + 2*std(surrs_cpl_strength2);
            plot(tm, surrs_cpl_strength2_mean, 'b-')
            plot(tm, surrs_cpl_strength2_2std, 'b--')
        end
    
        xlabel('Time'); ylabel('Coupling strength')
        title(['CPL strength | ' dbi_fig_title])
        saveas(gcf, path_to_save)
        close gcf;
    end
end

%% 
function upload_to_dbx(filename, result, config)   
    if config.save_to_dbx
        try
            filename_csv = strrep(filename, '.mat', '.csv');
            writematrix(result, filename_csv)
            resp = uploadToDropbox(config.dbx_access_token, filename_csv);
            if isstruct(resp); delete(filename_csv); end 
        catch e
            fprintf(1,'The identifier was:\n%s',e.identifier);
            fprintf(1,'There was an error! The message was:\n%s',e.message);
            disp('Upload to dbx didnt succeed.')
        end
    end
end


%% 
function labels = index2labels(EEMG, index)
    labels = {EEMG.chanlocs(1:128).labels};
end