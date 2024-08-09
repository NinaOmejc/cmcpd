
function data = channel_interpolation(data, config, data_type)
    
    if ~config.do_channel_interpolation,   return; end

    % ---------------------------- 
    % EEG data

    if strcmp(data_type, 'EEG')
        
        % automatic check
        chans_idx_to_interp_auto = auto_channel_interpolation_check(data, config);
        
        % dont auto exclude central channels 
        chans_to_interp_labels = {data.chanlocs(chans_idx_to_interp_auto).labels};
        central_chans_to_be_kept_mask = contains(chans_to_interp_labels, 'C');
        chans_idx_to_interp_auto = chans_idx_to_interp_auto(~central_chans_to_be_kept_mask);

        chans_to_interp_labels_new = {data.chanlocs(chans_idx_to_interp_auto).labels};
        fprintf('\n All automatically detected bad channels after centrality criteria:\n')
        disp(string(chans_to_interp_labels_new))

        % manual check of automatic check and manual interpolation
        if config.do_manual_inspection_of_channel_interpolation
            
            % manual check of automatic check
            chans_idx_to_interp_auto = check_of_auto_detected_bad_channels(data, chans_idx_to_interp_auto);
    
            if ~isempty(chans_idx_to_interp_auto)
                data = pop_interp(data, chans_idx_to_interp_auto);
                % data = pop_interp(data, chans_idx_to_interp_auto, 'invdist');
                data = eeg_checkset(data);
                data = add_interpolation_info(data, chans_idx_to_interp_auto);

                disp([newline 'Plotting auto interpolated data...' newline])
                plot_data_eeglab(data, 50, 64, 100, ['subject ' config.subject ' | ' config.task '. Check of automatically interpolated data'], 0)
            else
                disp("No channels were interpolated based on automatic detection method.")
            end
            
            % manual interpolation option
            data = manual_channel_interpolation(data, config, data_type);

        % if no manual check of automatic check, just interpolate
        else
            data = pop_interp(data, chans_idx_to_interp_auto, 'spherical');
            % data = pop_interp(data, chans_idx_to_interp_auto, 'invdist');
            data = eeg_checkset(data);
            data = add_interpolation_info(data, chans_idx_to_interp_auto);
        end
    
    % ---------------------------- 
    % EMG data (only manual check)
    elseif strcmp(data_type, 'EMG')

        if config.do_manual_inspection_of_channel_interpolation
            plot_data_eeglab(data, 50, 64, 500, ['subject ' config.subject ' | ' config.task '. Plotting interpolated data'], 0)
            data = manual_channel_interpolation(data, config, data_type);
        else
            disp('EMG channel interpolation skipped, cause it can only be manual.')
        end

    end

end

%%

function chans_idx_to_interp_auto = auto_channel_interpolation_check(data, config)

    % Find bad channels automatically
    % rereference and check bad channels on temporary file for better channel noise detection
    data_temp = data;
    data_temp.nbchan = data_temp.nbchan+1;  % + one for virtual
    data_temp.data(end+1,:) = zeros(1, data_temp.pnts);
    data_temp.chanlocs(1,data_temp.nbchan).labels = 'initialReference';
    data_temp = pop_reref(data_temp, []);
    data_temp = pop_select(data_temp,'nochannel',{'initialReference'});
    
    [~, ~, ~, chans_to_interp_auto] = clean_artifacts(data_temp,...
        'BurstCriterion', 'off', ...
        'WindowCriterion', 'off',...
        'FlatlineCriterion',config.flatLineCriterion,...
        'ChannelCriterion', config.channelCriterion,...
        'LineNoiseCriterion', 'off', ...
        'highpass_band',[0.25 0.75],...
        'ChannelCriterionMaxBadTime', config.channelCriterionMaxBadTime);
    
    chans_idx_to_interp_auto = find(chans_to_interp_auto); % transform logical array to indices
    chans_to_interp_labels = {data.chanlocs(chans_to_interp_auto).labels};

    fprintf('\n All automatically detected bad channels:\n')
    disp(string(chans_to_interp_labels))

    plot_data_eeglab(data, 50, 64, 100, ['subject ' config.subject ' | ' config.task '. Data before automatic interpolation'], 0);
    % win_length = 50
    % disp_chans = 64
    % fig_title = ['subject ' config.subject ' | ' config.task '. Data before automatic interpolation'];
    % wait_to_close = 1;
end


%%

function chans_to_interpolate_idx = check_of_auto_detected_bad_channels(data_temp, chans_to_interp_auto)
    chans_to_interpolate_labels = {data_temp.chanlocs(chans_to_interp_auto).labels};
    do_keep = 1;
    while do_keep
          do_keep = input('Do you want to KEEP any currently bad channels? yes-1/ no-0.  ');
          if do_keep
              if isnumeric(do_keep) 
                  channels_to_keep = input([newline 'Enter the channel names you want to KEEP as a cell of strings ' ...
                                          newline 'seperated by a comma (e.g. {''C0P2'', ''FFC5h''}):']);
                  if ischar(channels_to_keep)
                      channels_to_keep = cellstr(channels_to_keep);
                  end
                  
                  if isstring(channels_to_keep{1})
                      for ich = 1:length(channels_to_keep)
                          channels_to_keep{ich} = char(channels_to_keep{ich});
                      end
                  end
            
                  for chan_to_keep = channels_to_keep
                      chans_to_interpolate_labels = chans_to_interpolate_labels(~strcmp(chans_to_interpolate_labels, chan_to_keep{1}));
                  end
              end
              disp([newline 'These are the updated channels that will be interpolated: ' newline])
              disp(chans_to_interpolate_labels)
          end
    end
    chans_to_interpolate_idx = find(ismember({data_temp.chanlocs.labels}, chans_to_interpolate_labels));
    close gcf
end

%%

function data = manual_channel_interpolation(data, config, data_type)
    
    n_manual_intepolation_loops = 0;
    
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
