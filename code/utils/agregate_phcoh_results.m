% agregate over EMG channels

path_data_base = 'D:\Experiments\corticomuscular_analysis\data\real';

% load frequencies
fullfname_freqs = [path_data_base filesep 'timefreq' filesep 'time_freqs.mat'];
load(fullfname_freqs)

% find all conditions
path_phasecoh = [path_data_base filesep 'phasecoh'];
files_phasecoh = dir([path_phasecoh filesep 'IZO*']);
conds = {files_phasecoh(:).name};

% go through each condition
for icond = 1:length(conds)
    cond = conds{icond};
    if contains(cond, 'splitC') || ~contains(cond, 'split')
        continue
    end
    subjects = dir([path_phasecoh filesep cond filesep 'PD*']);

    % check what kind of surrogate data (short or long) and load appropriate
    if contains(cond, '_split')
        path_surr_phasecoh = [path_data_base filesep 'phasecoh_surr' filesep 'phasecoh_surr_split_n1_m1'];
        fullfname_surr_data = [path_surr_phasecoh filesep 'phasecoh_surr_split_n1_m1_all.mat'];
    else
        path_surr_phasecoh = [path_data_base filesep 'phasecoh_surr' filesep 'phasecoh_surr_n1_m1'];
        fullfname_surr_data = [path_surr_phasecoh filesep 'phasecoh_surr_n1_m1_all.mat'];
    end
    load(fullfname_surr_data)

    % go through each subjects
    for isubject = 1:length(subjects)
        subject = subjects(isubject).name;
        subject_files = dir([path_phasecoh filesep cond filesep subject filesep '*.mat']);

        path_phasecoh_agregated = [path_phasecoh filesep cond filesep subject filesep 'aggregated'];
        if ~exist(path_phasecoh_agregated, 'dir'); mkdir(path_phasecoh_agregated); end

        % get all eeg channel names
        chans_eeg = {};
        for ifile = 1:length(subject_files)
            ifile_split = strsplit(subject_files(ifile).name, '_');
            chans_eeg = [chans_eeg; ifile_split{end-1}];
        end
        chans_eeg = unique(chans_eeg);

        % go through each eeg channel, load phcoh results and group them by emg
        % channels
        for ichan_eeg = 1:length(chans_eeg)
            chan_eeg = chans_eeg{ichan_eeg};
            chan_eeg_files = subject_files(contains({subject_files(:).name}, ['_' chan_eeg '_']));

            % 1 - right leg, 2 - left leg
            phcoh_emg = {{}, {}};
            for ichan_eeg_file = 1:length(chan_eeg_files)
                load([chan_eeg_files(ichan_eeg_file).folder filesep chan_eeg_files(ichan_eeg_file).name])
                if contains(chan_eeg_files(ichan_eeg_file).name, 'II-')
                    phcoh_emg{2} = [phcoh_emg{2}; phcoh];
                else
                    phcoh_emg{1} = [phcoh_emg{1}; phcoh];
                end
                clear phcoh;
            end

            % plot coherence for all emg channels, also average and surrogates.
            phcoh_surr_mean = mean(phcoh_surrogates, 1);
            phcoh_surr_2std = phcoh_surr_mean + 2*std(phcoh_surrogates, 1);
            phcoh_emg_mean = {mean(cell2mat(phcoh_emg{1}), 1), mean(cell2mat(phcoh_emg{2}), 1)};

            for imuscle = 1:2
                if imuscle == 1; muscle_location = 'right'; else; muscle_location = 'left'; end
                iphcoh_emg = phcoh_emg{imuscle};
                iphcoh_emg_mean = phcoh_emg_mean{imuscle};

                figure,
                hold on;
                plot(f, phcoh_surr_mean, 'color', 'b', 'LineStyle', '-', 'LineWidth', 2)
                plot(f, phcoh_surr_2std, 'color', 'b', 'LineStyle', '--', 'LineWidth', 2)

                for i_emg = 1:length(iphcoh_emg)
                    plot(f, iphcoh_emg{i_emg, :}, 'color', '#A9A9A9', 'LineStyle', '-')
                end
                plot(f, iphcoh_emg_mean)
                shade_xy(f, iphcoh_emg_mean, f, phcoh_surr_2std, 'FillType',[1 2], 'FillColor', {'red'}, 'FillAlpha', 0.5)
                set(gca,'xscale','log','XTickLabel', {'1', '7', '15', '30', '50', '100', '150'}, 'XTick', [1 7 15 30 50 100 150])
                xlim([7, 150])
                ylim([0 0.15])
                for xl = [7 15 30 50 150]
                    l = xline(xl, 'k--');
                end
                fig_title = ['PhCoh IZO | ' subject ' | rect: 0 | csd: 0 | lf: 1 | hf: 150 | ' newline ...
                    chan_eeg ' | all emg channels ' muscle_location ];
                title(fig_title)
                saveas(gcf, [path_phasecoh_agregated filesep 'phcoh_aggregated_all_emg_' muscle_location '_' chan_eeg '.png'])
                close gcf

            end
        end
    end
end