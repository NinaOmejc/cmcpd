function main_connectivity(fullfname_config, idx_job)

    addpath(genpath('D:\Experiments\corticomuscular_analysis\cmcpd\code'))
    
    % load configuration file 
    load(fullfname_config)
    config = get_combinations(config);
    
    % either run locally with local matlab parallelization
    if idx_job == 0
        parfor idx_job = 1:size(config.combinations, 1)
            try
                cwt_pc_dbi_analysis(config, idx_job);
            catch e
                fprintf(1,'The identifier was:\n%s',e.identifier);
                fprintf(1,'There was an error! The message was:\n%s',e.message);
            end
        end
    % or call only one job (done in parallel on cluster)
    else
        try
            cwt_pc_dbi_analysis(config, idx_job);
        catch e
            fprintf(1,'The identifier was:\n%s',e.identifier);
            fprintf(1,'There was an error! The message was:\n%s',e.message);
        end
    end
end

%%

function config = get_combinations(config)
% Create grids and reshape the grids to get all combinations (each job does
% one combination)
    
    config.combination_names = {'task', 'task_split', 'subject', 'with_rectify_emg', 'with_csd', 'with_low_freq_cutoff', 'with_high_freq_cutoff', 'with_brain_ics'};
    
    if config.do_cwt
        [taskGrid, splitGrid, subjectGrid, rectifyGrid, csdGrid, lfGrid, hfGrid, icGrid] ...
            = ndgrid(config.tasks, config.task_splits, config.subjects, config.with_rectify_emgs, config.with_csds, ...
            config.with_low_freq_cutoffs, config.with_high_freq_cutoffs, config.with_brain_icss);
    
        config.combinations = [taskGrid(:), splitGrid(:), subjectGrid(:), rectifyGrid(:), csdGrid(:) ...
            lfGrid(:) hfGrid(:) icGrid(:)];
    end
    
    if config.do_pc && ~config.do_dbi
        [taskGrid, splitGrid, subjectGrid, rectifyGrid, csdGrid, lfGrid, hfGrid, icGrid, nGrid, mGrid] ...
            = ndgrid(config.tasks, config.task_splits, config.subjects, config.with_rectify_emgs, config.with_csds, ...
            config.with_low_freq_cutoffs, config.with_high_freq_cutoffs, config.with_brain_icss, config.pc_ns, config.pc_ms);
    
        config.combination_names = [config.combination_names, {'pc_n', 'pc_m'}];
        config.combinations = [taskGrid(:), splitGrid(:), subjectGrid(:), rectifyGrid(:), csdGrid(:) ...
            lfGrid(:) hfGrid(:) icGrid(:) nGrid(:) mGrid(:)];
    
    end
    
    if config.do_dbi && ~config.do_pc
        [taskGrid, splitGrid, subjectGrid, rectifyGrid, csdGrid, lfGrid, hfGrid, icGrid, winsizeGrid] ...
            = ndgrid(config.tasks, config.task_splits, config.subjects, config.with_rectify_emgs, config.with_csds, ...
            config.with_low_freq_cutoffs, config.with_high_freq_cutoffs, config.with_brain_icss, ...
            config.dbi_win_sizes);
    
        config.combination_names = [config.combination_names, {'dbi_win_size'}];
        config.combinations = [taskGrid(:), splitGrid(:), subjectGrid(:), rectifyGrid(:), csdGrid(:) ...
            lfGrid(:), hfGrid(:), icGrid(:), winsizeGrid(:)];
    end
    
    if config.do_dbi && config.do_pc
        [taskGrid, splitGrid, subjectGrid, rectifyGrid, csdGrid, lfGrid, hfGrid, icGrid, nGrid, mGrid, winsizeGrid] ...
            = ndgrid(config.tasks, config.task_splits, config.subjects, config.with_rectify_emgs, config.with_csds, ...
            config.with_low_freq_cutoffs, config.with_high_freq_cutoffs, config.with_brain_icss, ...
            config.pc_ns, config.pc_ms, config.dbi_win_sizes);
    
        config.combination_names = [config.combination_names, {'pc_n', 'pc_m', 'dbi_win_size'}];
        config.combinations = [taskGrid(:), splitGrid(:), subjectGrid(:), rectifyGrid(:), csdGrid(:) ...
            lfGrid(:), hfGrid(:), icGrid(:), nGrid(:), mGrid(:), winsizeGrid(:)];
    end
   
end

