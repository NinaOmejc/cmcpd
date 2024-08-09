% This code is used to remove single-subject artefactual components from
% data. 

function c_check_ica_comps_bysubs()
clear all;
eeglab nogui
path_data = 'D:\Experiments\corticomuscular_analysis\cmcpd\data\preproc\IZO';
config.subjects = ["PDH12", "PDP06"];
config.task_splits = {'SL', 'DL', 'SR', 'DR', 'C'};
config.save_brain_ICs = 0;
config.do_manual_check = 1;
config.plot_ic_comps = 1;

for isubject = 1:length(config.subjects)
    subject = config.subjects{isubject};
    eeglab nogui;
    disp(['Subject' subject])
    manual_rej_comps = [];
    for itask = 1:length(config.task_splits)
        try
            task_split = config.task_splits{itask};
            load([path_data '\' subject '\' subject '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_split' task_split '_mc_ica.mat']);
            EEMG = eeg_checkset(EEMG);
        catch
            disp(['File doesnt exist for ' subject ' task: ' task_split])
            continue
        end

        if exist([path_data '\' subject '\' subject '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_split' task_split '_mc_icacl_icatcl.mat'], 'file')
            disp('Already exist')
            continue
        else
            disp([  subject ' task: ' task_split 'Estimating components'' type using ICLabel package ...']);
        end

        [EEMG, manual_rej_comps] = decide_on_components(EEMG, config, manual_rej_comps);
       
        % save the EEMG struct
        save([path_data '\' subject '\' subject '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_split' task_split '_mc_icacl.mat'], "EEMG");
        close all;

    end
end
end

%%

function [EEMG, manual_rej_comps] = decide_on_components(EEMG, config, manual_rej_comps)
    if ~isfield(EEMG.etc, 'spatial_filter'); EEMG.etc.spatial_filter = []; end
    if isfield(EEMG.etc.spatial_filter, 'rejected_amica_comps')
        % print if there were some components already rejected.
        disp('Currently rejected AMICA comps:')
        disp(sort(cell2mat(EEMG.etc.spatial_filter.rejected_amica_comps(:, 1)')))
    else

    [EEMG, fig_handle] = ic_labeling(EEMG, config, manual_rej_comps);

    % Reject components (all but brain category above 85% certainty)
    % for label names check "EEG.etc.ic_classification.ICLabel.classes"
    % {'1: Brain', '2: Muscle', '3: Eye', '4: Heart','5: Line Noise', '6: Channel Noise', '7: Other'}

    rejComps = [];
    rejCompsLabs = [];
    for comp = 1:length(EEMG.etc.ic_classification.ICLabel.classifications)
        [maxVal, maxLab] = max(EEMG.etc.ic_classification.ICLabel.classifications(comp, :));
        if  maxLab == 3 && maxVal > 0.5
            rejComps(end+1) = comp;
            rejCompsLabs(end+1) = maxLab;
        elseif any(maxLab == [2, 4, 5, 6]) && maxVal > 0.75
            rejComps(end+1) = comp;
            rejCompsLabs(end+1) = maxLab;
        end
    end

    % manually check the automatically rejected components
    disp([newline 'These components would be automatically rejected: ' num2str(rejComps) newline])
    doKeep = 1;
    while config.do_manual_check && doKeep
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

if isfield(EEMG.etc.spatial_filter, 'rejected_amica_comps')
    [EEMG, fig_handle] = ic_labeling(EEMG, config, []);
    rejComps = [];
    rejCompsLabs = [];
end

% plot components time series

% additionally manually reject components
if ~isempty(manual_rej_comps)
    rejComps = manual_rej_comps;
    rejCompsLabs = [];
    for ICtoRemove = 1:length(rejComps)
        [~, maxLab] = max(EEMG.etc.ic_classification.ICLabel.classifications(ICtoRemove, :));
        rejCompsLabs(end+1) = maxLab;
    end
else
    pop_eegplot(EEMG, 0, 1, 1);
    nextManualRej = 1;
    while config.do_manual_check && nextManualRej
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
                        [~, maxLab] = max(EEMG.etc.ic_classification.ICLabel.classifications(ICtoRemove, :));
                        rejCompsLabs(end+1) = maxLab;
    
                    catch
                        disp(['There was an error for the rejection of IC label' num2str(ICtoRemove) ', skipping.'])
                    end
                end
            end
            disp(['These are the updated components that will be rejected: ' num2str(rejComps)])
    
        end
    end
end

% add info about which ICs were removed to the EEG struct
disp(['Final rejected components: ' num2str(rejComps)])
if isfield(EEMG.etc.spatial_filter, 'rejected_amica_comps')
    n_already_removed_comps = length(EEMG.etc.spatial_filter.rejected_amica_comps);
    for irc = 1: length(rejComps)
        EEMG.etc.spatial_filter.rejected_amica_comps(n_already_removed_comps+irc, 1) = num2cell(rejComps(irc))';
        EEMG.etc.spatial_filter.rejected_amica_comps(n_already_removed_comps+irc, 2) = EEMG.etc.ic_classification.ICLabel.classes(rejCompsLabs(irc));
    end
else
    EEMG.etc.spatial_filter.rejected_amica_comps = num2cell(rejComps)';
    for irc = 1: length(rejCompsLabs)
        EEMG.etc.spatial_filter.rejected_amica_comps(irc, 2) = EEMG.etc.ic_classification.ICLabel.classes(rejCompsLabs(irc));
    end
end

if isempty(manual_rej_comps);     manual_rej_comps = rejComps; end;


% brain components to save
if config.save_brain_ICs

    brain_comps = [];
    for comp = 1:length(EEMG.etc.ic_classification.ICLabel.classifications)
        [maxVal, maxLab] = max(EEMG.etc.ic_classification.ICLabel.classifications(comp, :));
        if  maxLab == 1
            brain_comps(end+1) = comp;
        end
    end
    disp(['Brain components: ' string(brain_comps)])

    additional_brain_comps = input('Add more brain comps, e.g. [1 2]: ');
    brain_comps = [brain_comps, additional_brain_comps];
    EEMG.etc.spatial_filter.brain_comps = brain_comps;

    % get component activations
    eeg_brain_ics = eeg_getdatact(EEMG, 'component', brain_comps);

    save([config.path_data_preproc filesep config.fname_data_preproc_after_csd '_brain_ICs.mat'], 'eeg_brain_ics')
end


% reject components (they will also be removed from the data)
EEMG = pop_subcomp( EEMG, rejComps, 0, 0);
EEMG.etc.spatial_filter.rejection_completed = 1;

for ifig = 1:size(fig_handle, 2)
    if isvalid(fig_handle(ifig)), close(fig_handle(ifig)); end
end

   screen_size = get(0, 'ScreenSize');
   width = screen_size(3);
   height = screen_size(4);

   eegplot(EEMG.data, ...
            'srate', EEMG.srate, ...
            'events', EEMG.event, ...
            'winlength', 20, ...
            'submean', 'on', ...
            'spacing', 30, ...
            'dispchans', 64, ...
            'position', [round(0.05*width), round(0.075*height), width-round(0.10*width), height-round(0.15*height)], ...
            'title', '');

% keep_cleaning = input("Would you like to do another round of data inspection and cleaning? [1-yes, 0-no]");
n_cleans = n_cleans + 1;
end


%%

function [EEG, fig_handle] = ic_labeling(EEG, config, manual_rej_comps)
    
    disp('Estimating components'' type using ICLabel package ...');

    EEG = iclabel(EEG, 'lite');

    %  Plot ICA components and optionally save figures
    if config.plot_ic_comps && isempty(manual_rej_comps)
        igroupStart = 1;
        numComps = size(EEG.icawinv,2);
        for igroup = 1 :1 % ceil(numComps/35)
            if igroup == ceil(numComps/35)
                pop_viewprops(EEG, 0, igroupStart:numComps);
            else
                pop_viewprops(EEG, 0, igroupStart:igroup*35);
            end
            igroupStart = igroupStart+35;
            fig_handle(igroup) = gcf;
        end
    else
        fig_handle = [];
    end

end

