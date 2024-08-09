function f_check_ica_comps_bytasks()
eeglab nogui
path_data = 'D:\Experiments\corticomuscular_analysis\cmcpd\data\preproc\IZO';
config.subjects = ["PDH12", "PDP06"];
%config.subjects = ["PDH04", "PDH05", "PDH06", "PDH07", "PDH09", "PDH10", "PDH12", ...
%    "PDH14", "PDH18", "PDH19", "PDH20", "PDH22", "PDH23", "PDH24", "PDH25", "PDH26",  ...
%    "PDP02", "PDP03", "PDP04", "PDP06", "PDP09", "PDP10", "PDP11", "PDP13", "PDP15", ...
%    "PDP17", "PDP18", "PDP19", "PDP20", "PDP21", "PDP22"];

% config.task_splits = {'SL', 'SR', 'DL', 'DR', 'C'};
config.task_splits = {'SL'};
config.save_brain_ICs = 1;
config.do_manual_check = 1;
config.plot_ic_comps = 1;

for itask = 1:length(config.task_splits)
    task_split = config.task_splits{itask};
    kept_comps = [];
    for isubject = 1:length(config.subjects)
        subject = config.subjects{isubject};
        disp(['Subject' subject])
        if exist([path_data '\' subject '\' subject '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_split' task_split '_mc_icacl_icatdic.mat'], 'file')
            disp('Already exist')
            continue
        else
            disp([  subject ' task: ' task_split ' | Estimating components'' type using ICLabel package ...']);
        end

        try
            load([path_data '\' subject '\' subject '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_split' task_split '_mc_icacl_icatd.mat']);
            EEMG = eeg_checkset(EEMG);
        catch
            disp(['File doesnt exist for ' subject ' task: ' task_split])
            continue
        end

        % EEMG = do_dipfit(EEMG);
        [EEMG, kept_comps] = decide_on_components(EEMG, config, kept_comps);

        % save the EEMG struct
        save([path_data '\' subject '\' subject '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_split' task_split '_mc_icacl_icatdic.mat'], 'EEMG')

        close all;

    end
end
end

%%

function EEMG = do_dipfit(EEMG)
n_eeg_chans = 128;
disp('Estimating single equivalent current dipoles...')
dipfitdefs;

[newlocs, transform] = coregister(EEMG.chanlocs(1:n_eeg_chans), template_models(2).chanfile, ...
    'warp', 'auto', 'manual', 'off');

EEMG = pop_dipfit_settings( EEMG, ...
    'hdmfile',['\\plugins\\dipfit4.2\\standard_BEM\\standard_vol.mat'],...
    'coordformat','MNI',...
    'mrifile',['\\plugins\\dipfit4.2\\standard_BEM\\standard_mri.mat'],...
    'chanfile',['\\plugins\\dipfit4.2\\standard_BEM\\elec\\standard_1005.elc'],...
    'coord_transform',transform ,...
    'chansel',[1:n_eeg_chans] );

EEMG.etc.dipfit.warping_channel_names = [];
EEMG.etc.dipfit.RV_threshold = 100;
EEMG.etc.dipfit.remove_outside_head = 'off';
EEMG.etc.dipfit.number_of_dipoles = 1;

EEMG = pop_multifit(EEMG, [1:size(EEMG.icaweights,1)] , ...
    'threshold', EEMG.etc.dipfit.RV_threshold ,...
    'dipoles', EEMG.etc.dipfit.number_of_dipoles, ...
    'rmout', EEMG.etc.dipfit.remove_outside_head);
end

%%

function [EEMG, kept_comps] = decide_on_components(EEMG, config, kept_comps)

if ~isfield(EEMG.etc, 'spatial_filter'); EEMG.etc.spatial_filter = []; end

[EEMG, fig_handle] = ic_labeling(EEMG, config,  kept_comps);

% brain components to save
brain_comps = [];
for comp = 1:length(EEMG.etc.ic_classification.ICLabel.classifications)
    [maxVal, maxLab] = max(EEMG.etc.ic_classification.ICLabel.classifications(comp, :));
    if  maxLab == 1 && maxVal > 0.5
        brain_comps(end+1) = comp;
    end
end
disp(['Brain components: ' string(brain_comps)])

if isempty(kept_comps)
    additional_brain_comps = input('Add more brain comps, e.g. [1 2]: ');
    brain_comps = unique([brain_comps, additional_brain_comps]);
    kept_comps = unique(brain_comps);
else
    brain_comps = kept_comps; 
    % brain_comps = unique([brain_comps, kept_comps]);
end
disp(['Final brain components: ' string(brain_comps)])

EEMG.etc.spatial_filter.brain_comps = sort(brain_comps);

% get component activations
% eeg_brain_ics = eeg_getdatact(EEMG, 'component', brain_comps);

% keep only brain components and project them
EEMG = pop_subcomp( EEMG, brain_comps, 0, 1);

for ifig = 1:size(fig_handle, 2)
    if isvalid(fig_handle(ifig)), close(fig_handle(ifig)); end
end

end


%%

function [EEG, fig_handle] = ic_labeling(EEG, config, manual_rej_comps)

disp('Estimating components'' type using ICLabel package ...');

EEG = iclabel(EEG, 'lite');

%  Plot ICA components and optionally save figures
if config.plot_ic_comps && isempty(manual_rej_comps)
    igroupStart = 1;
    numComps = size(EEG.icawinv,2);
    for igroup = 1 : ceil(numComps/35)
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

%%
% OUTEEG = eeg_compatlas(EEG);
