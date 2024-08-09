
addpath 'D:\Experiments\corticomuscular_analysis\code\cmcpd\preprocessing'
eeglab nogui
path_data = 'D:\Experiments\corticomuscular_analysis\cmcpd\data\preproc\IZO';

config.task_splits = {'SL', 'DL', 'SR', 'DR', 'C'};
config.subjects = ["PDH12", "PDP06"];
% config.subjects = ["PDH04", "PDH05", "PDH06", "PDH07", "PDH08", "PDH09", "PDH10", "PDH12", "PDH14", "PDH18", "PDH19", "PDH20", "PDH22", "PDH23", "PDH24", "PDH25", "PDH26",  ...
%                   "PDP02", "PDP03", "PDP04", "PDP06",  "PDP09", "PDP10", "PDP11", "PDP13", "PDP15", "PDP17", "PDP18", "PDP19", "PDP20", "PDP21", "PDP22"];
n_eeg_chans = 128;

for itask = 1:length(config.task_splits)
    task_split = config.task_splits{itask};
    for isubject = 1:length(config.subjects)
        subject = config.subjects{isubject};

        if exist([path_data '\' subject '\' subject '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_split' task_split '_mc_icacl_icatd.mat'], 'file')
            continue
        end

        try
            load([path_data '\' subject '\' subject '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_split' task_split '_mc_icacl_icat.mat']);
            EEMG = eeg_checkset(EEMG);
        catch
            disp('File doesnt exist. Skipping')
            continue
        end

        % DipFit
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

        save([path_data '\' subject '\' subject '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_split' task_split '_mc_icacl_icatd.mat'], "EEMG");

    end
end