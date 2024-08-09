clear all;
subject = 'PDH12';
DL = load(['D:\Experiments\corticomuscular_analysis\data\real\preproc\IZO\' subject '\' subject '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_splitDL_mc_icacl_icatdic.mat']);
SL = load(['D:\Experiments\corticomuscular_analysis\data\real\preproc\IZO\' subject '\' subject '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_splitSL_mc_icacl_icatdic.mat']);
DR = load(['D:\Experiments\corticomuscular_analysis\data\real\preproc\IZO\' subject '\' subject '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_splitDR_mc_icacl_icatdic.mat']);
SR = load(['D:\Experiments\corticomuscular_analysis\data\real\preproc\IZO\' subject '\' subject '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_splitSR_mc_icacl_icatdic.mat']);
C = load(['D:\Experiments\corticomuscular_analysis\data\real\preproc\IZO\' subject '\' subject '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_splitC_mc_icacl_icatdic.mat']);

eeglab nogui;

for task = {'SL', 'DL', 'SR', 'DR', 'C'}
    data = eval(task{:});
    EEG = eeg_compatlas(data.EEMG, 'components', 1:size(data.EEMG.icaweights, 1));

    pop_viewprops(data.EEMG, 0, 1:size(data.EEMG.icaweights, 1));
    sgtitle(task{:})
end

%%

load('D:\Experiments\corticomuscular_analysis\data\real\STUDY\ic_table_short.mat')
ic_table_manual = ic_table;

% save('D:\Experiments\corticomuscular_analysis\data\real\STUDY\ic_table_short.mat', "ic_table")

%%


