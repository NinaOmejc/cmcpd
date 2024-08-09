eeglab nogui 

path_data = 'D:\Experiments\corticomuscular_analysis\data\real\raw';
files = dir([path_data '\P*']);
new_fs = 500;

for ifile = 1:length(files)
    disp(['File : ' num2str(ifile)])
    if contains(files(ifile).name, 'onlydata')
        disp('not right file')
        continue
    end

    if exist([path_data filesep files(ifile).name(1:end-4) '_fs' num2str(new_fs) '_onlydata.mat'], 'file')
        disp('already done.')
        continue
    end

    loaded_variable = whos('-file', [path_data filesep files(ifile).name]);
    load([path_data filesep files(ifile).name])
    EEMG = eval(loaded_variable.name);
    EEMG = eeg_checkset(EEMG);
    EEMG = pop_resample(EEMG, new_fs);
    data = EEMG.data;
    save([path_data filesep files(ifile).name(1:end-4) '_fs' num2str(new_fs) '_onlydata.mat'], "data")

end