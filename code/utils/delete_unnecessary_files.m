% delete unnecessary files 
path_root = 'D:\Experiments\corticomuscular_analysis\data\real\timefreq';
data_folders = dir([path_root filesep 'IZO*']);

for ifolder = 1:length(data_folders)
    subjects_folders = dir([path_root filesep data_folders(ifolder).name filesep 'P*']);
    for isubject = 1:length(subjects_folders)
        mat_files = dir([subjects_folders(isubject).folder filesep subjects_folders(isubject).name filesep '*.mat']);
        for ifile = 1:length(mat_files)
            filename = mat_files(ifile).name;
            filename_split = strsplit(filename, '_');
            if ~(contains(filename_split{end}, 'EMuovi') || contains(filename_split{end}, 'C'))
                disp(['Deleting: ' filename])
                delete([mat_files(ifile).folder filesep filename])
            end
            
        end
    end
    
end
