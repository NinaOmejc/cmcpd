path_pc = 'D:\Experiments\corticomuscular_analysis\data\real\pc';
dirs1 = dir([path_pc filesep 'IZO*']);
for idir = 1:length(dirs1)
    path_dir1 = [dirs1(idir).folder filesep dirs1(idir).name];
    dirs2 = dir([path_dir1 filesep 'P*']);
    for jdir = 1:length(dirs2)
        files = dir([dirs2(jdir).folder filesep dirs2(jdir).name filesep '*.mat']);
        for ifile = 1:length(files)
            if files(ifile).bytes < 1
                delete([files(ifile).folder filesep files(ifile).name])
            end
        end
    end
end