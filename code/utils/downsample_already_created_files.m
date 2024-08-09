path_cwt = 'D:\Experiments\corticomuscular_analysis\data\real\cwt';
dirs1 = dir([path_cwt filesep 'IZO*']);
for idir = 1:length(dirs1)
    path_dir1 = [dirs1(idir).folder filesep dirs1(idir).name];
    dirs2 = dir([path_dir1 filesep 'P*']);
    for jdir = 1:length(dirs2)
        disp([path_dir1(end-8:end) '  ' dirs2(jdir).name])
        files = dir([dirs2(jdir).folder filesep dirs2(jdir).name filesep 'P*']);
        for ifile = 1:length(files)
            if contains(files(ifile).name, 'png')
                delete([files(ifile).folder filesep files(ifile).name])
            elseif contains(files(ifile).name, 'mat') && files(ifile).bytes > 50008597
                load([files(ifile).folder filesep files(ifile).name]);
                if ~isreal(WT)
                    WT = abs(WT.^2);
                    WT = single(downsample(WT', 6)');
                    save([files(ifile).folder filesep files(ifile).name], "WT", "f")
                end
                clear WT f
            end
        end
    end
end

