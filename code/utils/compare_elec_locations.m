


x_orig = cell2mat({EEG.chanlocs.X}') / 1000;
y_orig = cell2mat({EEG.chanlocs.Y}') / 1000;
z_orig = cell2mat({EEG.chanlocs.Z}') / 1000;



biosemi128 = readtable('D:\Experiments\corticomuscular_analysis\electrodes\channel_BioSemi_128_A1.csv');

figure, 
plot3(biosemi128.x, biosemi128.y, biosemi128.z, 'r.')
hold on;
plot3(x_orig, y_orig, z_orig, 'b.')


