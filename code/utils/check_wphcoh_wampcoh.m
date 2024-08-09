clear all;
sub = 'PDP13';
task = 'SL';
load(['D:\Experiments\corticomuscular_analysis\data\real\preproc\IZO\' sub '\' sub '_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_split' task '_mc_icacl_icatdic.mat']);
fs = EEMG.srate;
fmin = 4; fmax = 90; f0 = 1; nv = 30;
S1 = EEMG.data(18, :);
S2 = EEMG.data(130, :);

[WT1, ~] = wt_nina(S1, fs, 'fmin', fmin, 'fmax', fmax, 'Preprocess', 'off', 'Wavelet', 'Morlet', ...
                   'f0', f0, 'Plot', 'pow+tm', 'CutEdges', 'on', 'nv', nv, 'visible', 'on', ...
                   'event_tbl', [], 'Padding', 'none', 'fig_title', 'WT1 of S1', 'Display', 'notify');

[WT2, f] = wt_nina(S2, fs, 'fmin', fmin, 'fmax', fmax, 'Preprocess', 'off', 'Wavelet', 'Morlet', ...
                   'f0', f0, 'Plot', 'pow+tm', 'CutEdges', 'on', 'nv', nv, 'visible', 'on', ...
                   'event_tbl', [], 'Padding', 'none', 'fig_title', 'WT2 of S2', 'Display', 'notify');

% phase coherence
[TPC, PC] = tlphcoh_nina(WT1, WT2, f, 1, 1, fs, 1, [], [], 'TPC', 'on'); 
figure, plot(PC)

% phase coherence from the paper
a1 = real(WT1);
b1 = imag(WT1);
a2 = real(WT2);
b2 = imag(WT2);

cosdiff = (a1 .* a2 + b1 .* b2) ./ ( sqrt(a1.^2 + b1.^2) .* sqrt(a2.^2 + b2.^2) );
sindiff = (b1 .* a2 - a1 .* b2) ./ ( sqrt(a1.^2 + b1.^2) .* sqrt(a2.^2 + b2.^2) );
PC_paper = sqrt( nanmean(cosdiff, 2).^2 + nanmean(sindiff, 2).^2 );

figure, hold on
plot(PC, 'LineWidth', 4)
plot(PC_paper, 'LineWidth', 2)


% amplitude coherence
% [ampwcoh, ~, f]  = wcoherence(S1, S2, fs, 'FrequencyLimits', [fmin fmax], 'VoicesPerOctave', nv, 'PhaseDisplayThreshold', 1);

