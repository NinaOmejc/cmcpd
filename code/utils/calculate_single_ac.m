load('D:\Experiments\corticomuscular_analysis\data\real\preproc\IZO\PDH04\PDH04_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_splitDL_mc_icacl_icatdic.mat');

sig1 = EEMG.data(18, :);
sig2 = EEMG.data(129, :);
srate = EEMG.srate;

[WT1,freq, wopt] = wt_nina(sig1, srate, ...
                    'fmin', 4, 'fmax', 90, ...
                    'Preprocess', 'off', 'Wavelet', 'Morlet', 'f0', 1, ...
                    'Plot', 'off', 'CutEdges', 'on', 'nv', 30, ...
                    'event_tbl', [], 'Padding', 'none', 'fig_title', 'sig1');

[WT2] = wt_nina(sig2, srate, ...
                'fmin', 4, 'fmax', 90, ...
                'Preprocess', 'off', 'Wavelet', 'Morlet', 'f0', 1, ...
                'Plot', 'off', 'CutEdges', 'on', 'nv', 30, ...
                'event_tbl', [], 'Padding', 'none', 'fig_title', 'sig2');

[TPC, PC] = tlphcoh_nina(WT1, WT2, freq, 1, 1, srate, 1, [], [], '', 'on'); 


%%

% Parameters
srate = 1000;          % Sampling frequency
t = 0:1/srate:2-1/srate;  % Time vector
f1 = 10;            % Frequency of the first signal
f2 = 10;            % Frequency of the second signal
A = 1;              % Amplitude of the signals
noise_power = 0.2;  % Power of the additive noise

% Generate example signals
x1 = A*sin(2*pi*f1*t) + noise_power*randn(size(t));
x2 = A*sin(2*pi*f2*t) + noise_power*randn(size(t));

% Plot example signals
figure;
subplot(2,1,1);
plot(t, x1);
title('Signal 1');
xlabel('Time (s)');
ylabel('Amplitude');
subplot(2,1,2);
plot(t, x2);
title('Signal 2');
xlabel('Time (s)');
ylabel('Amplitude');

[WT1,freq] = wt_nina(x1, srate, ...
                    'fmin', 4, 'fmax', 90, ...
                    'Preprocess', 'off', 'Wavelet', 'Morlet', 'f0',  1, ...
                    'Plot', 'on', 'CutEdges', 'on', 'nv', 30, ...
                    'event_tbl', [], 'Padding', 'none', 'fig_title', 'sig1');

[WT2] = wt_nina(x2, srate, ...
                'fmin', 4, 'fmax', 90, ...
                'Preprocess', 'off', 'Wavelet', 'Morlet', 'f0', 1, ...
                'Plot', 'on', 'CutEdges', 'on', 'nv', 30, ...
                'event_tbl', [], 'Padding', 'none', 'fig_title', 'sig2');

[TPC, PC] = tlphcoh_nina(WT1, WT2, freq, 1, 1, srate, 1, [], [], '', 'on'); 

[TAC, AC] = tlampcoh_nina(WT1, WT2, freq, 1, 1, srate, 1, [], [], '', 'on'); 

wcoherence(x1, x2, srate);



% Perform wavelet transform on both signals
[wt_x1, f_x1] = cwt(x1, 'amor', Fs);
[wt_x2, f_x2] = cwt(x2, 'amor', Fs);

% Compute coherence
coherence = zeros(size(wt_x1));
for i = 1:size(wt_x1, 1)
    % Extract wavelet coefficients
    seg_x1 = wt_x1(i, :);
    seg_x2 = wt_x2(i, :);
    
    % Compute coherence
    CSD = conj(seg_x1) .* seg_x2;
    ASD_x1 = abs(seg_x1).^2;
    ASD_x2 = abs(seg_x2).^2;
    coherence(i, :) = abs(CSD).^2 ./ (ASD_x1 .* ASD_x2);
end

% Average coherence across segments (if necessary)
average_coherence = mean(coherence, 2);

% Plot coherence
figure;
plot(f_x1, average_coherence);
title('Coherence between Signals');
xlabel('Frequency (Hz)');
ylabel('Coherence');

