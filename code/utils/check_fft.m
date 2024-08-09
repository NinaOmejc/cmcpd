clear all; close all;
path_data = 'D:\Experiments\corticomuscular_analysis\data\real\preproc\IZO_dtrnd_lf0_hf0_onlytask_norect';
fname_data = 'PDH20_IZO_eemg_preproc_dtrnd_lf0_hf0_onlytask_norect';
load([path_data filesep fname_data])

emg = EEMG.data([1:128], :); % emg: [134:165 170:201

for ie = 1:2:size(emg, 1)
    ch_label = EEMG.chanlocs(ie).labels;
    S1 = emg(ie, :);
    Fs1 = 2000;
    Fs2 = 500;
    S2 = resampl_flow(S1, Fs1, Fs2); 
    L2 = length(S2);
    
    % Compute the Fourier transform for both signals
    Y2 = fft(S2, L2);
    
    % Compute the two-sided spectrum for both signals
    P2_2 = abs(Y2/L2);
    P1_2 = P2_2(1:L2/2+1);
    P1_2(2:end-1) = 2*P1_2(2:end-1);
    f2 = Fs2*(0:(L2/2))/L2;
    
    % Extend frequency range for plotting
    f2_extended = [f2, (Fs2+Fs2/2):Fs2/2:Fs2*10]; % Extend up to 10 times the sampling frequency
    P1_2_extended = [P1_2, zeros(1, length(f2_extended)-length(f2))]; % Pad with zeros
    
    % Plot the single-sided amplitude spectrum
    figure;
    hold on;
    plot(f2_extended,P1_2_extended,'LineWidth',1.5);
    hold off;
    title([num2str(ie) '-' ch_label]);
    xlabel('Frequency (Hz)');
    ylabel('|Y(f)|');
    legend('Fs = 500 Hz','Location','best');
    xlim([0 150]);
end

%%

% load('D:\Experiments\corticomuscular_analysis\data\real\preproc\IZO\PDH14\PDH14_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_splitDR_mc_ica.mat')
load('D:\Experiments\corticomuscular_analysis\data\real\preproc\IZO\PDH14\PDH14_IZO_eemg_onlytask_fs300_dtrnd_lf1_hf150_interp_splitDR_mc_icacl_icatdic.mat')
data = EEMG.data(18, :);
fs = EEMG.srate;
L = length(data);
Y = fft(data, L);

P2 = abs(Y/L); % Compute the two-sided spectrum P2
P1 = P2(1:L/2+1); % Compute the single-sided spectrum P1
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L; % Define the frequency domain f


% Plot the single-sided amplitude spectrum
figure;
plot(f(280:end), P1(280:end));
title('Single-Sided Amplitude Spectrum of X(t)');
xlabel('f (Hz)');
ylabel('|P1(f)|');

