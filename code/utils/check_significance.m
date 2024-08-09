
load('D:\Experiments\corticomuscular_analysis\analysis\dbi\v9_allchans\dbis_results_cc.mat')
term_idx = 38;
eeg_channel = 'avgcen';
freq_comb = 'eeg_beta_low_emg_beta_low';
dbi_h_emgl_freqcomb_ichan = dbis_results.S.([eeg_channel '_emgl']).(freq_comb).H.cc_avg(:, term_idx);
dbi_p_emgl_freqcomb_ichan = dbis_results.S.([eeg_channel '_emgl']).(freq_comb).P.cc_avg(:, term_idx);
dbi_h_emgr_freqcomb_ichan = dbis_results.S.([eeg_channel '_emgr']).(freq_comb).H.cc_avg(:, term_idx);
dbi_p_emgr_freqcomb_ichan = dbis_results.S.([eeg_channel '_emgr']).(freq_comb).P.cc_avg(:, term_idx);

if add_abs
    ranksum(abs(dbi_h_emgr_freqcomb_ichan), abs(dbi_p_emgr_freqcomb_ichan))
else
    ranksum(dbi_h_emgr_freqcomb_ichan, dbi_p_emgr_freqcomb_ichan)
end

dbi_h_emgl_freqcomb_ichan = dbis_results.D.([eeg_channel '_emgl']).(freq_comb).H.cc_avg(:, term_idx);
dbi_p_emgl_freqcomb_ichan = dbis_results.D.([eeg_channel '_emgl']).(freq_comb).P.cc_avg(:, term_idx);
dbi_h_emgr_freqcomb_ichan = dbis_results.D.([eeg_channel '_emgr']).(freq_comb).H.cc_avg(:, term_idx);
dbi_p_emgr_freqcomb_ichan = dbis_results.D.([eeg_channel '_emgr']).(freq_comb).P.cc_avg(:, term_idx);

if add_abs
    ranksum(abs(dbi_h_emgr_freqcomb_ichan), abs(dbi_p_emgr_freqcomb_ichan))
else
    ranksum(dbi_h_emgr_freqcomb_ichan, dbi_p_emgr_freqcomb_ichan)
end
