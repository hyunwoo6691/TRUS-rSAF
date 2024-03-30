snr_25mm_CON_mean = mean(10*log10(fele_25mm(:,:,1)));
snr_25mm_CON_std = std(10*log10(fele_25mm(:,:,1)),0,1);

snr_5mm_CON_mean = mean(10*log10(fele_5mm(:,:,1)));
snr_5mm_CON_std = std(10*log10(fele_5mm(:,:,1)),0,1);


num_syn = [3 5 7 9 11 13 15 17 19 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65];

syn_ = 3;
idx = find(num_syn == syn_);
snr_5mm_rSAF3_mean = mean(10*log10(fele_5mm(:,:,idx+1)));
snr_5mm_rSAF3_std = std(10*log10(fele_5mm(:,:,idx+1)),0,1);

syn_ = 25;
idx = find(num_syn == syn_);
snr_5mm_rSAF25_mean = mean(10*log10(fele_5mm(:,:,idx+1)));
snr_5mm_rSAF25_std = std(10*log10(fele_5mm(:,:,idx+1)),0,1);

syn_ = 41;
idx = find(num_syn == syn_);
snr_5mm_rSAF41_mean = mean(10*log10(fele_5mm(:,:,idx+1)));
snr_5mm_rSAF41_std = std(10*log10(fele_5mm(:,:,idx+1)),0,1);

figure;
errorbar((1:7)*10, snr_25mm_CON_mean, snr_25mm_CON_std, 'LineWidth', 3); hold on;
errorbar((1:7)*10, snr_5mm_CON_mean, snr_5mm_CON_std, 'LineWidth', 3);
errorbar((1:7)*10, snr_5mm_rSAF3_mean, snr_5mm_rSAF3_std, 'LineWidth', 3);
errorbar((1:7)*10, snr_5mm_rSAF25_mean, snr_5mm_rSAF25_std, 'LineWidth', 3);
errorbar((1:7)*10, snr_5mm_rSAF41_mean, snr_5mm_rSAF41_std, 'LineWidth', 3); hold off;
legend('CON(25mm)', 'CON(5mm)', 'rSAF(5mm)-3', 'rSAF(5mm)-25', 'rSAF(5mm)-41','location','northeast');