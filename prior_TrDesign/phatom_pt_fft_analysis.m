spatial_filter = 'gauss';
%%
sigma_ = 1;
aSclineTheta = linspace(-0.5*40, 0.5*40, 86); % Ground truth transmitted angle
d_theta = abs(aSclineTheta(1)-aSclineTheta(2));

%%
stMID.nTGC_Atten = 0.5;                % [dB]

stMID.nDCRType = 'high';
stMID.nDCRTap = 128;                   % BPF tap #
stMID.nDCRFcut = 1e6;
stMID.nDCRF1 = stRFInfo.nFc - 5e6;     % [Hz]            if BANDPASS
stMID.nDCRF2 = stRFInfo.nFc + 5e6;     % [Hz]           if BANDPASS

stMID.nDemodFreq = stRFInfo.nFc;
stMID.nDemodTap = 128;


%%
tmp = rsaf_;
tmp(isnan(tmp)) = 0;

% [tmp_tgc, fil] = spatial_filtering(tmp, sigma_, d_theta, spatial_filter); % 'gauss', 'boxcar', 'none'
% tmpp = tmp;
tmp_tgc = tmp;
% [tmp_dcr, Fil] = DCR(tmpp, stMID, stRFInfo);
% [tmp_tgc, aTGCCurve] = fDTGC(tmp_dcr, stMID, stRFInfo, stBFInfo, size(tmp,1), stRFInfo.nUnitDis);
tmp_env = abs(hilbert(tmp_tgc));

[tmp_env, fil] = spatial_filtering(tmp_env, sigma_, d_theta, spatial_filter); % 'gauss', 'boxcar', 'none'

%
% tmppp = tmp_env;
% [aYAxis, aZAxis, tmp_dsc] = ScanConverter_convex(tmppp, dr, da, stBFInfo.nRadius, nHeight, nWidth, dz, dy);
%%
filter_ = fil;
conv_psf = conv_(1727,:);
rsaf_psf = rsaf_(1727,:);

convol_conv = conv(conv_psf, filter_, 'same');
convol_rsaf = conv(rsaf_psf, filter_, 'same');

figure(1);
plot(abs(hilbert(conv_psf)),'LineWidth', 2); hold on;
plot(abs(hilbert(convol_conv)),'LineWidth', 2); 
plot(conv(abs(hilbert(conv_psf)),filter_, 'same'), 'LineWidth', 2);
plot(abs(hilbert(rsaf_psf)),'LineWidth', 2);  hold off;
legend('conv', 'conv filter', 'filter after env', 'rsaf');

figure(4);
plot(abs(fftshift(fft(abs(hilbert(conv_psf))))),'LineWidth', 2); hold on;
plot(abs(fftshift(fft(abs(hilbert(convol_conv))))),'LineWidth', 2); 
plot(abs(fftshift(fft(conv(abs(hilbert(conv_psf)),filter_, 'same')))), 'LineWidth', 2); 
plot(abs(fftshift(fft(abs(hilbert(rsaf_psf))))), 'LineWidth', 2);hold off;
legend('conv', 'filter after bf', 'filter after env', 'rsaf');

%%
psf_conv_ = bfed_conv_(1728,:);
psf_conv_rf = bfed_conv_rf(1728,:);
psf_conv_env = bfed_conv_env(1728,:);

psf_rsaf_ = bfed_rsaf_(1728,:);
psf_rsaf_rf = bfed_rsaf_rf(1728,:);
psf_rsaf_env = bfed_rsaf_env(1728,:);

figure(100); 
plot(psf_conv_); hold on;
plot(psf_conv_rf);
plot(psf_conv_env); hold off;
legend('conv', 'gauss at rf', 'gauss at env');

figure(101); 
plot(psf_rsaf_); hold on;
plot(psf_rsaf_rf);
plot(psf_rsaf_env); hold off;
legend('rsaf', 'gauss at rf', 'gauss at env');

figure(102);
plot(psf_rsaf_/max(psf_rsaf_),'LineWidth', 2, 'color', 'r'); hold on;
plot(psf_conv_/max(psf_rsaf_),'LineWidth', 2, 'color', 'k');
plot(psf_conv_rf/max(psf_rsaf_),'LineWidth', 2, 'color', 'b');hold off;
legend('rsaf', 'conv', 'conv filter');

%% spatial fft
fft_fil_ = abs(fftshift(fft(fil,86)));
figure(202);
plot(fft_fil_);

fft_conv_ = abs(fftshift(fft(psf_conv_)));
fft_conv_rf = abs(fftshift(fft(psf_conv_rf)));
fft_conv_env = abs(fftshift(fft(psf_conv_env)));

fft_rsaf_ = abs(fftshift(fft(psf_rsaf_)));
fft_rsaf_rf = abs(fftshift(fft(psf_rsaf_rf)));
fft_rsaf_env = abs(fftshift(fft(psf_rsaf_env)));

figure(200); 
% plot(fft_conv_/max(fft_conv_)); hold on;
% plot(fft_conv_rf/max(fft_conv_rf));
% plot(fft_conv_env/max(fft_conv_env)); hold off;
plot(fft_conv_); hold on;
plot(fft_conv_rf);
plot(fft_conv_env); hold off;
legend('conv', 'gauss at rf', 'gauss at env');

figure(201); 
% plot(fft_rsaf_/max(fft_rsaf_)); hold on;
% plot(fft_rsaf_rf/max(fft_rsaf_rf));
% plot(fft_rsaf_env/max(fft_rsaf_env)); hold off;
plot(fft_rsaf_); hold on;
plot(fft_rsaf_rf);
plot(fft_rsaf_env); hold off;
legend('rsaf', 'gauss at rf', 'gauss at env');

%%
figure(202);
plot(fft_rsaf_/max(fft_rsaf_),'LineWidth', 2, 'color', 'k');hold on;
plot(fft_fil_,'LineWidth', 2); 
plot(fft_conv_/max(fft_rsaf_),'LineWidth', 2, 'color', 'b');
plot(fft_conv_rf/max(fft_rsaf_),'LineWidth', 2, 'color', 'r'); hold  off;
legend('rsaf', 'filter', 'conv', 'conv filter');
%%
conv_d = abs(fftshift(fft(psf_conv_)));
conv_filt_d = abs(fftshift(fft(conv(psf_conv_,fil/sum(fil),'same'))));
conv_filt_dd = abs(fftshift(fft(conv(psf_conv_,fil,'same'))));
rsaf_d = abs(fftshift(fft(psf_rsaf_)));

figure;
plot(conv_d/max(conv_d)); hold on;
plot(conv_filt_d/max(conv_d));
plot(rsaf_d/max(conv_d)); hold off;
legend('conv', 'conv filternorm', 'rsaf');

%%
figure; plot(fft_fil_); hold on;
plot(conv_d/max(conv_d)); hold off;
%%
fft2_conv_ = abs(fftshift(fft2(bfed_conv_)));
figure; imagesc(fft2_conv_/max(fft2_conv_(:)));
caxis([0 1]);

fft2_conv_rf = abs(fftshift(fft2(bfed_conv_rf)));
figure; imagesc(fft2_conv_rf/max(fft2_conv_(:)));
caxis([0 1]);

fft2_conv_env = abs(fftshift(fft2(bfed_conv_env)));
figure; imagesc(fft2_conv_env/max(fft2_conv_(:)));
caxis([0 1]);

%%
[conv_filtered, fil] = spatial_filtering(conv_, sigma_, d_theta, 'gauss'); % 'gauss', 'boxcar', 'none'
figure; 
subplot(1,3,1); imagesc(db(abs(hilbert(conv_))/max(max(abs(hilbert(rsaf_)))))); title('conv'); caxis([-80 0]);
subplot(1,3,2); imagesc(db(abs(hilbert(conv_filtered))/max(max(abs(hilbert(rsaf_)))))); title('conv filtered'); caxis([-80 0]); 
subplot(1,3,3); imagesc(db(abs(hilbert(rsaf_))/max(max(abs(hilbert(rsaf_)))))); title('rsaf'); caxis([-80 0]);
