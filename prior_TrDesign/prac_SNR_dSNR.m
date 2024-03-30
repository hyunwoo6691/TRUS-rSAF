snr_REF = [104.8777818	91.3153037	74.8287561	60.36863128	47.91535036	37.43309533	31.62385683];

snr_6871 = [102.0758289	82.61848312	68.09641303	55.9742821	46.02091469	37.44129733	29.96730439];
snr_4724 = [106.2088466	85.23540699	70.9163758	59.35092596	49.56799539	40.18132826	33.78575628];
snr_3779 = [108.4303526	87.73063749	73.32009986	61.76613565	51.41225474	42.3453022	35.55494842];
snr_2362 = [112.0773217	91.29362208	76.80215687	65.21777971	54.31873095	45.10444252	39.64659234];
snr_1989 = [112.8998028	91.92386131	77.72562513	66.2060756	55.20025047	46.43828723	39.62734287];

%% SNR
axis_ = [10 20 30 40 50 60 70]*1e-3;
figure(1);
plot(axis_*1e3, snr_REF,'LineWidth',2,'color','k'); hold on;
plot(axis_*1e3, snr_6871,'LineWidth',2);
plot(axis_*1e3, snr_4724,'LineWidth',2);
plot(axis_*1e3, snr_3779,'LineWidth',2);
plot(axis_*1e3, snr_2362,'LineWidth',2);
plot(axis_*1e3, snr_1989,'LineWidth',2); hold off;
legend('TRUS-REF','6871','4724','3779','2362','1989');

%% Delta_SNR
dsnr_6871 = snr_6871 - snr_REF;
dsnr_4724 = snr_4724 - snr_REF;
dsnr_3779 = snr_3779 - snr_REF;
dsnr_2362 = snr_2362 - snr_REF;
dsnr_1989 = snr_1989 - snr_REF;

figure(2);
plot(axis_*1e3, dsnr_6871,'LineWidth',2); hold on;
plot(axis_*1e3, dsnr_4724,'LineWidth',2);
plot(axis_*1e3, dsnr_3779,'LineWidth',2);
plot(axis_*1e3, dsnr_2362,'LineWidth',2);
plot(axis_*1e3, dsnr_1989,'LineWidth',2); 
line([axis_(1) axis_(end)]*1e3, [0 0],'LineStyle','--','color','k');hold off;
legend('6871','4724','3779','2362','1989');