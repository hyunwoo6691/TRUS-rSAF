clc; close all; clear;

%%
num_syn = [3 5 7 9 11 13 15 17 19 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65];

depth_ = 5;
root_n = (num_syn);
figure(depth_); hold on;

rsaf_tmpM = mean(squeeze(fele_5mm(:,depth_,2:end))./fele_5mm(1,depth_,1));
rsaf_tmpStd = std(squeeze(fele_5mm(:,depth_,2:end))./fele_5mm(1,depth_,1),0,1);
errorbar(num_syn,rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(fele_10mm(:,depth_,2:end))./fele_10mm(1,depth_,1));
rsaf_tmpStd = std(squeeze(fele_10mm(:,depth_,2:end))./fele_10mm(1,depth_,1),0,1);
errorbar(num_syn,rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(fele_15mm(:,depth_,2:end))./fele_15mm(1,depth_,1));
rsaf_tmpStd = std(squeeze(fele_15mm(:,depth_,2:end))./fele_15mm(1,depth_,1),0,1);
errorbar(num_syn,rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(fele_20mm(:,depth_,2:end))./fele_20mm(1,depth_,1));
rsaf_tmpStd = std(squeeze(fele_20mm(:,depth_,2:end))./fele_20mm(1,depth_,1),0,1);
errorbar(num_syn,rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(fele_25mm(:,depth_,2:end))./fele_25mm(1,depth_,1));
rsaf_tmpStd = std(squeeze(fele_25mm(:,depth_,2:end))./fele_25mm(1,depth_,1),0,1);
errorbar(num_syn,rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);


plot(num_syn, root_n, 'LineWidth', 4, 'color', 'b'); hold off;

legend('fele = 5mm', 'fele = 10mm', 'fele = 15mm', 'fele = 20mm', 'fele = 25mm','ref');

%% f_elv
close all;
root_n = (num_syn);
% 5cm
depth_ = 3;
n_5mm = 65; idx_5mm = find(num_syn == n_5mm);
n_10mm = 33; idx_10mm = find(num_syn == n_10mm);
n_15mm = 17; idx_15mm = find(num_syn == n_15mm);
n_20mm = 9; idx_20mm = find(num_syn == n_20mm);
n_25mm = 5; idx_25mm = find(num_syn == n_25mm);

figure(depth_); hold on;

rsaf_tmpM = mean(squeeze(fele_5mm(:,depth_,2:idx_5mm+1)./fele_5mm(:,depth_,1)));
rsaf_tmpStd = std(squeeze(fele_5mm(:,depth_,2:idx_5mm+1)./fele_5mm(:,depth_,1)),0,1);
errorbar((num_syn(1:idx_5mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);
max_val = max(rsaf_tmpM);

rsaf_tmpM = mean(squeeze(fele_10mm(:,depth_,2:idx_10mm+1)./fele_10mm(1,depth_,1)));
rsaf_tmpStd = std(squeeze(fele_10mm(:,depth_,2:idx_10mm+1)./fele_10mm(1,depth_,1)),0,1);
errorbar((num_syn(1:idx_10mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(fele_15mm(:,depth_,2:idx_15mm+1)./fele_15mm(1,depth_,1)));
rsaf_tmpStd = std(squeeze(fele_15mm(:,depth_,2:idx_15mm+1)./fele_15mm(1,depth_,1)),0,1);
errorbar((num_syn(1:idx_15mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(fele_20mm(:,depth_,2:idx_20mm+1)./fele_20mm(1,depth_,1)));
rsaf_tmpStd = std(squeeze(fele_20mm(:,depth_,2:idx_20mm+1)./fele_20mm(1,depth_,1)),0,1);
errorbar((num_syn(1:idx_20mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(fele_25mm(:,depth_,2:idx_25mm+1)./fele_25mm(1,depth_,1)));
rsaf_tmpStd = std(squeeze(fele_25mm(:,depth_,2:idx_25mm+1)./fele_25mm(1,depth_,1)),0,1);
errorbar((num_syn(1:idx_25mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);

plot((num_syn), root_n, 'LineWidth', 4, 'color', 'b'); hold off;
ylim([0 1.2*max_val]);

legend('fele = 5mm', 'fele = 10mm', 'fele = 15mm', 'fele = 20mm', 'fele = 25mm','ref','location','southeast');

% 6cm
depth_ = 5;
n_5mm = 65; idx_5mm = find(num_syn == n_5mm);
n_10mm = 39; idx_10mm = find(num_syn == n_10mm);
n_15mm = 23; idx_15mm = find(num_syn == n_15mm);
n_20mm = 15; idx_20mm = find(num_syn == n_20mm);
n_25mm = 11; idx_25mm = find(num_syn == n_25mm);

figure(depth_); hold on;

rsaf_tmpM = mean(squeeze(fele_5mm(:,depth_,2:idx_5mm+1)./fele_5mm(1,depth_,1)));
rsaf_tmpStd = std(squeeze(fele_5mm(:,depth_,2:idx_5mm+1)./fele_5mm(1,depth_,1)),0,1);
errorbar((num_syn(1:idx_5mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);
max_val = max(rsaf_tmpM);

rsaf_tmpM = mean(squeeze(fele_10mm(:,depth_,2:idx_10mm+1)/fele_10mm(1,depth_,1)));
rsaf_tmpStd = std(squeeze(fele_10mm(:,depth_,2:idx_10mm+1)/fele_10mm(1,depth_,1)),0,1);
errorbar((num_syn(1:idx_10mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(fele_15mm(:,depth_,2:idx_15mm+1)/fele_15mm(1,depth_,1)));
rsaf_tmpStd = std(squeeze(fele_15mm(:,depth_,2:idx_15mm+1)/fele_15mm(1,depth_,1)),0,1);
errorbar((num_syn(1:idx_15mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(fele_20mm(:,depth_,2:idx_20mm+1)/fele_20mm(1,depth_,1)));
rsaf_tmpStd = std(squeeze(fele_20mm(:,depth_,2:idx_20mm+1)/fele_20mm(1,depth_,1)),0,1);
errorbar((num_syn(1:idx_20mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(fele_25mm(:,depth_,2:idx_25mm+1)/fele_25mm(1,depth_,1)));
rsaf_tmpStd = std(squeeze(fele_25mm(:,depth_,2:idx_25mm+1)/fele_25mm(1,depth_,1)),0,1);
errorbar((num_syn(1:idx_25mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);

plot((num_syn), root_n, 'LineWidth', 4, 'color', 'b'); hold off;
ylim([0 1.2*max_val]);

legend('fele = 5mm', 'fele = 10mm', 'fele = 15mm', 'fele = 20mm', 'fele = 25mm','ref','location','southeast');

% 7cm
depth_ = 7;
n_5mm = 65; idx_5mm = find(num_syn == n_5mm);
n_10mm = 45; idx_10mm = find(num_syn == n_10mm);
n_15mm = 27; idx_15mm = find(num_syn == n_15mm);
n_20mm = 19; idx_20mm = find(num_syn == n_20mm);
n_25mm = 15; idx_25mm = find(num_syn == n_25mm);

figure(depth_); hold on;

rsaf_tmpM = mean(squeeze(fele_5mm(:,depth_,2:idx_5mm+1)./fele_5mm(1,depth_,1)));
rsaf_tmpStd = std(squeeze(fele_5mm(:,depth_,2:idx_5mm+1)./fele_5mm(1,depth_,1)),0,1);
errorbar((num_syn(1:idx_5mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);
max_val = max(rsaf_tmpM);

rsaf_tmpM = mean(squeeze(fele_10mm(:,depth_,2:idx_10mm+1)/fele_10mm(1,depth_,1)));
rsaf_tmpStd = std(squeeze(fele_10mm(:,depth_,2:idx_10mm+1)/fele_10mm(1,depth_,1)),0,1);
errorbar((num_syn(1:idx_10mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(fele_15mm(:,depth_,2:idx_15mm+1)/fele_15mm(1,depth_,1)));
rsaf_tmpStd = std(squeeze(fele_15mm(:,depth_,2:idx_15mm+1)/fele_15mm(1,depth_,1)),0,1);
errorbar((num_syn(1:idx_15mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(fele_20mm(:,depth_,2:idx_20mm+1)/fele_20mm(1,depth_,1)));
rsaf_tmpStd = std(squeeze(fele_20mm(:,depth_,2:idx_20mm+1)/fele_20mm(1,depth_,1)),0,1);
errorbar((num_syn(1:idx_20mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(fele_25mm(:,depth_,2:idx_25mm+1)/fele_25mm(1,depth_,1)));
rsaf_tmpStd = std(squeeze(fele_25mm(:,depth_,2:idx_25mm+1)/fele_25mm(1,depth_,1)),0,1);
errorbar((num_syn(1:idx_25mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);

plot((num_syn), root_n, 'LineWidth', 4, 'color', 'b'); hold off;
ylim([0 1.2*max_val]);
legend('fele = 5mm', 'fele = 10mm', 'fele = 15mm', 'fele = 20mm', 'fele = 25mm','ref','location','southeast');



%% h
num_syn = [3 5 7 9 11 13 15 17 19 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65];

close all;
root_n = (num_syn);

h_3mm = snr_3mm;
h_4mm = snr_4mm;
h_5mm = snr_5mm;
h_6mm = snr_6mm;
h_7mm = snr_7mm;

% h_3mm = noise_3mm;
% h_4mm = noise_4mm;
% h_5mm = noise_5mm;
% h_6mm = noise_6mm;
% h_7mm = noise_7mm;

% h_3mm = signal_3mm;
% h_4mm = signal_4mm;
% h_5mm = signal_5mm;
% h_6mm = signal_6mm;
% h_7mm = signal_7mm;

% ref_ = 1./num_syn; % for noise
% ref_ = 1./(num_syn.^2); % for signal
ref_ = num_syn;

% 5cm
depth_ = 3;
n_3mm = 45; idx_3mm = find(num_syn == n_3mm);
n_4mm = 65; idx_4mm = find(num_syn == n_4mm);
n_5mm = 65; idx_5mm = find(num_syn == n_5mm);
n_6mm = 65; idx_6mm = find(num_syn == n_6mm);
n_7mm = 65; idx_7mm = find(num_syn == n_7mm);

figure(depth_); hold on;

rsaf_tmpM = mean(squeeze(h_3mm(:,depth_,2:idx_3mm+1)))./mean(squeeze(h_3mm(:,depth_,1)));
rsaf_tmpStd = std(squeeze(h_3mm(:,depth_,2:idx_3mm+1)./h_3mm(:,depth_,1)),0,1);
% errorbar((num_syn(1:idx_3mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);
plot(num_syn(1:idx_3mm), rsaf_tmpM, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(h_4mm(:,depth_,2:idx_4mm+1)))./mean(squeeze(h_4mm(:,depth_,1)));
rsaf_tmpStd = std(squeeze(h_4mm(:,depth_,2:idx_4mm+1)./h_4mm(1,depth_,1)),0,1);
% errorbar((num_syn(1:idx_4mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);
plot(num_syn(1:idx_4mm), rsaf_tmpM, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(h_5mm(:,depth_,2:idx_5mm+1)))./mean(squeeze(h_5mm(:,depth_,1)));
rsaf_tmpStd = std(squeeze(h_5mm(:,depth_,2:idx_5mm+1)./h_5mm(1,depth_,1)),0,1);
% errorbar((num_syn(1:idx_5mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);
plot(num_syn(1:idx_5mm), rsaf_tmpM, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(h_6mm(:,depth_,2:idx_6mm+1)))./mean(squeeze(h_6mm(:,depth_,1)));
rsaf_tmpStd = std(squeeze(h_6mm(:,depth_,2:idx_6mm+1)./h_6mm(1,depth_,1)),0,1);
% errorbar((num_syn(1:idx_6mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);
plot(num_syn(1:idx_6mm), rsaf_tmpM, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(h_7mm(:,depth_,2:idx_7mm+1)))./mean(squeeze(h_7mm(:,depth_,1)));
rsaf_tmpStd = std(squeeze(h_7mm(:,depth_,2:idx_7mm+1)./h_7mm(1,depth_,1)),0,1);
% errorbar((num_syn(1:idx_7mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);
plot(num_syn(1:idx_7mm), rsaf_tmpM, 'LineWidth', 3);
max_val = max(rsaf_tmpM);

plot((num_syn), ref_, 'LineWidth', 4, 'color', 'b'); hold off;
% ylim([0 1.2*max_val]);

legend('h = 3mm', 'h = 4mm', 'h = 5mm', 'h = 6mm', 'h = 7mm', 'ref','location','southeast');

% 6cm
depth_ = 5;
n_3mm = 53; idx_3mm = find(num_syn == n_3mm);
n_4mm = 65; idx_4mm = find(num_syn == n_4mm);
n_5mm = 65; idx_5mm = find(num_syn == n_5mm);
n_6mm = 65; idx_6mm = find(num_syn == n_6mm);
n_7mm = 65; idx_7mm = find(num_syn == n_7mm);

figure(depth_); hold on;

rsaf_tmpM = mean(squeeze(h_3mm(:,depth_,2:idx_3mm+1)))./mean(squeeze(h_3mm(:,depth_,1)));
rsaf_tmpStd = std(squeeze(h_3mm(:,depth_,2:idx_3mm+1)./h_3mm(:,depth_,1)),0,1);
% errorbar((num_syn(1:idx_3mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);
plot(num_syn(1:idx_3mm), rsaf_tmpM, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(h_4mm(:,depth_,2:idx_4mm+1)))./mean(squeeze(h_4mm(:,depth_,1)));
rsaf_tmpStd = std(squeeze(h_4mm(:,depth_,2:idx_4mm+1)./h_4mm(1,depth_,1)),0,1);
% errorbar((num_syn(1:idx_4mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);
plot(num_syn(1:idx_4mm), rsaf_tmpM, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(h_5mm(:,depth_,2:idx_5mm+1)))./mean(squeeze(h_5mm(:,depth_,1)));
rsaf_tmpStd = std(squeeze(h_5mm(:,depth_,2:idx_5mm+1)./h_5mm(1,depth_,1)),0,1);
% errorbar((num_syn(1:idx_5mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);
plot(num_syn(1:idx_5mm), rsaf_tmpM, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(h_6mm(:,depth_,2:idx_6mm+1)))./mean(squeeze(h_6mm(:,depth_,1)));
rsaf_tmpStd = std(squeeze(h_6mm(:,depth_,2:idx_6mm+1)./h_6mm(1,depth_,1)),0,1);
% errorbar((num_syn(1:idx_6mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);
plot(num_syn(1:idx_6mm), rsaf_tmpM, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(h_7mm(:,depth_,2:idx_7mm+1)))./mean(squeeze(h_7mm(:,depth_,1)));
rsaf_tmpStd = std(squeeze(h_7mm(:,depth_,2:idx_7mm+1)./h_7mm(1,depth_,1)),0,1);
% errorbar((num_syn(1:idx_7mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);
plot(num_syn(1:idx_7mm), rsaf_tmpM, 'LineWidth', 3);
max_val = max(rsaf_tmpM);

plot((num_syn), ref_, 'LineWidth', 4, 'color', 'b'); hold off;
% ylim([0 1.2*max_val]);

legend('h = 3mm', 'h = 4mm', 'h = 5mm', 'h = 6mm', 'h = 7mm', 'ref','location','southeast');

% 7cm
depth_ = 7;
n_3mm = 55; idx_3mm = find(num_syn == n_3mm);
n_4mm = 65; idx_4mm = find(num_syn == n_4mm);
n_5mm = 65; idx_5mm = find(num_syn == n_5mm);
n_6mm = 65; idx_6mm = find(num_syn == n_6mm);
n_7mm = 65; idx_7mm = find(num_syn == n_7mm);

figure(depth_); hold on;

rsaf_tmpM = mean(squeeze(h_3mm(:,depth_,2:idx_3mm+1)))./mean(squeeze(h_3mm(:,depth_,1)));
rsaf_tmpStd = std(squeeze(h_3mm(:,depth_,2:idx_3mm+1)./h_3mm(:,depth_,1)),0,1);
% errorbar((num_syn(1:idx_3mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);
plot(num_syn(1:idx_3mm), rsaf_tmpM, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(h_4mm(:,depth_,2:idx_4mm+1)))./mean(squeeze(h_4mm(:,depth_,1)));
rsaf_tmpStd = std(squeeze(h_4mm(:,depth_,2:idx_4mm+1)./h_4mm(1,depth_,1)),0,1);
% errorbar((num_syn(1:idx_4mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);
plot(num_syn(1:idx_4mm), rsaf_tmpM, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(h_5mm(:,depth_,2:idx_5mm+1)))./mean(squeeze(h_5mm(:,depth_,1)));
rsaf_tmpStd = std(squeeze(h_5mm(:,depth_,2:idx_5mm+1)./h_5mm(1,depth_,1)),0,1);
% errorbar((num_syn(1:idx_5mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);
plot(num_syn(1:idx_5mm), rsaf_tmpM, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(h_6mm(:,depth_,2:idx_6mm+1)))./mean(squeeze(h_6mm(:,depth_,1)));
rsaf_tmpStd = std(squeeze(h_6mm(:,depth_,2:idx_6mm+1)./h_6mm(1,depth_,1)),0,1);
% errorbar((num_syn(1:idx_6mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);
plot(num_syn(1:idx_6mm), rsaf_tmpM, 'LineWidth', 3);

rsaf_tmpM = mean(squeeze(h_7mm(:,depth_,2:idx_7mm+1)))./mean(squeeze(h_7mm(:,depth_,1)));
rsaf_tmpStd = std(squeeze(h_7mm(:,depth_,2:idx_7mm+1)./h_7mm(1,depth_,1)),0,1);
% errorbar((num_syn(1:idx_7mm)),rsaf_tmpM, rsaf_tmpStd, 'LineWidth', 3);
plot(num_syn(1:idx_7mm), rsaf_tmpM, 'LineWidth', 3);
max_val = max(rsaf_tmpM);

plot((num_syn), ref_, 'LineWidth', 4, 'color', 'b'); hold off;
% ylim([0 1.2*max_val]);
legend('h = 3mm', 'h = 4mm', 'h = 5mm', 'h = 6mm', 'h = 7mm', 'ref','location','southeast');







