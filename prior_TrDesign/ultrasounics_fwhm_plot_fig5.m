clc; clear; close all;
%%
max_syn = [53 49 45];
R_ = [5 10 15];

d_theta = 0.4724; % deg
d_rad = d_theta*pi/180;
aper_size = (max_syn*d_rad).*R_;

figure;hold on;
yyaxis right
bar(R_, aper_size);
yyaxis left
plot(R_, max_syn); 
hold off;
%%
num_syn = [3 5 7 9 11 13 15 17 19 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65];
depth_ = [1 2 3 4 5 6 7]*10;

%% 5mm
load('/Users/songhyunwoo/Documents/JohnsHopkins/Research/ongoing/rSAF/src/data/rsaf_optim_critical_param/rsaf_optim_R/r05mm_fele5_h5/FWHM/FWHM_rSAF.mat');
syn_ = max_syn(1);
syn_idx = find(num_syn == syn_);
fwhm_5mm = FWHM_rSAF(syn_idx,:);

%% 10mm
load('/Users/songhyunwoo/Documents/JohnsHopkins/Research/ongoing/rSAF/src/data/rsaf_optim_critical_param/rsaf_optim_R/r10mm_fele5_h5/FWHM/FWHM_rSAF.mat');
syn_ = max_syn(2);
syn_idx = find(num_syn == syn_);
fwhm_10mm = FWHM_rSAF(syn_idx,:);

%% 15mm
load('/Users/songhyunwoo/Documents/JohnsHopkins/Research/ongoing/rSAF/src/data/rsaf_optim_critical_param/rsaf_optim_R/r15mm_fele5_h5/FWHM/FWHM_rSAF.mat');
syn_ = max_syn(3);
syn_idx = find(num_syn == syn_);
fwhm_15mm = FWHM_rSAF(syn_idx,:);

%%
figure;
plot(depth_, fwhm_5mm); hold on;
plot(depth_, fwhm_10mm);
plot(depth_, fwhm_15mm); hold off;
legend('R = 5mm', 'R = 10mm', 'R = 15mm', 'location', 'northwest');
