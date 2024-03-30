clc; clear; close all;
%%
max_syn = [23 37 51 61 65];
h_ = [3 4 5 6 7];

figure;
plot(h_, max_syn);
%%
num_syn = [3 5 7 9 11 13 15 17 19 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65];
depth_ = [1 2 3 4 5 6 7]*10;

%% 3mm
load('/Users/songhyunwoo/Documents/JohnsHopkins/Research/ongoing/rSAF/src/data/rsaf_optim_critical_param/rsaf_optim_h/r10mm_h3mm/FWHM/FWHM_rSAF.mat')
syn_ = max_syn(1);
syn_idx = find(num_syn == syn_);
fwhm_3mm = FWHM_rSAF(syn_idx,:);

%% 4mm
load('/Users/songhyunwoo/Documents/JohnsHopkins/Research/ongoing/rSAF/src/data/rsaf_optim_critical_param/rsaf_optim_h/r10mm_h4mm/FWHM/FWHM_rSAF.mat')
syn_ = max_syn(2);
syn_idx = find(num_syn == syn_);
fwhm_4mm = FWHM_rSAF(syn_idx,:);

%% 5mm
load('/Users/songhyunwoo/Documents/JohnsHopkins/Research/ongoing/rSAF/src/data/rsaf_optim_critical_param/rsaf_optim_h/r10mm_h5mm/FWHM/FWHM_rSAF.mat')
syn_ = max_syn(3);
syn_idx = find(num_syn == syn_);
fwhm_5mm = FWHM_rSAF(syn_idx,:);

%% 6mm
load('/Users/songhyunwoo/Documents/JohnsHopkins/Research/ongoing/rSAF/src/data/rsaf_optim_critical_param/rsaf_optim_h/r10mm_h6mm/FWHM/FWHM_rSAF.mat')
syn_ = max_syn(4);
syn_idx = find(num_syn == syn_);
fwhm_6mm = FWHM_rSAF(syn_idx,:);

%% 7mm
load('/Users/songhyunwoo/Documents/JohnsHopkins/Research/ongoing/rSAF/src/data/rsaf_optim_critical_param/rsaf_optim_h/r10mm_h7mm/FWHM/FWHM_rSAF.mat')
syn_ = max_syn(5);
syn_idx = find(num_syn == syn_);
fwhm_7mm = FWHM_rSAF(syn_idx,:);

%%
figure;
plot(depth_, fwhm_3mm); hold on;
plot(depth_, fwhm_4mm);
plot(depth_, fwhm_5mm);
plot(depth_, fwhm_6mm);
plot(depth_, fwhm_7mm); hold off;
legend('h = 3mm', 'h = 4mm', 'h = 5mm', 'h = 6mm', 'h = 7mm', 'location', 'northwest');

%%
