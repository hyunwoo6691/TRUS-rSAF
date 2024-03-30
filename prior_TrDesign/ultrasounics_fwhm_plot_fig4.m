clc; clear; close all;
%%
max_syn = [49 23 11 7 5];
fele = [5 10 15 20 25];

figure;
plot(fele, max_syn);

%%
num_syn = [3 5 7 9 11 13 15 17 19 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65];
depth_ = [1 2 3 4 5 6 7]*10;

%% 5mm
load('/Users/songhyunwoo/Documents/JohnsHopkins/Research/ongoing/rSAF/src/data/rsaf_optim_critical_param/rsaf_optim_fele/r10mm_lens_focus_05mm/FWHM/FWHM_rSAF.mat')
syn_ = 49;
syn_idx = find(num_syn == syn_);
fwhm_5mm = FWHM_rSAF(syn_idx,:);

%% 10mm
load('/Users/songhyunwoo/Documents/JohnsHopkins/Research/ongoing/rSAF/src/data/rsaf_optim_critical_param/rsaf_optim_fele/r10mm_lens_focus_10mm/FWHM/FWHM_rSAF.mat')
syn_ = 23;
syn_idx = find(num_syn == syn_);
fwhm_10mm = FWHM_rSAF(syn_idx,:);

%% 15mm
load('/Users/songhyunwoo/Documents/JohnsHopkins/Research/ongoing/rSAF/src/data/rsaf_optim_critical_param/rsaf_optim_fele/r10mm_lens_focus_15mm/FWHM/FWHM_rSAF.mat')
syn_ = 11;
syn_idx = find(num_syn == syn_);
fwhm_15mm = FWHM_rSAF(syn_idx,:);

%% 20mm
load('/Users/songhyunwoo/Documents/JohnsHopkins/Research/ongoing/rSAF/src/data/rsaf_optim_critical_param/rsaf_optim_fele/r10mm_lens_focus_20mm/FWHM/FWHM_rSAF.mat')
syn_ = 7;
syn_idx = find(num_syn == syn_);
fwhm_20mm = FWHM_rSAF(syn_idx,:);

%% 25mm
load('/Users/songhyunwoo/Documents/JohnsHopkins/Research/ongoing/rSAF/src/data/rsaf_optim_critical_param/rsaf_optim_fele/r10mm_lens_focus_25mm/FWHM/FWHM_rSAF.mat')
syn_ = 19;
syn_idx = find(num_syn == syn_);
fwhm_25mm = FWHM_rSAF(syn_idx,:);

%%
figure;
plot(depth_, fwhm_5mm); hold on;
plot(depth_, fwhm_10mm);
plot(depth_, fwhm_15mm);
plot(depth_, fwhm_20mm);
plot(depth_, fwhm_25mm); hold off;
legend('f_{ele} = 5mm', 'f_{ele} = 10mm', 'f_{ele} = 15mm', 'f_{ele} = 20mm', 'f_{ele} = 25mm', 'location', 'northwest');

%%
