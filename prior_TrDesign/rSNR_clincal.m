clc; clear; close all;

%%
dir_ = uigetdir('./data','');

folder_list = dir(dir_);
flag = 0;
if(strcmp(folder_list(3).name,'.DS_Store')), flag = 1; end
folder_list = folder_list(3+flag:end);

%% clinical transducer (i.e., f_ele = 25mm)
clinical_case = folder_list(1).name;
dir_clinical = [dir_ '/' clinical_case '/SNR'];
% load SNR
load([dir_clinical '/snr_var.mat']);

snr_clinic_CON = mean(snr_(:,:,1));
snr_clinic_rSAF = squeeze(mean(snr_(:,:,2:end)));

%% optimized transducer (i.e., f_ele = 5mm, 10mm, 15mm, 20mm, 25mm)
snr_optim_rSAF = {};
for case_idx = 2:numel(folder_list)
    optim_case = folder_list(case_idx).name;
    disp(optim_case);
    dir_optim = [dir_ '/' optim_case '/SNR'];
    % load SNR
    load([dir_optim '/snr_var.mat']);
    
    snr_tmp = squeeze(mean(snr_(:,:,2:end)));
    snr_optim_rSAF = cat(3, snr_optim_rSAF, snr_tmp);
end

%% num_syn
num_syn_c = {};

num_syn_tmp = [3 5 7 9 11 13 15 17 19 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65];
num_syn_c = cat(3, num_syn_c, num_syn_tmp);

num_syn_tmp = [3 5 7 9 11 13 15 17 19 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65 67 69 71 73 ...
    75 77 79 81 83];
num_syn_c = cat(3, num_syn_c, num_syn_tmp);

num_syn_tmp = [3 5 7 9 11 13 15 17 19 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65 67 69 71 73 ...
    75 77 79 81 83 85 87 89 91 93 95 97 99];
num_syn_c = cat(3, num_syn_c, num_syn_tmp);

num_syn_tmp = [3 5 7 9 11 13 15 17 19 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65 67 69 71 73 ...
    75 77 79 81 83 85 87 89 91 93 95 97 99 101 103 105 107 109 111 113 115];
num_syn_c = cat(3, num_syn_c, num_syn_tmp);

num_syn_tmp = [3 5 7 9 11 13 15 17 19 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65 67 69 71 73 ...
    75 77 79 81 83 85 87 89 91 93 95 97 99 101 103 105 107 109 111 113 115 117 119 121 123 125 127 129 131];
num_syn_c = cat(3, num_syn_c, num_syn_tmp);
%% plot rSNR
depth_ = [1 3 5 7];


for d_idx = 1:numel(depth_)
    
    clinic_tmp = snr_clinic_CON(depth_(d_idx));
    
    leg_ = {};
    figure(depth_(d_idx));
    for case_idx = 1:numel(snr_optim_rSAF)
        num_syn = num_syn_c{case_idx};
        
        optim_tmp = snr_optim_rSAF{case_idx};
        optim_tmp = optim_tmp(depth_(d_idx),:);
        
        rSNR_tmp = optim_tmp / clinic_tmp;
        leg_ = cat(3,leg_,folder_list(case_idx+1).name);
        plot(num_syn,rSNR_tmp, 'LineWidth', 2, 'color', [0.12 0.12 0.12]*case_idx); hold on;
    end
    plot(num_syn, ones(1,numel(num_syn)), 'LineWidth', 2, 'color', 'k','LineStyle','--'); hold off;
    ylim([0 3]); xlim([0 num_syn(end)]);
    legend(leg_,'Interpreter','none','location','southeast');
end
%% single depth
depth_ = 1;
clinic_tmp = snr_clinic_CON(depth_);

leg_ = {};
figure(depth_);
for case_idx = 1:numel(snr_optim_rSAF)
    num_syn = num_syn_c{case_idx};
    
    optim_tmp = snr_optim_rSAF{case_idx};
    optim_tmp = optim_tmp(depth_,:);
    
    rSNR_tmp = (optim_tmp / clinic_tmp);
    leg_ = cat(3,leg_,folder_list(case_idx+1).name);
    plot(num_syn,rSNR_tmp, 'LineWidth', 2, 'color', [0.12 0.12 0.12]*case_idx); hold on;
end
plot(num_syn, ones(1,numel(num_syn)), 'LineWidth', 2, 'color', 'k','LineStyle','--'); hold off;
ylim([0 3]); xlim([0 num_syn(end)]);
legend(leg_,'Interpreter','none','location','southeast');