%%% Performance degradation plot
%%% error vs e_max
%%% worst case
clc; clear all; close all;
%% Parameter setting
bCutoff = 0;
cut_off_value = 20;

%%%%%%%
idx = 2;       % 1: conv, 2: SA
depth = 6;     % interested depth;
size_Z = 10e-3;
size_Y = 35e-3;
%%%%%%%

%% Select master path
master_path = uigetdir('', '');
path_ground_truth = [master_path '/Ground truth/Saved data'];
path_error = [master_path '/Errors'];

%% Folder load (ground truth)
disp('load ground truth data');
gt_beamforming_list = dir(path_ground_truth);
load([path_ground_truth '/' gt_beamforming_list(idx+2).name]);
img_ground_truth = stSaveInfo.mImg_DSC;
y_axis = stSaveInfo.aYAxis_DSC;
z_axis = stSaveInfo.aZAxis_DSC;
radius = stSaveInfo.nRadius;
clear stSaveInfo;

yIdx = (y_axis >= -0.5*size_Y) & (y_axis <= 0.5*size_Y);
zIdx = (z_axis-radius >= depth*1e-2 - 0.5*size_Z) & (z_axis-radius <= depth*1e-2 + 0.5*size_Z);
gt_roi = img_ground_truth(zIdx, yIdx);

%% Folder load (error)
disp('load error data');
error_list = dir(path_error);
no_of_error_case = numel(error_list)-2;
no_of_batch = numel(dir([path_error '/' error_list(3).name]))-2;

max_xcorr_map = zeros(no_of_batch, no_of_error_case);
for eIdx = 3:numel(error_list)
    case_name = error_list(eIdx).name;
    path_error_case = [path_error '/' case_name];
    batch_list = dir(path_error_case);
    for bIdx = 3:numel(batch_list)
        sample_name = [batch_list(bIdx).name '/Saved data'];
        error_beamforming_method = dir([path_error_case '/' sample_name]);
        load([path_error_case '/' sample_name '/' error_beamforming_method(idx+2).name]);
        disp(['Case: ', error_list(eIdx).name '<<<>>> Sample: ', num2str(bIdx-2) ', ' error_beamforming_method(idx+2).name]);
        img_error = stSaveInfo.mImg_DSC;
        error_roi = img_error(zIdx, yIdx);
        % cross correlation
        correlation_map = normxcorr2(gt_roi, error_roi);
        max_xcorr_map(bIdx-2, eIdx-2) = max(correlation_map(:));
    end
end

%% Plot
figure(1); boxplot(max_xcorr_map, 'Labels', {'0.1', '0.2', '0.3', '0.4', '0.5'});
xlabel('E_max [deg]', 'Interpreter', 'None'); ylabel('XCORR score'); grid on; ylim([0.96 1]);
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf,'Position',[379 53 1500 676]);