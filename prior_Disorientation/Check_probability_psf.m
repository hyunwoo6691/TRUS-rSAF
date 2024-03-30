clc; clear; close all;
addpath('02.Functions');

%%%%%%%%%%%%%%%%%%%%
dB = -6;
%%%%%%%%%%%%%%%%%%%%

% axis setting
a_depth = 1:10; % cm
% a_errors = 0.1*(1:10);
a_errors = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 2.0 5.0];
% a_focal_points = 8:1:25; % mm
a_focal_points = 25;

%% get master directory
dir_master = uigetdir('01.Data','');
dir_figure = [dir_master '/Figures'];
dir_master = [dir_master '/Data'];
focal_point_list = dir(dir_master);
flag = 0;
if(strcmp(focal_point_list(3).name, '.DS_Store')), flag = 1; end
focal_point_list = focal_point_list(3+flag:end);

%% Get the DSC image data

for f_idx = 1:numel(focal_point_list)
    focal_point_tmp = focal_point_list(f_idx).name;
    disp(['>>> Folder: ' focal_point_tmp]);
    %% get directory & data loading
    dir_ = [dir_master '/' focal_point_tmp];
    
    dir_resol = [dir_ '/Resolutions_' num2str(dB) 'dB'];
    dir_err = [dir_ '/Errors'];
    
    err_list = dir(dir_err);
    flag = 0;
    if(strcmp(err_list(3).name, '.DS_Store')), flag = 1; end
    err_list = err_list(4+flag:end);
    
    total_conv = {};
    total_sa = {};
    for e_idx = 1:numel(err_list)
        err_case = err_list(e_idx).name;
        
        disp(err_case);
        
        dir_tmp = [dir_err '/' err_case];
        sample_list = dir(dir_tmp);
        flag = 0;
        if(strcmp(sample_list(3).name, '.DS_Store')), flag = 1; end
        sample_list = sample_list(3+flag:end);
        
        v_img_conv = [];
        v_img_sa = [];
        for s_idx = 1:numel(sample_list)
            sample_tmp = sample_list(s_idx).name;
            tmp_list = dir([dir_tmp '/' sample_tmp]);
            tmp_conv = tmp_list(end-1).name;
            tmp_sa = tmp_list(end).name;
            
            load([dir_tmp '/' sample_tmp '/' tmp_conv]);
            img_tmp = stSaveInfo.mImg_DSC;
            v_img_conv = cat(3, v_img_conv, img_tmp);
            load([dir_tmp '/' sample_tmp '/' tmp_sa]);
            img_tmp = stSaveInfo.mImg_DSC;
            v_img_sa = cat(3, v_img_sa, img_tmp);
            clear stSaveInfo;
        end
        total_conv = cat(3, total_conv, v_img_conv);
        total_sa = cat(3, total_sa, v_img_conv);
    end
    
   
end

%% take mean error SA < Conv
mask_mean = v_err_sa < v_err_conv;

statistics_mean = [];
for d_idx = 1:numel(a_depth)
    h_tmp = squeeze(h_(d_idx, :, :));
    mask_tmp = squeeze(mask_mean(d_idx, :, :));
    
    result_tmp = h_tmp .* mask_tmp;
    statistics_mean = cat(3, statistics_mean, result_tmp);
end

statistics_sum_mean = sum(statistics_mean, 3);

figure(1);
imagesc(a_focal_points, a_errors, statistics_sum_mean); colorbar; caxis([0 numel(a_depth)]);
ylabel('E_max [deg]', 'Interpreter', 'None'); xlabel('Focal point [mm]'); title('Mean error');
set(gcf, 'Position', [169 167 905 296]);

%% take std SA < Conv
mask_std = v_std_sa < v_std_conv;

statistics_std = [];
for d_idx = 1:numel(a_depth)
    h_tmp = squeeze(h_(d_idx, :, :));
    mask_tmp = squeeze(mask_std(d_idx, :, :));
    
    result_tmp = h_tmp .* mask_tmp;
    statistics_std = cat(3, statistics_std, result_tmp);
end

statistics_sum_std = sum(statistics_std, 3);

figure(2);
imagesc(a_focal_points, a_errors, statistics_sum_std); colorbar; caxis([0 numel(a_depth)]);
ylabel('E_max [deg]', 'Interpreter', 'None'); xlabel('Focal point [mm]'); title('STD');
set(gcf, 'Position', [632 77 905 296]);

%% case when sa is better than conventional (ground truth image)
mask_enhancement = enhancement_ < 1;
mask_vs_behind = zeros(size(enhancement_));
for d_idx = 1: numel(a_depth)
    for f_idx = 1:numel(focal_point_list)
        if (a_depth(d_idx)*10 > a_focal_points(f_idx))
            mask_vs_behind(d_idx, f_idx) = 1;
        end
    end
end
mask_p = mask_enhancement.*mask_vs_behind;

figure(3);
sgtitle('Resolution enhancement [ Resol_SA / Resol_Conv ] ', 'Interpreter', 'None');
subplot(1,2,1); imagesc(a_focal_points, a_depth*10, enhancement_); colorbar;
ylabel('Depth [mm]'); xlabel('Focal point [mm]');
subplot(1,2,2); imagesc(a_focal_points, a_depth*10, mask_p); colorbar;
ylabel('Depth [mm]'); xlabel('Focal point [mm]');
set(gcf, 'Position', [172 599 905 296]);

%% spatial resolution of ground truth
figure(4)
sgtitle('Spatial resolution of ground truth');
subplot(1,2,1); imagesc(a_focal_points, a_depth*10, gt_int_conv); colorbar;
ylabel('Depth [mm]'); xlabel('Focal point [mm]'); title('Conventional'); caxis([0 max(gt_int_conv(:))]);
subplot(1,2,2); imagesc(a_focal_points, a_depth*10, gt_int_sa); colorbar;
ylabel('Depth [mm]'); xlabel('Focal point [mm]'); title('Synthetic Aperture'); caxis([0 max(gt_int_conv(:))]);

set(gcf, 'Position', [682 481 905 296]);

%%  Get the optimized focusing point
best_mean_case= (statistics_sum_mean == max(statistics_sum_mean(:)));
best_std_case = (statistics_sum_std == max(statistics_sum_std(:)));
% best_mean_case= (statistics_sum_mean >= 8);
% best_std_case = (statistics_sum_std >= 8);

best_performance = best_mean_case.*best_std_case;

figure(5);
imagesc(a_focal_points, a_errors, best_performance); colorbar; caxis([0 1]);
ylabel('E_max [deg]', 'Interpreter', 'None'); xlabel('Focal point [mm]'); title('Best performance at');
set(gcf, 'Position', [632 77 905 296]);

%% Plot best performance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
best_focal_point = 25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

focal_idx = find(a_focal_points == best_focal_point);

best_focal_point_ = focal_point_list(focal_idx).name;
% get directory & data loading
dir_best = [dir_master '/' best_focal_point_];

dir_resol_best = [dir_best '/Resolutions_' num2str(dB) 'dB'];
dir_err_best = [dir_best '/Errors'];

load([dir_resol_best '/error_0/Sample001.mat']);
resol_gt_best = mLateralResol;
clear mLateralResol;

% ground truth FWHM along depth
figure(6);
plot(a_depth*10, resol_gt_best(1,:), 'LineWidth', 3, 'LineStyle', '--', 'Marker', '+', 'MarkerSize', 10); hold on;
plot(a_depth*10, resol_gt_best(2,:), 'LineWidth', 3, 'LineStyle', '--', 'Marker', 'o', 'MarkerSize', 10); hold off;
grid on; ylim([0 25]);
xlabel('depth [mm]'); ylabel('FWHM [mm]'); title('Ground truth FWHM');
legend('Conventional', 'SA', 'location', 'Northwest');
set(gcf, 'Position', [509 498 990 355]);


% std of (ground truth - error)
std_best_conv = v_std_conv(:,:,focal_idx);
std_best_sa = v_std_sa(:,:,focal_idx);

figure(7); sgtitle('Standard deviation');
subplot(1,2,1); imagesc(a_errors, a_depth*10, std_best_conv); colorbar;
ylabel('Depth [mm]'); xlabel('E_max [deg]', 'Interpreter', 'None'); title('Conventional'); 
caxis([0 max(max(std_best_conv(:)),max(std_best_sa(:)))]);
subplot(1,2,2); imagesc(a_errors, a_depth*10, std_best_sa); colorbar;
ylabel('Depth [mm]'); xlabel('E_max [deg]', 'Interpreter', 'None'); title('Synthetic Aperture'); 
caxis([0 max(max(std_best_conv(:)),max(std_best_sa(:)))]);

set(gcf, 'Position', [682 481 905 296]);

figure(8);
imagesc(a_errors, a_depth*10, std_best_conv./std_best_sa); colorbar;
ylabel('Depth [mm]'); xlabel('E_max [deg]', 'Interpreter', 'None'); caxis([1 2]);

%%
err_axis = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 5];
err_axiss = [0.1 0.2 0.5 1 2 5];
% err_axiss = err_axis;

mask_ = zeros(1, numel(err_axis));
for k = 1:numel(err_axiss)
    idx_tmp = find(err_axiss(k) == err_axis);
    mask_(idx_tmp) = 1;
end
mask_ = logical(mask_);

mean_of_error_conv = mean(err_conv_all);
mean_of_error_sa = mean(err_sa_all);

std_of_error_conv = std(err_conv_all);
std_of_error_sa = std(err_sa_all);

% plot
figure(9);
plot(err_axiss, mean_of_error_conv(mask_), 'LineWidth', 2.5, 'Marker', '*', 'MarkerSize', 10, 'LineStyle', '--'); hold on;
plot(err_axiss, mean_of_error_sa(mask_), 'LineWidth', 2.5, 'Marker', '*', 'MarkerSize', 10, 'LineStyle', '--'); hold off;
grid on; legend('CON', 'eSAF', 'location', 'northeast');
xlabel('\epsilon _{max}'); ylabel('Distortion parameter [mm]');
set(gca,'Fontsize',14,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf, 'Position', [-240 2049 697 266]);

figure(10);
plot(err_axiss, (std_of_error_conv(mask_)), 'LineWidth', 2.5, 'Marker', '*', 'MarkerSize', 10, 'LineStyle', '--'); hold on;
plot(err_axiss, (std_of_error_sa(mask_)), 'LineWidth', 2.5, 'Marker', '*', 'MarkerSize', 10, 'LineStyle', '--'); hold off;
grid on; legend('CON', 'eSAF', 'location', 'northeast');
xlabel('\epsilon _{max}'); ylabel('Robustness score');
set(gca,'Fontsize',14,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf, 'Position', [-238 1689 697 266]);

