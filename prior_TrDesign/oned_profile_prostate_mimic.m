clc; close all; clear;

addpath('functions');

%%
dir_ = uigetdir('./data','');

folder_list = dir(dir_);
flag = 0;
if(strcmp(folder_list(3).name,'.DS_Store')), flag = 1; end
folder_list = folder_list(3+flag:end);

%% load data
imgs_ = {};
y_ = {}; z_ = {};
radius_ = []; fov_ = [];
legend_ = {};

for f_idx = 1:numel(folder_list)
    folder_tmp = folder_list(f_idx).name;
    disp(['>>> ' folder_tmp ]);
    dir_tmp = [dir_ '/' folder_tmp];
    
    legend_ = cat(1, legend_, folder_tmp);
    %
    load([dir_tmp '/Parameters.mat']);
    acoustic_ = stParam.stRFInfo;
    bf_ = stParam.stBFInfo;
    trans_ = stParam.stTRInfo;
    
    imgpos_y = stParam.mImgY;
    imgpos_z = stParam.mImgZ;
    
    clear stParam;
    
    % mid processing paramete
    mid_.nTGC_Atten = 0.5;                % [dB]
    
    mid_.nDCRType = 'high';
    mid_.nDCRTap = 128;                   % BPF tap #
    mid_.nDCRFcut = 1e6;
    
    % dsc parameter
    scanline_theta = linspace(-0.5*bf_.nFOV, 0.5*bf_.nFOV, bf_.nScline); % Ground truth transmitted angle
    depth_ = linspace(bf_.nRadius, bf_.nRadius+bf_.nDth, bf_.nDthSpl);
    
    da = abs(scanline_theta(1)-scanline_theta(2));
    dr = abs(depth_(1)-depth_(2));
    view_depth = bf_.nDth + (bf_.nRadius*(1-cosd(0.5*bf_.nFOV)));
    view_width = 2* (bf_.nRadius+bf_.nDth)*sind(0.5*bf_.nFOV);
    
    dz = 1e-4;
    dy = 1e-4;
    height = round(view_depth / dz);
    width = round(view_width / dy);
    
    % load data
    dir_data = [dir_tmp '/errors_bf/error_0/Sample001/Element_64/'];
    file_list = dir(dir_data);
    
    load([dir_data '/' file_list(end).name]);
    
    beamformed_data = stSaveInfo.mBFedData;
    beamformed_data = InterpNan(beamformed_data);
    env_data = mid_proc(beamformed_data, mid_, acoustic_, bf_);
    
    [axis_y, axis_z, dsc_data] = ScanConverter_convex(env_data, dr, da, bf_.nRadius, height, width, dz, dy);
    dsc_data(dsc_data == 50) = 0;
    imgs_ = cat(3, imgs_, dsc_data);
    y_ = cat(3, y_, axis_y);
    z_ = cat(3, z_, axis_z);
    radius_ = cat(2, radius_, bf_.nRadius);
    fov_ = cat(2, fov_, bf_.nFOV);
end

%% get roi
depth_pos = [12.5 27.5 42.5 57.5]*1e-3; % point target

% depth_pos = [20 35 50]*1e-3;

dor = 1e-3; % point target
% dor = 6e-3;

roi_ = {};
for f_idx = 1:numel(folder_list)
    roi_d = [];
    img_tmp = imgs_{f_idx};
    y_tmp = y_{f_idx};
    z_tmp = z_{f_idx};
    radius_tmp = radius_(f_idx);
    reject_idx = find(img_tmp == 0);
   
    for d_idx = 1:numel(depth_pos)
        depth_tmp = depth_pos(d_idx) + radius_tmp;
        
        mask_map = getROIMask(y_tmp, z_tmp, depth_tmp, radius_tmp, dor, reject_idx, 0);
        
        roi_tmp = zeros(size(mask_map));
        roi_tmp(mask_map) = img_tmp(mask_map);
        
        % mean projection
        mean_prj = mean(roi_tmp,1);
        
        roi_d = cat(1, roi_d, mean_prj);
        figure(10201);
        subplot(numel(depth_pos), numel(folder_list), (d_idx-1)*numel(folder_list) + f_idx);
        imagesc(y_tmp, z_tmp-radius_tmp,db(roi_tmp/max(img_tmp(:))));
        axis equal; axis tight; caxis([-60 0]); colormap gray;
    end
    roi_ = cat(3, roi_, roi_d);
end
set(gcf,'Position',[75 202 1541 745]);
%% plot
gt_posY = [-7.6 -2.74 1.25 4.21 6.18 7.13]*1e-3;
gt_posZ = [12.4 12.3 12.5 12.1 12.4 12.3;
           27.0 27.4 27.5 27.3 27.0 26.8;
           41.9 42.4 42.5 42.3 42.1 42.0;
           57.1 57.4 57.5 57.4 57.2 57.1]*1e-3;

close all;
figure(1010); imagesc(axis_y, axis_z-15e-3, db(dsc_data/max(dsc_data(:)))); colormap gray; caxis([-60 0]);
for f_idx = 1:numel(folder_list)
    roi_tmp = roi_{f_idx};
    y_tmp = y_{f_idx};
    for d_idx = 1:numel(depth_pos)
        psf_tmp = roi_tmp(d_idx,:);
        idx_tmp = find(psf_tmp,1,'first'):find(psf_tmp,1,'last');
        psf_tmp = psf_tmp(idx_tmp);
        
%         axis_ = atand(y_tmp(idx_tmp)/depth_pos(d_idx));
        axis_ = y_tmp(idx_tmp);
        
        figure(d_idx);
        plot(axis_*1e3, (psf_tmp/max(psf_tmp)), 'LineWidth',2, 'color', [0.12 0.12 0.12]*f_idx); hold on;
    end
end

for d_idx = 1:numel(depth_pos)
    gt_z_tmp = gt_posZ(d_idx,:);
%     gt_ang = atand(gt_posY./gt_z_tmp);
%     gt_ang = gt_posY;
    
    figure(d_idx);
%     for l_idx = 1:numel(gt_ang)
%         line([gt_ang(l_idx) gt_ang(l_idx)]*1e3, [0 1],'color','r');
%     end
    
%     xlim([-35 35]);
    xlim([-15 15]); % point target
%     xlim([-30 30]);
    ylim([0 1]);
    set(gcf,'Position',[384 592 886 314]);
    legend(legend_,'Interpreter','none','location','northeast');
end








