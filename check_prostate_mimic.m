clc; clear; close all;
addpath('functions');
%% set directory
dir_ = uigetdir('./data','');

%% load parameter
load([dir_ '/Parameters.mat']);
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

%% set data directory
dir_data = [dir_ '/errors_bf'];

total_case = dir(dir_data);
flag = 0;
if(strcmp(total_case(3).name,'.DS_Store')), flag = 1; end
total_case = total_case(3+flag:end);

dir_gt = [dir_data '/' total_case(1).name];

err_case = total_case(2:end);

%% GROUND TRUTH : show image
dir_gt_ = [dir_gt '/Sample001/Element_64'];
bf_cases = dir(dir_gt_);

disp('>>> load ground truth data');
load([dir_gt_ '/' bf_cases(end).name]);
env_data = mid_proc(stSaveInfo.mBFedData, mid_, acoustic_, bf_);
[axis_y, axis_z, dsc_gt, reject_idx] = dsc(env_data, dr, da, bf_, height, width, dz, dy);

figure(1);
imagesc(axis_y*1e3, (axis_z-bf_.nRadius)*1e3, db(dsc_gt/max(dsc_gt(:))));
axis equal; axis tight; colormap gray; title('GROUND TRUTH');
caxis([-50 0]);
set(gcf, 'Position', [561 431 936 516]);
disp('done'); pause(1);

%% ERROR : load & dsc
disp('>>> load error data');
dsc_err = {}; dsc_err_avg = {};

for e_idx = 1:numel(err_case)
    disp(['    ' err_case(e_idx).name]);
    dir_tmp = [dir_data '/' err_case(e_idx).name];
    sample_list = dir(dir_tmp);
    
    flag = 0;
    if(strcmp(sample_list(3).name,'.DS_Store')), flag = 1; end
    sample_list = sample_list(3+flag:end);
    
    dsc_tmp = [];
    for s_idx = 1:numel(sample_list)
        if(mod(s_idx, numel(sample_list)/4) == 0), disp(['      ' num2str(round(100*s_idx/numel(sample_list))) '%...']); end
        dir_tmpp = [dir_tmp '/' sample_list(s_idx).name '/Element_64'];
        
        load([dir_tmpp '/' bf_cases(end).name]);
%         env_data = mid_proc(stSaveInfo.mBFedData, mid_, acoustic_, bf_);
%         [axis_y, axis_z, mOutput, reject_idx] = dsc(env_data, dr, da, bf_, height, width, dz, dy);
        
%         dsc_tmp = cat(3, dsc_tmp, mOutput);
        
        % average filter
        [beamformed_data_fil, fil] = spatial_filtering(stSaveInfo.mBFedData, 0.1, bf_.nDTheta, 'boxcar');
        env_data = mid_proc(beamformed_data_fil, mid_, acoustic_, bf_);
        [axis_y, axis_z, mOutput, reject_idx] = dsc(env_data, dr, da, bf_, height, width, dz, dy);
        dsc_tmp = cat(3, dsc_tmp, mOutput);
    end
    
    dsc_err = cat(3, dsc_err, dsc_tmp);
end
disp('done');

%% BOTH : close up 2d-psf
dim_ = [70 70]*1e-3; % horizontal, vertical
center_ = [0, 35]*1e-3; % horizontal, vertial

% get axis
range_hor = [center_(1)-0.5*dim_(1) center_(1)+0.5*dim_(1)];
range_ver = [center_(2)-0.5*dim_(2) center_(2)+0.5*dim_(2)];

idx_y = [find(abs(range_hor(1)-axis_y) == min(abs(range_hor(1)-axis_y))) find(abs(range_hor(2)-axis_y) == min(abs(range_hor(2)-axis_y)))];
idx_z = [find(abs(range_ver(1)-axis_z+bf_.nRadius) == min(abs(range_ver(1)-axis_z+bf_.nRadius))) find(abs(range_ver(2)-axis_z+bf_.nRadius) == min(abs(range_ver(2)-axis_z+bf_.nRadius)))];

axis_yy = axis_y(idx_y(1):idx_y(2));
axis_zz = axis_z(idx_z(1):idx_z(2)) - bf_.nRadius;

% ground truth
roi_gt = dsc_gt(idx_z(1):idx_z(2),idx_y(1):idx_y(2),1);

figure(1231); imagesc(axis_yy,axis_zz,db(roi_gt/max(max(dsc_gt(:,:,1))))); colormap gray; caxis([-50 0]); axis tight; axis equal; set(gcf,'Position',[360 420 370 370]);

% error
indices_ = [1 1 1 1 1 1];
for e_idx = 1:numel(err_case)
    dsc_err_tmp = dsc_err{e_idx};
    
    roi_tmp = dsc_err_tmp(idx_z(1):idx_z(2),idx_y(1):idx_y(2),indices_(e_idx));
    
    figure(12314+e_idx); % CON
    imagesc(axis_yy,axis_zz,db(roi_tmp/max(dsc_err_tmp(:)))); colormap gray; caxis([-50 0]); axis tight; axis equal; set(gcf,'Position',[360+400*e_idx 420 370 370]);
end

%% SET ROI POSITION
depth_pos = [12.5 27.5 42.5 57.5]*1e-3; % point target

% depth_pos = [20 35 50]*1e-3;

dor = 1e-3; % point target
% dor = 6e-3; % mass target
%% GROUND TRUTH : extract 1d psf (mean projection)
disp('>>> extract ground truth data');
roi_gt = [];
for d_idx = 1:numel(depth_pos)
    depth_tmp = depth_pos(d_idx);
    
    mask_map = getROIMask(axis_y, axis_z, depth_tmp+bf_.nRadius, bf_.nRadius, dor, reject_idx,1);
    
    roi_tmp = zeros(size(mask_map));
    roi_tmp(mask_map) = dsc_gt(mask_map);
    
    % mean projection
    mean_prj = mean(roi_tmp,1);
    
    roi_gt = cat(1, roi_gt, mean_prj);
end

%% GROUND TRUTH : plot 1d psf
for d_idx = 1:numel(depth_pos)
    roi_tmpp = roi_gt(d_idx,:);
    roi_tmp = roi_tmpp / max(roi_tmpp);
    
    figure(d_idx);
    plot(axis_y*1e3, roi_tmp - mean(roi_tmp), 'LineWidth',2,'color','k');
    
    title(['GROUND TRUTH - Depth: ' num2str(depth_pos(d_idx)*1e3) 'mm']);
    xlim([-15 15]); ylim([0 1]);
    set(gcf,'Position',[150 1400-300*d_idx 400 160]);
end

%% ERROR : extract 1d psf
disp('>>> extract error data');
roi_err = {};
for d_idx = 1:numel(depth_pos)
    depth_tmp = depth_pos(d_idx);
    
    mask_map = getROIMask(axis_y, axis_z, depth_tmp+bf_.nRadius, bf_.nRadius, dor, reject_idx,1);
    
    roi_err_tmp = [];
    for e_idx = 1:numel(err_case)
        err_tmp = dsc_err{e_idx};
        roi_err_tmpp = [];
        for s_idx = 1:numel(sample_list)
            err_tmpp = err_tmp(:,:,s_idx);
            
            roi_tmp = zeros(size(mask_map));
            roi_tmp(mask_map) = err_tmpp(mask_map);
            
            % mean projection
            mean_prj_tmp = mean(roi_tmp,1);
            
            roi_err_tmpp = cat(1, roi_err_tmpp, mean_prj_tmp);
        end
        roi_err_tmp = cat(3, roi_err_tmp, roi_err_tmpp);
    end
    
    roi_err = cat(3, roi_err, roi_err_tmp);
end

%% ERROR : plot 1d psf
axis_yy = axis_y * 1e3;
for d_idx = 1:numel(depth_pos)
    roi_tmppp = roi_err{d_idx};
    
    for e_idx = 1:numel(err_case)
        roi_tmpp = roi_tmppp(:,:,e_idx);
        
        roi_tmpp = roi_tmpp / max(roi_tmpp(:));
        
        mean_tmpp = mean(roi_tmpp,1);
        std_tmpp = std(roi_tmpp,0,1);
        
        figure(10*d_idx+e_idx);
        plot(axis_yy, mean_tmpp-mean(mean_tmpp),'k','LineWidth',2); hold on;
        patch([axis_yy flipud(axis_yy)], [mean_tmpp-std_tmpp-mean(mean_tmpp) flipud(mean_tmpp+std_tmpp)-mean(mean_tmpp)], 'k', 'FaceAlpha',0.2); hold off;
        
        title(err_case(e_idx).name,'Interpreter','none');
        xlim([-15 15]); ylim([0 1]);
        
        set(gcf,'Position',[150+450*e_idx 1400-300*d_idx 400 160]);
    end
end

%% XCORR2 : calculation
disp('>>> xcorr2 calculation');
xcorr_vals = zeros(numel(err_case), numel(sample_list));

for e_idx = 1:numel(err_case)
    disp(['    ' err_case(e_idx).name]);
    dsc_err_tmp = dsc_err{e_idx};
    
    for s_idx = 1:numel(sample_list)
        dsc_err_tmpp = dsc_err_tmp(:,:,s_idx);
        xcorrmap = normxcorr2(db(dsc_err_tmpp), db(dsc_gt));
        
        xcorr_vals(e_idx, s_idx) = xcorrmap((size(xcorrmap,1)+1)/2,(size(xcorrmap,2)+1)/2);
    end
end

%% XCORR2 : plot
axis_ = [0.1 0.2 0.5 1 2 5];

xcorr_mean = mean(xcorr_vals,2);
xcorr_std = std(xcorr_vals,0,2);

figure(3149);
errorbar(axis_,xcorr_mean, xcorr_std, 'LineWidth',2,'color','k');
xlabel('\sigma'); ylabel('XCORR2');

%% CONTRAST : calculation
radius_cyst = 4e-3;
offset_ = 0;
deg_ = [-18 0 18];
depth_cysts = [20 35 50]*1e-3;

bg_deg = [11 9 8];
fraction_ = 1; % REF

CNR_GT = zeros(numel(depth_cysts),1); CNR_ERR = zeros(numel(depth_cysts),numel(sample_list),numel(err_case));
Contrast_GT = zeros(numel(depth_cysts),1); Contrast_ERR = zeros(numel(depth_cysts),numel(sample_list),numel(err_case));

for c_idx = 1:numel(depth_cysts)
    disp(['>>> Cyst at ' num2str(depth_cysts(c_idx)*1e3) 'cm']);
    %% set roi region
    cyst_y = (bf_.nRadius + depth_cysts(c_idx)) * sind(deg_(3));
    cyst_z = (bf_.nRadius + depth_cysts(c_idx)) * cosd(deg_(3)) - bf_.nRadius;
    circle_cyst = [cyst_y-radius_cyst+offset_ cyst_z-radius_cyst+offset_ 2*fraction_*radius_cyst 2*fraction_*radius_cyst]*1e3;
    
    %     bg_y = (bf_.nRadius + depth_cysts(c_idx)) * sind(deg_(3)+bg_deg(c_idx));
    %     bg_z = (bf_.nRadius + depth_cysts(c_idx)) * cosd(deg_(3)+bg_deg(c_idx)) - bf_.nRadius;
    bg_y = (bf_.nRadius + depth_cysts(c_idx)) * sind(deg_(1));
    bg_z = (bf_.nRadius + depth_cysts(c_idx)) * cosd(deg_(1)) - bf_.nRadius;
    circle_bg = [bg_y-radius_cyst+offset_ bg_z-radius_cyst+offset_ 2*fraction_*radius_cyst 2*fraction_*radius_cyst]*1e3;
    %% check roi region
    figure(1);
    subplot(numel(depth_cysts),2,(c_idx-1)*2+1);
    imagesc(axis_y*1e3,(axis_z-bf_.nRadius)*1e3, db(dsc_gt/max(dsc_gt(:)))); caxis([-60 0]); hold on;
    rectangle('Position',circle_cyst,'Curvature', [1 1],'EdgeColor','y');
    rectangle('Position',circle_bg,'Curvature', [1 1],'EdgeColor','y'); hold off;
    axis tight; axis equal;
    colormap gray;
    
    %% get roi mask
    mask_map = zeros(numel(axis_z), numel(axis_y));
    
    [grid_z, grid_y] = ndgrid(axis_z, axis_y);
    
    dist_map_cyst = sqrt((grid_z-bf_.nRadius-cyst_z).^2 + (grid_y-cyst_y).^2);
    dist_map_bg = sqrt((grid_z-bf_.nRadius-bg_z).^2 + (grid_y-bg_y).^2);
    
    mask_cyst = dist_map_cyst<fraction_*radius_cyst;
    mask_bg = dist_map_bg<fraction_*radius_cyst;
    
    mask_map(mask_cyst) = 1;
    mask_map(mask_bg) = 1;
    
    figure(1);
    subplot(numel(depth_cysts),2,c_idx*2);
    imagesc(axis_y*1e3, (axis_z-bf_.nRadius)*1e3,mask_map); hold on;
    rectangle('Position',circle_cyst,'Curvature', [1 1],'EdgeColor','r');
    rectangle('Position',circle_bg,'Curvature', [1 1],'EdgeColor','r'); hold off;
    axis equal; axis tight;
    
    %% extract roi - ground truth
    cyst_env = dsc_gt(mask_cyst);
    bg_env = dsc_gt(mask_bg);
    
    CNR_GT(c_idx) = abs(mean(cyst_env) - mean(bg_env))/sqrt(std(cyst_env)^2 + std(bg_env)^2);
    Contrast_GT(c_idx) = 10*log10(sum(bg_env.^2)/sum(cyst_env.^2));
    
    %% extract roi - error
    for e_idx = 1:numel(err_case)
        err_tmp = dsc_err{e_idx};
        for s_idx = 1:numel(sample_list)
            err_tmpp = err_tmp(:,:,s_idx);
            
            cyst_env = err_tmpp(mask_cyst);
            bg_env = err_tmpp(mask_bg);
            
            CNR_ERR(c_idx, s_idx, e_idx) = abs(mean(cyst_env) - mean(bg_env))/sqrt(std(cyst_env)^2 + std(bg_env)^2);
            Contrast_ERR(c_idx, s_idx, e_idx) = 10*log10(sum(bg_env.^2)/sum(cyst_env.^2));
        end
    end
end

%% CONTRAST : plot
axiss_ = [0 0.1 0.2 0.5 1 2 5];
plot_cnr_mean = zeros(numel(depth_cysts), 1+numel(err_case));
plot_contrast_mean = zeros(numel(depth_cysts), 1+numel(err_case));
plot_cnr_std = zeros(numel(depth_cysts), 1+numel(err_case));
plot_contrast_std = zeros(numel(depth_cysts), 1+numel(err_case));

for c_idx = 1:numel(depth_cysts)
    plot_cnr_mean(c_idx,1) = CNR_GT(c_idx);
    plot_contrast_mean(c_idx,1) = Contrast_ERR(c_idx);
    
    for e_idx = 1:numel(err_case)
        cnr_tmp = CNR_ERR(:,:,e_idx);
        contrast_tmp = Contrast_ERR(:,:,e_idx);
        
        plot_cnr_mean(:,e_idx+1) = mean(cnr_tmp,2);
        plot_cnr_std(:,e_idx+1) = std(cnr_tmp,0,2);
        
        plot_contrast_mean(:,e_idx+1) = mean(contrast_tmp,2);
        plot_contrast_std(:,e_idx+1) = std(contrast_tmp,0,2);
    end
end

% plot
for c_idx = 1:numel(depth_cysts)
%     % CNR
%     figure(800+c_idx);
%     errorbar(axiss_, plot_cnr_mean(c_idx,:), plot_cnr_std(c_idx,:),'LineWidth',2);
%     set(gcf,'Position',[390+300*(c_idx-1) 600 202 299]);
    
    % Contrast
    figure(900+c_idx);
    errorbar(axiss_, plot_contrast_mean(c_idx,:), plot_contrast_std(c_idx,:),'LineWidth',2,'color','k');
    title(['Contrast @ ' num2str(depth_cysts(c_idx)*1e3) ' mm']);
%     ylim([13 35]);
    set(gcf,'Position',[390+300*(c_idx-1) 200 202 299]);
end




