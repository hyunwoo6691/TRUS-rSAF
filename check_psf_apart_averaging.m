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
dir_data = [dir_ '/errors_bf_AWGN_20'];

total_case = dir(dir_data);
flag = 0;
if(strcmp(total_case(3).name,'.DS_Store')), flag = 1; end
total_case = total_case(3+flag:end);

dir_gt = [dir_data '/' total_case(1).name];

err_case = total_case(2:end);

%% GROUND TRUTH : show image
dir_gt_ = [dir_gt '/Sample001/Element_64'];
bf_cases = dir(dir_gt_);
flag = 0;
if(strcmp(bf_cases(3).name, '.DS_Store')), flag = 1; end
bf_cases = bf_cases(3+flag:end);

disp('>>> load ground truth data');
dsc_gt = [];
for k = 1:2
    load([dir_gt_ '/' bf_cases(k).name]);
    env_data = mid_proc(stSaveInfo.mBFedData, mid_, acoustic_, bf_);
    [axis_y, axis_z, mOutput, reject_idx] = dsc(env_data, dr, da, bf_, height, width, dz, dy);
    
    dsc_gt = cat(3, dsc_gt, mOutput);
    
    figure(1);
    subplot(1,2,k);
    imagesc(axis_y*1e3, (axis_z-bf_.nRadius)*1e3, db(mOutput/max(mOutput(:))));
    axis equal; axis tight; colormap gray; title(bf_cases(k).name,'Interpreter','None');
    caxis([-50 0]);
end
set(gcf, 'Position', [561 550 1069 397]);
disp('done'); pause(1);

%% averaging
load([dir_gt_ '/' bf_cases(1).name]);
[beamformed_data_fil, fil] = spatial_filtering(stSaveInfo.mBFedData, 0, bf_.nDTheta, 'boxcar');
env_data = mid_proc(beamformed_data_fil, mid_, acoustic_, bf_);
[axis_y, axis_z, mOutput, reject_idx] = dsc(env_data, dr, da, bf_, height, width, dz, dy);

dsc_gt = cat(3, dsc_gt, mOutput);

%% ERROR : load & dsc
disp('>>> load error data');
dsc_err_CON = {}; dsc_err_rSAF = {}; 
dsc_err_CON_fil = {};

sigma_case = [0.1 0.2 0.5 1 2 5];

for e_idx = 1:numel(err_case)
    disp(['    ' err_case(e_idx).name]);
    dir_tmp = [dir_data '/' err_case(e_idx).name];
    sample_list = dir(dir_tmp);
    
    flag = 0;
    if(strcmp(sample_list(3).name,'.DS_Store')), flag = 1; end
    sample_list = sample_list(3+flag:end);
    
    dsc_tmp_CON = []; dsc_tmp_rSAF = []; 
    dsc_tmp_CON_fil = [];
    for s_idx = 1:numel(sample_list)
        if(mod(s_idx, numel(sample_list)/4) == 0), disp(['      ' num2str(round(100*s_idx/numel(sample_list))) '%...']); end
        dir_tmpp = [dir_tmp '/' sample_list(s_idx).name '/Element_64'];
        
        % load TRUS-CON
        load([dir_tmpp '/' bf_cases(1).name]);
        env_data = mid_proc(stSaveInfo.mBFedData, mid_, acoustic_, bf_);
        [axis_y, axis_z, mOutput, reject_idx] = dsc(env_data, dr, da, bf_, height, width, dz, dy);
%         
        dsc_tmp_CON = cat(3, dsc_tmp_CON, mOutput);
        
        % average filter
        [beamformed_data_fil, fil] = spatial_filtering(stSaveInfo.mBFedData, sigma_case(e_idx), bf_.nDTheta, 'boxcar');
        env_data = mid_proc(beamformed_data_fil, mid_, acoustic_, bf_);
        [axis_y, axis_z, mOutput, reject_idx] = dsc(env_data, dr, da, bf_, height, width, dz, dy);
        
        dsc_tmp_CON_fil = cat(3, dsc_tmp_CON_fil, mOutput);
        
        % load TRUS-rSAF
        load([dir_tmpp '/' bf_cases(2).name]);
        env_data = mid_proc(stSaveInfo.mBFedData, mid_, acoustic_, bf_);
        [axis_y, axis_z, mOutput, reject_idx] = dsc(env_data, dr, da, bf_, height, width, dz, dy);
%         
        dsc_tmp_rSAF = cat(3, dsc_tmp_rSAF, mOutput);
    end
    
    dsc_err_CON = cat(3, dsc_err_CON, dsc_tmp_CON);
    dsc_err_rSAF = cat(3, dsc_err_rSAF, dsc_tmp_rSAF);
    dsc_err_CON_fil = cat(3, dsc_err_CON_fil, dsc_tmp_CON_fil);
end
disp('done');
%% SET MASK REGION
depth = 30e-3;
dor = 1e-3;
mask_map = getROIMask(axis_y, axis_z, depth+bf_.nRadius, bf_.nRadius, dor, reject_idx,1);

%% BOTH : close up 2d-psf
dim_ = [30 5]*1e-3; % horizontal, vertical
center_ = [0, depth]; % horizontal, vertial

% get axis
range_hor = [center_(1)-0.5*dim_(1) center_(1)+0.5*dim_(1)];
range_ver = [center_(2)-0.5*dim_(2) center_(2)+0.5*dim_(2)];

idx_y = [find(abs(range_hor(1)-axis_y) == min(abs(range_hor(1)-axis_y))) find(abs(range_hor(2)-axis_y) == min(abs(range_hor(2)-axis_y)))];
idx_z = [find(abs(range_ver(1)-axis_z+bf_.nRadius) == min(abs(range_ver(1)-axis_z+bf_.nRadius))) find(abs(range_ver(2)-axis_z+bf_.nRadius) == min(abs(range_ver(2)-axis_z+bf_.nRadius)))];

axis_yy = axis_y(idx_y(1):idx_y(2));
axis_zz = axis_z(idx_z(1):idx_z(2)) - bf_.nRadius;

% ground truth
roi_gt_CON = dsc_gt(idx_z(1):idx_z(2),idx_y(1):idx_y(2),1);
roi_gt_rSAF = dsc_gt(idx_z(1):idx_z(2),idx_y(1):idx_y(2),2);
roi_gt_CON_fil = dsc_gt(idx_z(1):idx_z(2),idx_y(1):idx_y(2),3);

figure(1231); imagesc(axis_yy,axis_zz,db(roi_gt_CON/max(max(dsc_gt(:,:,1))))); colormap gray; caxis([-60 0]); axis tight; axis equal; set(gcf,'Position',[30 1237 410 90]);
figure(1232); imagesc(axis_yy,axis_zz,db(roi_gt_rSAF/max(max((dsc_gt(:,:,2)))))); colormap gray; caxis([-60 0]); axis tight; axis equal; set(gcf,'Position',[450 1237 410 90]);
figure(1233); imagesc(axis_yy,axis_zz,db(roi_gt_CON_fil/max(max((dsc_gt(:,:,3)))))); colormap gray; caxis([-60 0]); axis tight; axis equal; set(gcf,'Position',[870 1237 410 90]);

% error
indices_ = [1 1 1 1 1 1];
for e_idx = 1:numel(err_case)
    dsc_err_CON_tmp = dsc_err_CON{e_idx};
    dsc_err_rSAF_tmp = dsc_err_rSAF{e_idx};
    dsc_err_CON_tmp_fil = dsc_err_CON_fil{e_idx};
    
    roi_CON_tmp = dsc_err_CON_tmp(idx_z(1):idx_z(2),idx_y(1):idx_y(2),indices_(e_idx));
    roi_rSAF_tmp = dsc_err_rSAF_tmp(idx_z(1):idx_z(2),idx_y(1):idx_y(2),indices_(e_idx));
    roi_CON_tmp_fil = dsc_err_CON_tmp_fil(idx_z(1):idx_z(2),idx_y(1):idx_y(2),indices_(e_idx));
    
    figure(12314+e_idx); % CON
    imagesc(axis_yy,axis_zz,db(roi_CON_tmp/max(dsc_err_CON_tmp(:)))); colormap gray; caxis([-60 0]); axis tight; axis equal; set(gcf,'Position',[30 1237-200*e_idx 410 90]);
    figure(123140+e_idx); % rSAF
    imagesc(axis_yy,axis_zz,db(roi_rSAF_tmp/max(dsc_err_rSAF_tmp(:)))); colormap gray; caxis([-60 0]); axis tight; axis equal; set(gcf,'Position',[450 1237-200*e_idx 410 90]);
    figure(1231400+e_idx); % CON-Average
    imagesc(axis_yy,axis_zz,db(roi_CON_tmp_fil/max(dsc_err_CON_tmp_fil(:)))); colormap gray; caxis([-60 0]); axis tight; axis equal; set(gcf,'Position',[870 1237-200*e_idx 410 90]);
end
%% GROUND TRUTH : extract 1d psf (mean projection)
% offset_ = 0.1; % 60mm
% offset_ = 0.06; % 30mm
% offset_ = 0.05;
offset_ = 0;

gt_psf = [];
for k = 1:size(dsc_gt,3)
    dsc_tmp = dsc_gt(:,:,k);
    psf_tmp = mean(dsc_tmp.*mask_map,1);
    
    gt_psf = cat(1, gt_psf, psf_tmp);
end

norm_val = max(gt_psf(1,:)); % Maximum value of TRUS-rSAF
figure(3221*depth*1e3); 
for k = 1:size(dsc_gt,3)
    gt_tmp = gt_psf(k,:)/norm_val;
    plot(atand(axis_y/depth), gt_tmp - mean(gt_tmp) + offset_,'LineWidth',3);hold on;
%     plot(atand(axis_y/depth), gt_tmp - mean(gt_tmp), 'LineWidth',3);hold on;
end
hold off; 
title(['GROUND TRUTH - Depth: ' num2str(depth*1e3) 'mm']);
xlim([-30 30]); ylim([0 1]);
legend('TRUS-CON', 'TRUS-rSAF');
set(gcf,'Position',[1483 456 832 339]);

%% ERROR : extract 1d psf
disp(['>>> Depth : ' num2str(depth*1e3) 'mm']);
err_psf_CON = []; err_psf_rSAF = []; err_psf_CON_fil = [];

for e_idx = 1:numel(err_case)
    dsc_tmp_CON = dsc_err_CON{e_idx};
    dsc_tmp_rSAF = dsc_err_rSAF{e_idx};
    dsc_tmp_CON_fil = dsc_err_CON_fil{e_idx};
    
    err_psf_CON_tmp = []; err_psf_rSAF_tmp = []; err_psf_CON_fil_tmp = [];
    for s_idx = 1:numel(sample_list)
        dsc_tmp = dsc_tmp_CON(:,:,s_idx);
        psf_tmp = mean(dsc_tmp.*mask_map,1);
        
        err_psf_CON_tmp = cat(1, err_psf_CON_tmp, psf_tmp);
        
        dsc_tmp = dsc_tmp_rSAF(:,:,s_idx);
        psf_tmp = mean(dsc_tmp.*mask_map,1);
        
        err_psf_rSAF_tmp = cat(1, err_psf_rSAF_tmp, psf_tmp);
        
        dsc_tmp = dsc_tmp_CON_fil(:,:,s_idx);
        psf_tmp = mean(dsc_tmp.*mask_map,1);
        
        err_psf_CON_fil_tmp = cat(1, err_psf_CON_fil_tmp, psf_tmp);
    end
    err_psf_CON = cat(3, err_psf_CON, err_psf_CON_tmp);
    err_psf_rSAF = cat(3, err_psf_rSAF, err_psf_rSAF_tmp);
    err_psf_CON_fil = cat(3, err_psf_CON_fil, err_psf_CON_fil_tmp);
end

%% ERROR : plot 1d psf
% offset_ = 0.05; % 60mm
% offset_ = 0.04; % 30mm
offset_ = 0;

for e_idx = 1:numel(err_case)
    err_psf_CON_tmp = err_psf_CON(:,:,e_idx);
    err_psf_rSAF_tmp = err_psf_rSAF(:,:,e_idx);
    err_psf_CON_fil_tmp = err_psf_CON_fil(:,:,e_idx);
    
    mean_psf_CON = mean(err_psf_CON_tmp);
    mean_psf_rSAF = mean(err_psf_rSAF_tmp);
    mean_psf_CON_fil = mean(err_psf_CON_fil_tmp);
    
    mean_psf_CON_norm = mean_psf_CON/norm_val;
    mean_psf_rSAF_norm = mean_psf_rSAF/norm_val;
    mean_psf_CON_fil_norm = mean_psf_CON_fil/norm_val;
    
    figure(331*depth*1e3+e_idx); hold on;
    for s_idx = 1:numel(sample_list)
        err_CON_tmp_norm = err_psf_CON_tmp(s_idx,:)/norm_val;
        plot(atand(axis_y/depth),err_CON_tmp_norm-mean(err_CON_tmp_norm)+offset_,'color',[0.12 0.12 0.12]*7); 
    end
    
    for s_idx = 1:numel(sample_list)
        err_rSAF_tmp_norm = err_psf_rSAF_tmp(s_idx,:)/norm_val;
        plot(atand(axis_y/depth),err_rSAF_tmp_norm-mean(err_rSAF_tmp_norm)+offset_,'color',[0.58 0.82 0.99]);
    end
    
    for s_idx = 1:numel(sample_list)
        err_CON_fil_tmp_norm = err_psf_CON_fil_tmp(s_idx,:)/norm_val;
        plot(atand(axis_y/depth),err_CON_fil_tmp_norm-mean(err_CON_fil_tmp_norm)+offset_,'color',[0.47 0.67 0.19]); 
    end
    
    plot(atand(axis_y/depth), mean_psf_CON_norm - mean(mean_psf_CON_norm)+offset_,'LineWidth',3,'color','k'); 
    plot(atand(axis_y/depth), mean_psf_rSAF_norm - mean(mean_psf_rSAF_norm)+offset_,'LineWidth',3,'color','b');
    plot(atand(axis_y/depth), mean_psf_CON_fil_norm - mean(mean_psf_CON_fil_norm)+offset_,'LineWidth',3,'color','g'); hold off;
    
    title(err_case(e_idx).name,'Interpreter','None');
    xlim([-30 30]); ylim([0 1]);
    
    set(gcf,'Position',[1 456 832 339]);
end

%% axis setting
axis_ = [0.1 0.2 0.5 1 2 5];

%% EROR : Coefficient of variance
CV_CON = zeros(1, numel(err_case));
CV_rSAF = zeros(1, numel(err_case));
CV_CON_fil = zeros(1, numel(err_case));

figure(13121209);
plot(gt_tmp);

% idx_ = 361; % 60mm
idx_ = 390; % 30mm

for e_idx = 1:numel(err_case)
    err_psf_CON_tmp = err_psf_CON(:,:,e_idx);
    err_psf_rSAF_tmp = err_psf_rSAF(:,:,e_idx);
    err_psf_CON_fil_tmp = err_psf_CON_fil(:,:,e_idx);
    
    CV_CON(e_idx) = std(err_psf_CON_tmp(:,idx_))/mean(err_psf_CON_tmp(:,idx_));
    CV_rSAF(e_idx) = std(err_psf_rSAF_tmp(:,idx_))/mean(err_psf_rSAF_tmp(:,idx_));
    CV_CON_fil(e_idx) = std(err_psf_CON_fil_tmp(:,idx_))/mean(err_psf_CON_fil_tmp(:,idx_));
end
figure(13125);
plot(axis_, CV_CON); hold on;
plot(axis_, CV_rSAF);
plot(axis_, CV_CON_fil); hold off;
set(gcf,'Position',[1868 326 396 337]);





