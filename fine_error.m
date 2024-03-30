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
dsc_gt = [];
for k = 1
    load([dir_gt_ '/' bf_cases(end).name]);
    env_data = mid_proc(stSaveInfo.mBFedData, mid_, acoustic_, bf_);
    [axis_y, axis_z, mOutput, reject_idx] = dsc(env_data, dr, da, bf_, height, width, dz, dy);
    
    dsc_gt = cat(3, dsc_gt, mOutput);
    
    figure(1);
    imagesc(axis_y*1e3, (axis_z-bf_.nRadius)*1e3, db(mOutput/max(mOutput(:))));
    axis equal; axis tight; colormap gray; title(bf_cases(k).name,'Interpreter','None');
    caxis([-50 0]);
end
set(gcf,'Position',[100 100 330 392]);
disp('done'); pause(1);

%% ERROR : load & dsc
disp('>>> load error data');
dsc_err_rSAF = {};

for e_idx = 1:numel(err_case)
    disp(['    ' err_case(e_idx).name]);
    dir_tmp = [dir_data '/' err_case(e_idx).name];
    
    sample_list = dir(dir_tmp);
    flag = 0;   
    if(strcmp(sample_list(3).name, '.DS_Store')), flag = 1; end
    sample_list = sample_list(3+flag:end);
    
    dsc_tmp_rSAF = [];
    for s_idx = 1:numel(sample_list)
        dir_tmpp = [dir_tmp '/' sample_list(s_idx).name '/Element_64'];
        
        % load TRUS-rSAF
        load([dir_tmpp '/' bf_cases(end).name]);
        env_data = mid_proc(stSaveInfo.mBFedData, mid_, acoustic_, bf_);
        [axis_y, axis_z, mOutput, reject_idx] = dsc(env_data, dr, da, bf_, height, width, dz, dy);
        
        dsc_tmp_rSAF = cat(3, dsc_tmp_rSAF, mOutput);
    end
    
    dsc_err_rSAF = cat(3, dsc_err_rSAF, dsc_tmp_rSAF);
end
disp('done');
%% SET MASK REGION
depth = 50e-3;
dor = 5e-3;
mask_map = getROIMask(axis_y, axis_z, depth+bf_.nRadius, bf_.nRadius, dor, reject_idx,1);

sampleROI = dsc_gt .* mask_map;
sampleROI(reject_idx) = 1e-25;

figure;
imagesc(axis_y*1e3, (axis_z-bf_.nRadius)*1e3, db(sampleROI/max(dsc_gt(:))));
colormap gray; axis tight; axis equal;
caxis([-50 0]);
%% GROUND TRUTH : extract 1d psf (mean projection)
offset_ = 0.1; % 60mm
% offset_ = 0.06; % 30mm

gt_psf = [];
for k = 1
    psf_tmp = mean(dsc_gt.*mask_map,1);
    
    gt_psf = cat(1, gt_psf, psf_tmp);
end

norm_val = max(gt_psf(1,:)); % Maximum value of TRUS-rSAF

figure(322141551);
for k = 1
    gt_tmp = gt_psf(k,:)/norm_val;
    plot(atand(axis_y/depth), gt_tmp - mean(gt_tmp) + offset_,'LineWidth',2);hold on;
end
hold off;
title(['GROUND TRUTH - Depth: ' num2str(depth*1e3) 'mm']);
xlim([-65 65]); ylim([0 1]);
legend('TRUS-CON', 'TRUS-rSAF');
set(gcf,'Position',[1483 456 832 339]);

%% ERROR : extract 1d psf
disp(['>>> Depth : ' num2str(depth*1e3) 'mm']);
err_psf_rSAF = [];

for e_idx = 1:numel(err_case)
    dsc_tmp_rSAF = dsc_err_rSAF{e_idx};
    
    err_psf_rSAF_tmp = [];
    for s_idx = 1
        dsc_tmp = dsc_tmp_rSAF(:,:,s_idx);
        psf_tmp = squeeze(mean(dsc_tmp.*mask_map,1));
        
        err_psf_rSAF_tmp = cat(1, err_psf_rSAF_tmp, psf_tmp);
    end
    err_psf_rSAF = cat(1, err_psf_rSAF, err_psf_rSAF_tmp);
end

% ERROR : plot 1d psf
offset_ = 0.1; % 60mm
% offset_ = 0.04; % 30mm

leg_ = {};

figure(112314);
plot(atand(axis_y/depth), db(gt_tmp),'LineWidth',3,'color','k');hold on;
leg_ = cat(1, leg_, 'GROUND TRUTH');

for e_idx = 1:numel(err_case)
    err_psf_rSAF_tmp = err_psf_rSAF(e_idx,:);
    
    mean_psf_rSAF_norm = err_psf_rSAF_tmp/norm_val;
    
    plot(atand(axis_y/depth), db(mean_psf_rSAF_norm),'LineWidth',3);
    
    leg_ = cat(1, leg_, err_case(e_idx).name);
end

hold off;
legend(leg_,'Interpreter','none');
xlim([-65 65]); ylim([-50 0]);
set(gcf,'Position',[410 560 832 339]);

%% subtraction from ground truth
legg_ = {};

figure(112315);

for e_idx = 1:numel(err_case)
    err_psf_rSAF_tmp = err_psf_rSAF(e_idx,:);
    
    mean_psf_rSAF_norm = err_psf_rSAF_tmp/norm_val;
    
    subtractedTmp = db(mean_psf_rSAF_norm) - db(gt_tmp);
    
    plot(atand(axis_y/depth), (subtractedTmp),'LineWidth',3); hold on;
    
    legg_ = cat(1, legg_, err_case(e_idx).name);
end

hold off;
legend(legg_,'Interpreter','none');
xlim([-65 65]); ylim([0 20]);
set(gcf,'Position',[410 100 832 339]);

%% b-mode subtraction
figure(3215);
imagesc(axis_y*1e3, (axis_z-bf_.nRadius)*1e3, db(dsc_gt/max(dsc_gt(:))));
axis equal; axis tight; colormap gray; title(err_case(e_idx).name,'Interpreter','None');
caxis([-50 0]);
set(gcf,'Position',[100 100 330 392]);
ylim([0 70]); xlim([-35 35]);

for e_idx = 1:numel(err_case)
    errTmp = dsc_err_rSAF{e_idx};
    subTmp = dsc_gt - errTmp(:,:,1);
    
    figure(e_idx);
    subplot(1,2,1);
    imagesc(axis_y*1e3, (axis_z-bf_.nRadius)*1e3, db(errTmp/max(dsc_gt(:))));
    axis equal; axis tight; colormap gray;
    caxis([-50 0]);
    ylim([0 70]); xlim([-35 35]);
    subplot(1,2,2);
    imagesc(axis_y*1e3, (axis_z-bf_.nRadius)*1e3, db(subTmp/max(dsc_gt(:))));
    axis equal; axis tight; colormap gray;
    caxis([-50 0]);
    set(gcf,'Position',[440 100 760 392]);
    ylim([0 70]); xlim([-35 35]);
end

%% axis setting
axis_ = [0.01 0.02 0.05 0.1 0.2 0.5];
%% SSIM - entire image
ssimValsEntire_ = zeros(numel(err_case), numel(sample_list));

dscGTNorm = dsc_gt / max(dsc_gt(:));

for e_idx = 1:numel(err_case)
    disp(err_case(e_idx).name);
    err_rSAF_tmp = dsc_err_rSAF{e_idx}/max(dsc_gt(:));
    for s_idx = 1:numel(sample_list)
        dscErrTmp = err_rSAF_tmp(:,:,s_idx);
        
        [ssimVal, ssimMap] = ssim(db(dscGTNorm),db(dscErrTmp));
        
        ssimValsEntire_(e_idx,s_idx) = ssimVal;
    end
end
%%
ssimValsMean = mean(ssimValsEntire_,2);
ssimValsSTD = std(ssimValsEntire_,0,2);

figure(2031040);
errorbar(axis_,ssimValsMean,ssimValsSTD,'LineWidth',2);
xlabel('\sigma'); ylabel('SSIM');

