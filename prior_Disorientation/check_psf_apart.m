clear; close all; clc;

%% set directory
dir_ = uigetdir('', '');

dir_data = [dir_ '/Errors'];
dir_save = [dir_ '/Saved']; mkdir(dir_save);

dir_gt = [dir_data '/error_0/Sample001/Element_64'];

gt_list = dir(dir_gt);

%% load parameter
load([dir_ '/Parameters.mat']);
stRFInfo = stParam.stRFInfo;
stBFInfo = stParam.stBFInfo;
stTRInfo = stParam.stTRInfo;
stTxInfo = stParam.stTxInfo;
mImgY = stParam.mImgY;
mImgZ = stParam.mImgZ;
aYAxis_DSC = stParam.aYAxis_DSC;
aZAxis_DSC = stParam.aZAxis_DSC;

clear stParam

%% load ground truth
name_rsaf = gt_list(end).name;
name_con = gt_list(end-1).name;

load([dir_gt '/' name_rsaf]);
img_rsaf_gt = stSaveInfo.mImg_DSC;

load([dir_gt '/' name_con]);
img_con_gt = stSaveInfo.mImg_DSC;

clear stSaveInfo;

%% Image
deg = 3.75;
line_one = [-stBFInfo.nRadius*sind(deg), -stBFInfo.nRadius*(1-cosd(deg));...
    -(stBFInfo.nRadius+stBFInfo.nDth)*sind(deg), (stBFInfo.nRadius + stBFInfo.nDth)*cosd(deg)-stBFInfo.nRadius]*1e3; % y, z

line_two = [stBFInfo.nRadius*sind(deg), -stBFInfo.nRadius*(1-cosd(deg));...
    (stBFInfo.nRadius+stBFInfo.nDth)*sind(deg), (stBFInfo.nRadius + stBFInfo.nDth)*cosd(deg)-stBFInfo.nRadius]*1e3; % y, z

figure(1);
subplot(1,2,1);imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, img_con_gt);
colormap gray; colorbar; xlabel('Elevational'); ylabel('Axial'); title('CON'); %caxis([-50 0]);
subplot(1,2,2);imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, img_rsaf_gt);
colormap gray; colorbar; xlabel('Elevational'); ylabel('Axial'); title('rSAF'); %caxis([-50 0]);
set(gcf, 'Position', [-2933 1039 1533 556]);

figure(2);
imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, img_con_gt); hold on;
line([line_one(1,1) line_one(2,1)], [line_one(1,2) line_one(2,2)], 'color', 'r', 'LineWidth', 1);
line([line_two(1,1) line_two(2,1)], [line_two(1,2) line_two(2,2)], 'color', 'r', 'LineWidth', 1); hold off;
colormap gray; colorbar; xlabel('Elevational'); ylabel('Axial'); %caxis([-50 0]);
set(gcf, 'Position', [-2287 856 788 640]);

%% load data with error
err_list = dir(dir_data);
flag = 0; if(strcmp(err_list(3).name, '.DS_Store')), flag = 1; end
err_list = err_list(4+flag:end);

disp('>>> Data loading....');
CON = {}; rSAF = {};
for e_idx = 1:numel(err_list)
    err_tmp = err_list(e_idx).name;
    
    disp(['    Error: ' err_tmp]);
    
    dir_tmp = [dir_data '/' err_tmp];
    
    sample_list = dir(dir_tmp);
    flag = 0; if(strcmp(sample_list(3).name, '.DS_Store')), flag = 1; end
    sample_list = sample_list(3+flag:end);
    
    con_ = []; rsaf_ = [];
    for s_idx = 1:numel(sample_list)
        sample_tmp = sample_list(s_idx).name;
        dir_tmpp = [dir_tmp '/' sample_tmp '/Element_64'];
        
        bf_cases = dir(dir_tmpp);
        con_tmp = bf_cases(end-1).name;
        rsaf_tmp = bf_cases(end).name;
        
        load([dir_tmpp '/' con_tmp]);
        con_ = cat(3, con_, stSaveInfo.mImg_DSC);
        load([dir_tmpp '/' rsaf_tmp]);
        rsaf_ = cat(3, rsaf_, stSaveInfo.mImg_DSC);
    end
    
    CON = cat(3, CON, con_);
    rSAF = cat(3, rSAF, rsaf_);
end

%% get mask for each radial depth
depth_ = (3:7) * 1e-2; % radial direction
% offset = 0.03e-3;% 3cm
% offset = 0.09e-3;
% depth = depth_ + offset;
depth = depth_;
% depth(2:5) = depth(2:5);% + 0.1e-3;
% depth = depth_;
depth(1) = depth_(1) + 0.03e-3;
depth(4) = depth_(4) + 0.1e-3;

y_axis = aYAxis_DSC;
z_axis = aZAxis_DSC - stBFInfo.nRadius;

[z_grid, y_grid] = ndgrid(z_axis, y_axis);

dist_map = sqrt(z_grid.^2 + y_grid.^2);

mask_map = logical(zeros(size(dist_map)));
for d_idx = 1:numel(depth)
    for k = 1:numel(y_axis)
        tmp = dist_map(:,k);
        idx = find(abs(tmp - depth(d_idx)) == min(abs(tmp-depth(d_idx))));
        mask_map(idx(end), k) = true;
    end
end

reject_area = (img_con_gt == 0);
mask_map(reject_area) = false;
figure(3);imagesc(mask_map); title('mask map');


% ground truth compare
con_1d_gt = {};
rsaf_1d_gt = {};
for d_idx = 1:numel(depth)
    % create 1-D mask map for radial depth
    mask_map = logical(zeros(size(dist_map)));
    for k = 1:size(dist_map,2)
        tmp = dist_map(:,k);
        idx = find(abs(tmp - depth(d_idx)) == min(abs(tmp-depth(d_idx))));
        mask_map(idx(1), k) = true;
    end
%     idx = find(abs(z_axis - depth(d_idx)) == min(abs(z_axis - depth(d_idx))));
%     mask_map(idx(1),:) = true;
    
    con_1d = img_con_gt(mask_map);
    con_1d = con_1d(con_1d ~= 50); % only take imaging region
%     con_1d = con_1d + abs(max(con_1d)); % normalize
    
    rsaf_1d = img_rsaf_gt(mask_map);
    rsaf_1d = rsaf_1d(rsaf_1d ~= 50); % only take imaging region
%     rsaf_1d = rsaf_1d + abs(max(rsaf_1d)); % normalize
   
    radial_axis = linspace(-0.5*stBFInfo.nFOV, 0.5*stBFInfo.nFOV, numel(con_1d));
    
    figure(4); subplot(numel(depth), 1, d_idx);
    plot(radial_axis, con_1d, 'LineWidth', 2, 'color', 'k'); hold on;
    plot(radial_axis, rsaf_1d, 'LineWidth', 2, 'color', 'b'); hold off;
    legend('CON', 'rSAF');title(['depth: ' num2str(depth_(d_idx)*1e2) 'cm']);
    xlabel('Radial axis[deg]'); ylabel('PSF'); 
    xlim([-10 10]); %ylim([-30 0]);
    
    con_1d_gt = cat(3, con_1d_gt, con_1d);
    rsaf_1d_gt = cat(3, rsaf_1d_gt, rsaf_1d);
end
set(gcf, 'Position', [1830 140 1135 954]);

%% ground truth image - at each depth
% set roi
% Coordinate setting
nNoPT = 7;
nInitialPosY = 0e-3;
nInitialPosZ = 10e-3;
nSpacing = 10e-3;

% ROI position
aTargetPosY = linspace(nInitialPosY, nInitialPosY, nNoPT);
aTargetPosZ = linspace(nInitialPosZ, nInitialPosZ+(nNoPT-1)*nSpacing, nNoPT);

nROISize_Hor = 20e-3;
nROISize_Ver = 4e-3;

mROIPosY = zeros(nNoPT, 2);
mROIPosZ = zeros(nNoPT,2);

mROIPosY(:,1) = aTargetPosY - 0.5*nROISize_Hor; % left
mROIPosY(:,2) = aTargetPosY + 0.5*nROISize_Hor; % right

mROIPosZ(:,1) = aTargetPosZ - 0.5*nROISize_Ver; % upper
mROIPosZ(:,2) = aTargetPosZ + 0.5*nROISize_Ver; % bottom

for p_idx = 3:nNoPT
    nLft = mROIPosY(p_idx,1);
    nRgt = mROIPosY(p_idx,2);
    nUp = mROIPosZ(p_idx,1);
    nDn = mROIPosZ(p_idx,2);
    
    lIdx = find(abs(aYAxis_DSC-nLft)==min(abs(aYAxis_DSC-nLft)));
    rIdx = find(abs(aYAxis_DSC-nRgt)==min(abs(aYAxis_DSC-nRgt)));
    uIdx = find(abs(aZAxis_DSC-stBFInfo.nRadius-nUp)==min(abs(aZAxis_DSC-stBFInfo.nRadius-nUp)));
    dIdx = find(abs(aZAxis_DSC-stBFInfo.nRadius-nDn)==min(abs(aZAxis_DSC-stBFInfo.nRadius-nDn)));
    
    y_roi = aYAxis_DSC(lIdx:rIdx);
    z_roi = aZAxis_DSC(uIdx:dIdx);
    
    % create 1-D mask map for radial depth
    mask_map = logical(zeros(size(dist_map)));
    for k = 1:numel(y_axis)
        tmp = dist_map(:,k);
        idx = find(abs(tmp - depth(p_idx-2)) == min(abs(tmp-depth(p_idx-2))));
        mask_map(idx(1), k) = true;
    end
%     idx = find(abs(z_axis - depth(p_idx-2)) == min(abs(z_axis - depth(p_idx-2))));
%     mask_map(idx(1),:) = true;
    
    mask_roi = mask_map(uIdx:dIdx, lIdx:rIdx);
    %%%%%% CON
    roi_con = img_con_gt(uIdx:dIdx, lIdx:rIdx);
    idx_roi = find(roi_con ~= 0);
%     roi_con(idx_roi) = roi_con(idx_roi) + abs(max(roi_con(idx_roi))); % normalize
    roi_con(idx_roi) = roi_con(idx_roi)/abs(max(roi_con(idx_roi)));
    
    %%%%%% rsaf
    roi_rsaf = img_rsaf_gt(uIdx:dIdx, lIdx:rIdx);
    idx_roi = find(roi_con ~= 0);
%     roi_rsaf(idx_roi) = roi_rsaf(idx_roi) + abs(max(roi_rsaf(idx_roi))); % normalize
    roi_rsaf(idx_roi) = roi_rsaf(idx_roi)/abs(max(roi_rsaf(idx_roi)));
    
    figure(5);
    subplot(5, 2, (p_idx-1-2)*2+1);
    imagesc(y_roi*1e3, (z_roi-stBFInfo.nRadius)*1e3, roi_con); hold on;
    [ContourLine_bg, h] = contour(y_roi*1e3, (z_roi-stBFInfo.nRadius)*1e3, mask_roi, [1 1],'ShowText','off', 'LineColor','r', 'LineWidth', 1);
    title(['Depth: ' num2str(p_idx) 'cm | CON']); colormap default; 
%     caxis([-50 0]); 
    hold off; colormap gray;
    subplot(5, 2, (p_idx-2)*2);
    imagesc(y_roi*1e3, (z_roi-stBFInfo.nRadius)*1e3, roi_rsaf); hold on;
    [ContourLine_bg, h] = contour(y_roi*1e3, (z_roi-stBFInfo.nRadius)*1e3, mask_roi, [1 1],'ShowText','off', 'LineColor','r', 'LineWidth', 1);
    title(['Depth: ' num2str(p_idx) 'cm | rSAF']); colormap default; 
%     caxis([-50 0]); 
    hold off; colormap gray;
end
set(gcf, 'Position', [1830 140 1135 954]);

%% error case
% take mean of 100 samples
error_CON = {}; error_rSAF = {};
for d_idx = 1:numel(depth) % at each depth
    disp(['>>> depth: ' num2str(depth(d_idx)*1e2) 'cm']);
    mask_map = logical(zeros(size(dist_map)));
    
    % make mask at corresponding depth
    for k = 1:numel(y_axis)
        tmp = dist_map(:,k);
        idx = find((abs(tmp - depth(d_idx)) == min(abs(tmp-depth(d_idx)))));
        mask_map(idx(1),k) = true;
    end
%     idx = find(abs(z_axis - depth(d_idx)) == min(abs(z_axis - depth(d_idx))));
%     mask_map(idx(1),:) = true;
    
    con_mean = []; rsaf_mean = [];
    for e_idx = 1:numel(err_list)
        con_tmp = cell2mat(CON(:,:,e_idx));
        rsaf_tmp = cell2mat(rSAF(:,:,e_idx));
        
        con_samples = []; rsaf_samples = [];
        for s_idx = 1:numel(sample_list) % at each sample
            con_tmpp = con_tmp(:,:,s_idx);
            rsaf_tmpp = rsaf_tmp(:,:,s_idx);
            
            con_1d = con_tmpp(mask_map);
%                         idx = (con_1d ~= 0);
%                         con_1d(idx) = con_1d(idx) + abs(max(con_1d(idx)));
                        con_1d = con_1d(con_1d ~= 50); % only take the imaging region
%                         con_1d = con_1d + abs(max(con_1d)); % normalize
            
            rsaf_1d = rsaf_tmpp(mask_map);
            %             idx = (con_1d ~= -30);
            %             rsaf_1d(idx) = rsaf_1d(idx) + abs(max(rsaf_1d(idx)));
                        rsaf_1d = rsaf_1d(rsaf_1d ~= 50); % only take the imaging region
            %             rsaf_1d = rsaf_1d + abs(max(rsaf_1d)); % normalize
            
            con_samples = cat(2, con_samples, con_1d);
            rsaf_samples = cat(2, rsaf_samples, rsaf_1d);
        end
        con_mean = cat(2, con_mean, mean(con_samples, 2));
        rsaf_mean = cat(2, rsaf_mean, mean(rsaf_samples, 2));
    end
    
    error_CON = cat(3, error_CON, con_mean);
    error_rSAF = cat(3, error_rSAF, rsaf_mean);
end

%% specific depth
depth_roi = 3;
d_idx = depth_roi - (depth_(1)*1e2-1);

% con_ = cell2mat(error_CON(:,:,d_idx));
% rsaf_ = cell2mat(error_rSAF(:,:,d_idx));

for e_idx = 1:numel(err_list)
    con_gt_tmp = cell2mat(con_1d_gt(:,:,d_idx));
    rsaf_gt_tmp = cell2mat(rsaf_1d_gt(:,:,d_idx));

    assert(numel(con_gt_tmp) == numel(rsaf_gt_tmp) , 'Dimension conflicts');
    radial_axis = linspace(-0.5*stBFInfo.nFOV, 0.5*stBFInfo.nFOV, numel(con_gt_tmp));
    
    % ground truth
    figure(1);
    plot(radial_axis, con_gt_tmp, 'LineWidth', 2, 'color', [0 0 0]); hold on;
    figure(2);
    plot(radial_axis, rsaf_gt_tmp, 'LineWidth', 2, 'color', [0 0 0]); hold on;
    
    % error
    % normalize
    con_norm = con_(:,e_idx);
    con_norm = con_norm(con_norm ~= 0);
    con_norm = con_norm + abs(max(con_norm));
    
    rsaf_norm = rsaf_(:,e_idx);
    rsaf_norm = rsaf_norm(rsaf_norm ~= 0);
    rsaf_norm = rsaf_norm + abs(max(rsaf_norm));
    
    figure(1);
    plot(radial_axis, con_norm, 'LineWidth', 2, 'LineStyle', '-.', 'color', [0.12 0.12 0.12]*e_idx);
    figure(2);
    plot(radial_axis, rsaf_norm, 'LineWidth', 2, 'LineStyle', '-.', 'color', [0.12 0.12 0.12]*e_idx); 
end
figure(1);
grid on; hold off;
title(['CON | Depth: ' num2str(depth_(d_idx)*1e2) 'cm | ' err_list(e_idx).name], 'Interpreter', 'None');
% xlim([-5 5]);
ylim([0 5e-23]);
xlabel('radial axis [deg]'); ylabel('PSF [dB]');
legend('GT', '\sigma = 0.1\circ', '\sigma = 0.2\circ', '\sigma = 0.5\circ',...
        '\sigma = 1\circ', '\sigma = 2\circ', '\sigma = 5\circ');
% set(gcf, 'Position', [1824 628 1080 413]);


figure(2);
grid on; hold off;
title(['rSAF | Depth: ' num2str(depth_(d_idx)*1e2) 'cm | ' err_list(e_idx).name], 'Interpreter', 'None');
% xlim([-5 5]);
ylim([0 5e-23]);
xlabel('radial axis [deg]'); ylabel('PSF [dB]');
legend('GT', '\sigma = 0.1\circ', '\sigma = 0.2\circ', '\sigma = 0.5\circ',...
        '\sigma = 1\circ', '\sigma = 2\circ', '\sigma = 5\circ');
% set(gcf, 'Position', [1824 628 1080 413]);
