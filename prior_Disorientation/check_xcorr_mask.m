clc; clear; close all;

%% set directory, load paramter & ground truth images
dir_ = uigetdir('','');

load([dir_ '/Parameters.mat']);
stBFInfo = stParam.stBFInfo;
aYAxis_DSC = stParam.aYAxis_DSC;
aZAxis_DSC = stParam.aZAxis_DSC;
clear stParam;

dir_err = [dir_ '/Errors'];
err_list = dir(dir_err);
flag = 0;
if(strcmp(err_list(3).name, '.DS_Store')), flag = 1; end
err_list = err_list(4+flag:end);

dir_gt = [dir_err '/' err_list(3+flag).name '/Sample001/Element_64'];
img_list = dir(dir_gt);

load([dir_gt '/' img_list(end-1).name]);
img_gt_con = stSaveInfo.mImg_DSC;
load([dir_gt '/' img_list(end).name]);
img_gt_rsaf = stSaveInfo.mImg_DSC;


%% set roi
% Coordinate setting
nNoPT = 7;
nInitialPosY = 0e-3;
nInitialPosZ = 10e-3;
nSpacing = 10e-3;

% ROI position
aTargetPosY = linspace(nInitialPosY, nInitialPosY, nNoPT);
aTargetPosZ = linspace(nInitialPosZ, nInitialPosZ+(nNoPT-1)*nSpacing, nNoPT);

nROISize_Hor = 20e-3;
nROISize_Ver = 8e-3;

mROIPosY = zeros(nNoPT, 2);
mROIPosZ = zeros(nNoPT,2);

mROIPosY(:,1) = aTargetPosY - 0.5*nROISize_Hor; % left
mROIPosY(:,2) = aTargetPosY + 0.5*nROISize_Hor; % right

mROIPosZ(:,1) = aTargetPosZ - 0.5*nROISize_Ver; % upper
mROIPosZ(:,2) = aTargetPosZ + 0.5*nROISize_Ver; % bottom

dB = [-6 -6];
interp_ratio = 5;

%% get error images

xcorr_vals_con = {};
xcorr_vals_rsaf = {};
for err_idx = 1:numel(err_list)
    err_tmp = err_list(err_idx).name;
    disp(['>>> ' err_tmp]);
    dir_tmp = [dir_err '/' err_tmp];
    
    spl_list = dir(dir_tmp);
    flag = 0;
    if(strcmp(spl_list(3).name, '.DS_Store')), flag = 1; end
    spl_list = spl_list(3+flag:end);
    
    xcorr_vals_con_tmp = zeros(nNoPT, numel(spl_list));
    xcorr_vals_rsaf_tmp = zeros(nNoPT, numel(spl_list));
    for s_idx = 1:numel(spl_list)
        spl_tmp = spl_list(s_idx).name;
        disp(['    ' spl_tmp]); tic;
        dir_tmpp = [dir_tmp '/' spl_tmp '/Element_64'];
        
        img_list = dir(dir_tmpp);
        
        load([dir_tmpp '/' img_list(end-1).name]);
        img_con = stSaveInfo.mImg_DSC;
        load([dir_tmpp '/' img_list(end).name]);
        img_rsaf = stSaveInfo.mImg_DSC;
        
        for p_idx = 1:nNoPT
            nLft = mROIPosY(p_idx,1);
            nRgt = mROIPosY(p_idx,2);
            nUp = mROIPosZ(p_idx,1);
            nDn = mROIPosZ(p_idx,2);
            
            lIdx = find(abs(aYAxis_DSC-nLft)==min(abs(aYAxis_DSC-nLft)));
            rIdx = find(abs(aYAxis_DSC-nRgt)==min(abs(aYAxis_DSC-nRgt)));
            uIdx = find(abs(aZAxis_DSC-stBFInfo.nRadius-nUp)==min(abs(aZAxis_DSC-stBFInfo.nRadius-nUp)));
            dIdx = find(abs(aZAxis_DSC-stBFInfo.nRadius-nDn)==min(abs(aZAxis_DSC-stBFInfo.nRadius-nDn)));
            
            %%%%%% CON
            roi_con = img_gt_con(uIdx:dIdx, lIdx:rIdx);
            idx_roi = find(roi_con ~= -30);
            roi_con(idx_roi) = roi_con(idx_roi) + abs(max(roi_con(idx_roi)));
            interp_con = interp2(roi_con, interp_ratio);
            idx = find(interp_con > dB(1));
            mask_con = zeros(size(interp_con));
            mask_con(idx) = 1;
            
            roi_con_err = img_con(uIdx:dIdx, lIdx:rIdx);
            idx_roi = find(roi_con_err ~= -30);
            roi_con_err(idx_roi) = roi_con_err(idx_roi) + abs(max(roi_con_err(idx_roi)));
            interp_con_err = interp2(roi_con_err, interp_ratio);
            idx = find(interp_con_err > dB(1));
            mask_con_err = zeros(size(interp_con_err));
            mask_con_err(idx) = 1;
            
            % resize to smaller roi (for lighter computation of normxcorr2)
            first_idx = idx(1);
            end_idx = idx(end);
            mask_center_ver = round(idx(round(0.5*numel(idx)))/size(interp_con,2));
            
            resize_left = max(round(first_idx/size(interp_con,1)) - 100, 1);
            resize_right = min(round(end_idx/size(interp_con,1)) + 100, size(interp_con,2));
            resize_up = max(mask_center_ver - round(0.5 * (resize_right - resize_left)), 1);
            resize_down = min(mask_center_ver + round(0.5 * (resize_right - resize_left)),size(interp_con,1));
            
            mask_con_rs = mask_con(resize_up:resize_down, resize_left:resize_right);
            mask_con_err_rs = mask_con_err(resize_up:resize_down, resize_left:resize_right);
            
            % do xcorr and take the center value
            xcorr_map = normxcorr2(mask_con_rs, mask_con_err_rs);
            center_coord = ceil(size(xcorr_map)/2);
            xcorr_vals_con_tmp(p_idx, s_idx) = xcorr_map(center_coord(1), center_coord(2));
            
            
            %%%%%% rsaf
            roi_rsaf = img_gt_rsaf(uIdx:dIdx, lIdx:rIdx);
            idx_roi = find(roi_con ~= -30);
            roi_rsaf(idx_roi) = roi_rsaf(idx_roi) + abs(max(roi_rsaf(idx_roi)));
            interp_rsaf = interp2(roi_rsaf, interp_ratio);
            idx = find(interp_rsaf > dB(1));
            mask_rsaf = zeros(size(interp_rsaf));
            mask_rsaf(idx) = 1;
            
            roi_rsaf_err = img_rsaf(uIdx:dIdx, lIdx:rIdx);
            idx_roi = find(roi_rsaf_err ~= -30);
            roi_rsaf_err(idx_roi) = roi_rsaf_err(idx_roi) + abs(max(roi_rsaf_err(idx_roi)));
            interp_rsaf_err = interp2(roi_rsaf_err, interp_ratio);
            idx = find(interp_rsaf_err > dB(1));
            mask_rsaf_err = zeros(size(interp_rsaf_err));
            mask_rsaf_err(idx) = 1;
            
            % resize to smaller roi (for lighter computation of normxcorr2)
            first_idx = idx(1);
            end_idx = idx(end);
            mask_center_ver = round(idx(round(0.5*numel(idx)))/size(interp_rsaf,2));
            
            resize_left = round(first_idx/size(interp_con,1)) - 100;
            resize_right = round(end_idx/size(interp_con,1)) + 100;
            resize_up = max(mask_center_ver - round(0.5 * (resize_right - resize_left)), 1);
            resize_down = min(mask_center_ver + round(0.5 * (resize_right - resize_left)),size(interp_con,1));
            
            mask_rsaf_rs = mask_rsaf(resize_up:resize_down, resize_left:resize_right);
            mask_rsaf_err_rs = mask_rsaf_err(resize_up:resize_down, resize_left:resize_right);
            
            % do xcorr and take the center value
            xcorr_map = normxcorr2(mask_rsaf_rs, mask_rsaf_err_rs);
            center_coord = ceil(size(xcorr_map)/2);
            xcorr_vals_rsaf_tmp(p_idx, s_idx) = xcorr_map(center_coord(1), center_coord(2));
        end
        toc;
    end
    xcorr_vals_con = cat(3, xcorr_vals_con, xcorr_vals_con_tmp);
    xcorr_vals_rsaf = cat(3, xcorr_vals_rsaf, xcorr_vals_rsaf_tmp);toc;
end

%% statistics
mean_con = zeros(nNoPT, numel(err_list));
std_con = zeros(nNoPT, numel(err_list));

mean_rsaf = zeros(nNoPT, numel(err_list));
std_rsaf = zeros(nNoPT, numel(err_list));

for err_idx = 1:numel(err_list)
    % CON
    xcorr_con_tmp = cell2mat(xcorr_vals_con(err_idx));
    mean_con(:, err_idx) = mean(xcorr_con_tmp, 2);
    std_con(:,err_idx) = std(xcorr_con_tmp, 0, 2);
    
    % rSAF
    xcorr_rsaf_tmp = cell2mat(xcorr_vals_rsaf(err_idx));
    mean_rsaf(:, err_idx) = mean(xcorr_rsaf_tmp, 2);
    std_rsaf(:,err_idx) = std(xcorr_rsaf_tmp, 0, 2);
end

% pad a case for ground truth (mean = 1, std = 0)
mean_CON = cat(2, ones(nNoPT, 1), mean_con);
mean_rSAF = cat(2, ones(nNoPT, 1), mean_rsaf);
std_CON = cat(2, zeros(nNoPT, 1), std_con);
std_rSAF = cat(2, zeros(nNoPT, 1), std_rsaf);
%% plot
e_axis = [0 0.1 0.2 0.5 1 2 5];

% figure(1);
for d_idx = 1:nNoPT
    %     subplot(nNoPT, 1, d_idx);
    figure(d_idx);
    errorbar(e_axis, mean_CON(d_idx, :), std_CON(d_idx, :), 'LineWidth', 2, 'color', 'k'); hold on;
    errorbar(e_axis, mean_rSAF(d_idx, :), std_rSAF(d_idx, :), 'LineWidth', 2, 'color', 'b'); hold off;
    xlabel('\epsilon'); ylabel('XCORR val'); ylim([0.4 1]);
    legend('CON', 'rSAF');
end

%% take mean along depth
mean_conn = zeros(1, numel(err_list));
std_conn = zeros(1, numel(err_list));

mean_rsaff = zeros(1, numel(err_list));
std_rsaff = zeros(1, numel(err_list));
for err_idx = 1:numel(err_list)
    xcorr_con_tmp = cell2mat(xcorr_vals_con(err_idx));
    mean_conn(err_idx) = mean(xcorr_con_tmp(:));
    std_conn(err_idx) = std(xcorr_con_tmp(:));
    
    xcorr_rsaf_tmp = cell2mat(xcorr_vals_rsaf(err_idx));
    mean_rsaff(err_idx) = mean(xcorr_rsaf_tmp(:));
    std_rsaff(err_idx) = std(xcorr_rsaf_tmp(:));
end

mean_CONN = cat(2, 1, mean_conn);
mean_rSAFF = cat(2, 1, mean_rsaff);
std_CONN = cat(2, 0, std_conn);
std_rSAFF = cat(2, 0, std_rsaff);

figure;
errorbar(e_axis, mean_CONN, std_CONN, 'LineWidth', 4, 'color', 'k'); hold on;
errorbar(e_axis, mean_rSAFF, std_rSAFF, 'LineWidth', 4, 'color', 'b'); hold off;
xlabel('\epsilon'); ylabel('XCORR val'); ylim([0.4 1]);
legend('CON', 'rSAF');

%%
roi__ = roi_rsaf;
figure;
imagesc(aYAxis_DSC(lIdx:rIdx)*1e3,(aZAxis_DSC(uIdx:dIdx)-stBFInfo.nRadius)*1e3,roi__); colormap gray; hold on; caxis([-80 0]);
[ContourLine, h] = contour(aYAxis_DSC(lIdx:rIdx)*1e3,(aZAxis_DSC(uIdx:dIdx)-stBFInfo.nRadius)*1e3,roi__, dB,'ShowText','off', 'LineColor','r', 'LineWidth', 1);
