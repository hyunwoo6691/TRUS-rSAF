clear; close all; clc;

%% set directory
dir_ = uigetdir('', '');

dir_data = [dir_ '/errors'];
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
colormap gray; colorbar; xlabel('Elevational'); ylabel('Axial'); title('CON'); caxis([-50 0]);
subplot(1,2,2);imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, img_rsaf_gt);
colormap gray; colorbar; xlabel('Elevational'); ylabel('Axial'); title('rSAF'); caxis([-50 0]);
set(gcf, 'Position', [-2933 1039 1533 556]);

figure(2);
imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, img_con_gt); hold on;
line([line_one(1,1) line_one(2,1)], [line_one(1,2) line_one(2,2)], 'color', 'r', 'LineWidth', 1);
line([line_two(1,1) line_two(2,1)], [line_two(1,2) line_two(2,2)], 'color', 'r', 'LineWidth', 1); hold off;
colormap gray; colorbar; xlabel('Elevational'); ylabel('Axial'); caxis([-50 0]);
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
nROISize_Ver = 8e-3;

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
    
    %%%%%% CON
    roi_con = img_con_gt(uIdx:dIdx, lIdx:rIdx);
    idx_roi = find(roi_con ~= -30);
    roi_con(idx_roi) = roi_con(idx_roi) + abs(max(roi_con(idx_roi))); % normalize
    
    %%%%%% rsaf
    roi_rsaf = img_rsaf_gt(uIdx:dIdx, lIdx:rIdx);
    idx_roi = find(roi_con ~= -30);
    roi_rsaf(idx_roi) = roi_rsaf(idx_roi) + abs(max(roi_rsaf(idx_roi))); % normalize
    
    figure(5);
    subplot(5, 2, (p_idx-1-2)*2+1);
    imagesc(roi_con); title(['Depth: ' num2str(p_idx) 'cm | CON']); colormap default; colorbar; caxis([-50 0]);
    subplot(5, 2, (p_idx-2)*2);
    imagesc(roi_rsaf); title(['Depth: ' num2str(p_idx) 'cm | rSAF']); colormap default; colorbar; caxis([-50 0]);
end
set(gcf, 'Position', [-2393 1060 1381 797]);

%% compare image at each depth
% CON
for e_idx = 1:numel(err_list)
    con_err = cell2mat(CON(e_idx));
    rand_idx = round((rand(1)+0.001)*numel(sample_list));
    con_err_tmp = con_err(:,:,rand_idx);
    for p_idx = 3:nNoPT
        nLft = mROIPosY(p_idx,1);
        nRgt = mROIPosY(p_idx,2);
        nUp = mROIPosZ(p_idx,1);
        nDn = mROIPosZ(p_idx,2);
        
        lIdx = find(abs(aYAxis_DSC-nLft)==min(abs(aYAxis_DSC-nLft)));
        rIdx = find(abs(aYAxis_DSC-nRgt)==min(abs(aYAxis_DSC-nRgt)));
        uIdx = find(abs(aZAxis_DSC-stBFInfo.nRadius-nUp)==min(abs(aZAxis_DSC-stBFInfo.nRadius-nUp)));
        dIdx = find(abs(aZAxis_DSC-stBFInfo.nRadius-nDn)==min(abs(aZAxis_DSC-stBFInfo.nRadius-nDn)));
        
        %%%%%% ground truth
        roi_con_gt = img_con_gt(uIdx:dIdx, lIdx:rIdx);
        idx_roi = find(roi_con_gt ~= -30);
        roi_con_gt(idx_roi) = roi_con_gt(idx_roi) + abs(max(roi_con_gt(idx_roi))); % normalize
        
        %%%%%% error - randomly sample one image
        roi_con_err = con_err_tmp(uIdx:dIdx, lIdx:rIdx);
        idx_roi = find(roi_con_err ~= -30);
        roi_con_err(idx_roi) = roi_con_err(idx_roi) + abs(max(roi_con_err(idx_roi))); % normalize
        
        figure(e_idx*10); sgtitle(['CON | ' err_list(e_idx).name], 'Interpreter', 'None');
        subplot(5, 2, (p_idx-1-2)*2+1);
        imagesc(roi_con_gt); title(['Depth: ' num2str(p_idx) 'cm | ground truth']); colormap default; caxis([-40 0]);
        subplot(5, 2, (p_idx-2)*2);
        imagesc(roi_con_err); title(['Depth: ' num2str(p_idx) 'cm | error']); colormap default; caxis([-40 0]);
    end
    set(gcf, 'Position', [-2393 1060 1381 797]);
end

% rSAF
for e_idx = 1:numel(err_list)
    rSAF_err = cell2mat(rSAF(e_idx));
    rand_idx = round(rand(1)*numel(sample_list));
    rSAF_err_tmp = rSAF_err(:,:,rand_idx);
    for p_idx = 3:nNoPT
        nLft = mROIPosY(p_idx,1);
        nRgt = mROIPosY(p_idx,2);
        nUp = mROIPosZ(p_idx,1);
        nDn = mROIPosZ(p_idx,2);
        
        lIdx = find(abs(aYAxis_DSC-nLft)==min(abs(aYAxis_DSC-nLft)));
        rIdx = find(abs(aYAxis_DSC-nRgt)==min(abs(aYAxis_DSC-nRgt)));
        uIdx = find(abs(aZAxis_DSC-stBFInfo.nRadius-nUp)==min(abs(aZAxis_DSC-stBFInfo.nRadius-nUp)));
        dIdx = find(abs(aZAxis_DSC-stBFInfo.nRadius-nDn)==min(abs(aZAxis_DSC-stBFInfo.nRadius-nDn)));
        
        %%%%%% ground truth
        roi_rSAF_gt = img_rsaf_gt(uIdx:dIdx, lIdx:rIdx);
        idx_roi = find(roi_rSAF_gt ~= -30);
        roi_rSAF_gt(idx_roi) = roi_rSAF_gt(idx_roi) + abs(max(roi_rSAF_gt(idx_roi))); % normalize
        
        %%%%%% error - randomly sample one image
        roi_rSAF_err = rSAF_err_tmp(uIdx:dIdx, lIdx:rIdx);
        idx_roi = find(roi_rSAF_err ~= -30);
        roi_rSAF_err(idx_roi) = roi_rSAF_err(idx_roi) + abs(max(roi_rSAF_err(idx_roi))); % normalize
        
        figure(e_idx*100); sgtitle(['rSAF | ' err_list(e_idx).name], 'Interpreter', 'None');
        subplot(5, 2, (p_idx-1-2)*2+1);
        imagesc(roi_rSAF_gt); title(['Depth: ' num2str(p_idx) 'cm | ground truth']); colormap default; caxis([-40 0]);
        subplot(5, 2, (p_idx-2)*2);
        imagesc(roi_rSAF_err); title(['Depth: ' num2str(p_idx) 'cm | error']); colormap default; caxis([-40 0]);
    end
    set(gcf, 'Position', [-2393 800 1381 797]);
end

%% SSIM at each depth

ssim_CON = zeros(numel(nNoPT)-2, numel(err_list), numel(sample_list));
ssim_rSAF = zeros(numel(nNoPT)-2, numel(err_list), numel(sample_list));
for e_idx = 1:numel(err_list)
    disp(['>>> Error: ' err_list(e_idx).name ]);
    rSAF_err = cell2mat(rSAF(e_idx));
    con_err = cell2mat(CON(e_idx));
    
    for s_idx = 1:numel(sample_list)
        rSAF_err_tmp = rSAF_err(:,:,s_idx);
        con_err_tmp = con_err(:,:,s_idx);
        
        for p_idx = 3:nNoPT
            nLft = mROIPosY(p_idx,1);
            nRgt = mROIPosY(p_idx,2);
            nUp = mROIPosZ(p_idx,1);
            nDn = mROIPosZ(p_idx,2);
            
            lIdx = find(abs(aYAxis_DSC-nLft)==min(abs(aYAxis_DSC-nLft)));
            rIdx = find(abs(aYAxis_DSC-nRgt)==min(abs(aYAxis_DSC-nRgt)));
            uIdx = find(abs(aZAxis_DSC-stBFInfo.nRadius-nUp)==min(abs(aZAxis_DSC-stBFInfo.nRadius-nUp)));
            dIdx = find(abs(aZAxis_DSC-stBFInfo.nRadius-nDn)==min(abs(aZAxis_DSC-stBFInfo.nRadius-nDn)));
            
            %%%%%% ground truth
            % CON
            roi_CON_gt = img_con_gt(uIdx:dIdx, lIdx:rIdx);
            idx_roi = find(roi_CON_gt ~= -30);
            roi_CON_gt(idx_roi) = roi_CON_gt(idx_roi) + abs(max(roi_CON_gt(idx_roi))); % normalize
            % rSAF
            roi_rSAF_gt = img_rsaf_gt(uIdx:dIdx, lIdx:rIdx);
            idx_roi = find(roi_rSAF_gt ~= -30);
            roi_rSAF_gt(idx_roi) = roi_rSAF_gt(idx_roi) + abs(max(roi_rSAF_gt(idx_roi))); % normalize
            
            %%%%%% error - randomly sample one image
            % CON
            roi_con_err = con_err_tmp(uIdx:dIdx, lIdx:rIdx);
            idx_roi = find(roi_con_err ~= -30);
            roi_con_err(idx_roi) = roi_con_err(idx_roi) + abs(max(roi_con_err(idx_roi))); % normalize
            % rSAF
            roi_rSAF_err = rSAF_err_tmp(uIdx:dIdx, lIdx:rIdx);
            idx_roi = find(roi_rSAF_err ~= -30);
            roi_rSAF_err(idx_roi) = roi_rSAF_err(idx_roi) + abs(max(roi_rSAF_err(idx_roi))); % normalize
            
            %%%%%% SSIM calculation
            % CON
            [ssimval, ssimmap] = ssim(roi_CON_gt, roi_con_err);
            ssim_CON(p_idx-2, e_idx, s_idx) = ssimval;
            % rSAF
            [ssimval, ssimmap] = ssim(roi_rSAF_gt, roi_rSAF_err);
            ssim_rSAF(p_idx-2, e_idx, s_idx) = ssimval;
        end
    end
end
%% get mean & std at each depth
err_axis = [0.1 0.2 0.5 1 2 5];

ssim_CON_mean = mean(ssim_CON, 3);
ssim_rSAF_mean = mean(ssim_rSAF, 3);

ssim_CON_std = size(ssim_CON_mean);
ssim_rSAF_std = size(ssim_rSAF_mean);
for d_idx = 1:size(ssim_CON,1)
    for e_idx = 1:size(ssim_CON,2)
        ssim_CON_std(d_idx, e_idx) = std(squeeze(ssim_CON(d_idx, e_idx, :)));
        ssim_rSAF_std(d_idx, e_idx) = std(squeeze(ssim_rSAF(d_idx, e_idx, :)));
    end
end

for d_idx = 1:size(ssim_CON_mean, 1)
    figure(10000); subplot(1, size(ssim_CON_mean,1), d_idx);
    errorbar(err_axis, ssim_CON_mean(d_idx, :), ssim_CON_std(d_idx, :), 'LineWidth', 2, 'color', 'k'); hold on;
    errorbar(err_axis, ssim_rSAF_mean(d_idx, :), ssim_rSAF_std(d_idx, :), 'LineWidth', 2, 'color', 'b'); hold off;
    title(['Depth: ' num2str(d_idx+2) 'cm']); xlabel('/varepsilon [deg]'); ylabel('SSIM'); grid on;
    legend('CON', 'rSAF'); ylim([min(min(ssim_CON_mean(:)-ssim_CON_std(:)), min(ssim_rSAF_mean(:)-ssim_rSAF_std(:))) 1]);
end
set(gcf, 'Position', [-2807 953 1243 510]);

%% get overall mean && std
ssim_CON_overall_mean = zeros(1, numel(err_list)); ssim_CON_overall_std = zeros(1, numel(err_list));
ssim_rSAF_overall_mean = zeros(1, numel(err_list)); ssim_rSAF_overall_std = zeros(1, numel(err_list));

for e_idx = 1:numel(err_list)
    ssim_con_tmp = squeeze(ssim_CON(:,e_idx,:));
    ssim_CON_overall_mean(e_idx) = mean(ssim_con_tmp(:));
    ssim_CON_overall_std(e_idx) = std(ssim_con_tmp(:));
    
    ssim_rsaf_tmp = squeeze(ssim_rSAF(:,e_idx,:));
    ssim_rSAF_overall_mean(e_idx) = mean(ssim_rsaf_tmp(:));
    ssim_rSAF_overall_std(e_idx) = std(ssim_rsaf_tmp(:));
end

figure(10001);
errorbar(err_axis, ssim_CON_overall_mean, ssim_CON_overall_std, 'LineWidth', 2, 'color', 'k'); hold on;
errorbar(err_axis, ssim_rSAF_overall_mean, ssim_rSAF_overall_std, 'LineWidth', 2, 'color', 'b'); hold off;
xlabel('/varepsilon [deg]'); ylabel('SSIM'); grid on;
legend('CON', 'rSAF');


%% compute SSIM for entire image
ssim_CON_entire = zeros(numel(err_list), numel(sample_list));
ssim_rSAF_entire = zeros(numel(err_list), numel(sample_list));

xcorr_vals_rsaf_tmp = zeros(numel(err_list), numel(sample_list));
xcorr_vals_con_tmp = zeros(numel(err_list), numel(sample_list));

for e_idx = 1:numel(err_list)
    disp(['>>> Error: ' err_list(e_idx).name ]);
    rSAF_err = cell2mat(rSAF(e_idx));
    con_err = cell2mat(CON(e_idx));
    
    for s_idx = 1:numel(sample_list)
        rSAF_err_tmp = rSAF_err(:,:,s_idx);
        con_err_tmp = con_err(:,:,s_idx);
        
        %%%%%% ground truth
        % CON
        roi_CON_gt = img_con_gt;
        idx_roi = find(roi_CON_gt ~= -30);
        roi_CON_gt(idx_roi) = roi_CON_gt(idx_roi) + abs(max(roi_CON_gt(idx_roi))); % normalize
        % rSAF
        roi_rSAF_gt = img_rsaf_gt;
        idx_roi = find(roi_rSAF_gt ~= -30);
        roi_rSAF_gt(idx_roi) = roi_rSAF_gt(idx_roi) + abs(max(roi_rSAF_gt(idx_roi))); % normalize
        
        %%%%%% error - randomly sample one image
        % CON
        roi_con_err = con_err_tmp;
        idx_roi = find(roi_con_err ~= -30);
        roi_con_err(idx_roi) = roi_con_err(idx_roi) + abs(max(roi_con_err(idx_roi))); % normalize
        % rSAF
        roi_rSAF_err = rSAF_err_tmp;
        idx_roi = find(roi_rSAF_err ~= -30);
        roi_rSAF_err(idx_roi) = roi_rSAF_err(idx_roi) + abs(max(roi_rSAF_err(idx_roi))); % normalize
        
        %%%%%% SSIM calculation
        % CON
        [ssimval, ssimmap] = ssim(roi_CON_gt(roi_CON_gt ~= -30), roi_con_err(roi_con_err ~= -30));
        ssim_CON_entire(e_idx, s_idx) = ssimval;
        % rSAF
        [ssimval, ssimmap] = ssim(roi_rSAF_gt(roi_rSAF_gt ~= -30), roi_rSAF_err(roi_rSAF_err ~= -30));
        ssim_rSAF_entire(e_idx, s_idx) = ssimval;
        
%         %%%%%% XCORR
%         xcorr_map = normxcorr2(roi_CON_gt, roi_con_err);
%         center_coord = ceil(size(xcorr_map)/2);
%         xcorr_vals_con_tmp(e_idx, s_idx) = xcorr_map(center_coord(1), center_coord(2));
%         xcorr_map = normxcorr2(roi_rSAF_gt, roi_rSAF_err);
%         center_coord = ceil(size(xcorr_map)/2);
%         xcorr_vals_rsaf_tmp(e_idx, s_idx) = xcorr_map(center_coord(1), center_coord(2));
    end
end


%% get mean & std
ssim_CON_mean_entire = mean(ssim_CON_entire, 2);
ssim_CON_std_entire = std(ssim_CON_entire, 0, 2);

ssim_rSAF_mean_entire = mean(ssim_rSAF_entire, 2);
ssim_rSAF_std_entire = std(ssim_rSAF_entire, 0, 2);

figure(10001);
errorbar(err_axis, ssim_CON_mean_entire, ssim_CON_std_entire, 'LineWidth', 2, 'color', 'k'); hold on;
errorbar(err_axis, ssim_rSAF_mean_entire, ssim_rSAF_std_entire, 'LineWidth', 2, 'color', 'b'); hold off;
xlabel('\varepsilon [deg]'); ylabel('SSIM'); grid on;
legend('CON', 'rSAF');

xcorr_CON_mean_entire = mean(xcorr_vals_con_tmp, 2);
xcorr_CON_std_entire = std(xcorr_vals_con_tmp, 0, 2);

xcorr_rSAF_mean_entire = mean(xcorr_vals_rsaf_tmp, 2);
xcorr_rSAF_std_entire = std(xcorr_vals_rsaf_tmp, 0, 2);

figure(10002);
errorbar(err_axis, xcorr_CON_mean_entire, xcorr_CON_std_entire, 'LineWidth', 2, 'color', 'k'); hold on;
errorbar(err_axis, xcorr_rSAF_mean_entire, xcorr_rSAF_std_entire, 'LineWidth', 2, 'color', 'b'); hold off;
xlabel('\varepsilon [deg]'); ylabel('SSIM'); grid on;
legend('CON', 'rSAF');
