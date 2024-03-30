%% figure for TMI manuscript
depth_p = 30e-3;
idx_ = find(depth_p == depth_);

%%
dir_save_matfig = '/Users/songhyunwoo/Documents/JohnsHopkins/Research/rSAF/manuscript/MATLAB_fig';
%% 2D psf
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

roi_depths = [3 6];
for p_idx = roi_depths
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
    idx_roi = find(roi_con ~= -30);
    roi_con(idx_roi) = roi_con(idx_roi) + abs(max(roi_con(idx_roi))); % normalize
    
    fig = figure(p_idx*10);
    imagesc(y_roi*1e3, (z_roi-stBFInfo.nRadius)*1e3, roi_con); hold on;
    %     [ContourLine_bg, h] = contour(y_roi*1e3, (z_roi-stBFInfo.nRadius)*1e3, mask_roi, [1 1],'ShowText','off', 'LineColor','r', 'LineWidth', 1);
    caxis([-50 0]); hold off; colormap gray; axis tight; axis equal;
    set(gcf, 'Position', [-2052 1639 371 117]);
    savefig(fig, [dir_save_matfig '/psf_apart_gt_con_' num2str(p_idx) 'cm.fig']);
    
    for e_idx = 1:numel(err_list)
        roi_con = CON{e_idx};
        roi_con = roi_con(uIdx:dIdx, lIdx:rIdx, 1);
        idx_roi = find(roi_con ~= -30);
        roi_con(idx_roi) = roi_con(idx_roi) + abs(max(roi_con(idx_roi))); % normalize
        fig = figure(p_idx*10 + e_idx);
        imagesc(y_roi*1e3, (z_roi-stBFInfo.nRadius)*1e3, roi_con); hold on;
        %     [ContourLine_bg, h] = contour(y_roi*1e3, (z_roi-stBFInfo.nRadius)*1e3, mask_roi, [1 1],'ShowText','off', 'LineColor','r', 'LineWidth', 1);
        caxis([-50 0]); hold off; colormap gray; title(err_list(e_idx).name, 'Interpreter','none');colormap gray; axis tight; axis equal;
        set(gcf, 'Position', [-2052 1639 371 117]);
        savefig(fig, [dir_save_matfig '/psf_apart_con_' num2str(p_idx) 'cm_' err_list(e_idx).name '.fig']);
    end
    
    %%%%%% rsaf
    roi_rsaf = img_rsaf_gt(uIdx:dIdx, lIdx:rIdx);
    idx_roi = find(roi_con ~= -30);
    roi_rsaf(idx_roi) = roi_rsaf(idx_roi) + abs(max(roi_rsaf(idx_roi))); % normalize
    
    fig = figure(p_idx*100);
    imagesc(y_roi*1e3, (z_roi-stBFInfo.nRadius)*1e3, roi_rsaf); hold on;
    %     [ContourLine_bg, h] = contour(y_roi*1e3, (z_roi-stBFInfo.nRadius)*1e3, mask_roi, [1 1],'ShowText','off', 'LineColor','r', 'LineWidth', 1);
    caxis([-50 0]); hold off; colormap gray;colormap gray; axis tight; axis equal;
    set(gcf, 'Position', [-2052 1639 371 117]);
    savefig(fig, [dir_save_matfig '/psf_apart_gt_rsaf_' num2str(p_idx) 'cm.fig']);
    
    for e_idx = 1:numel(err_list)
        roi_rsaf = rSAF{e_idx};
        roi_rsaf = roi_rsaf(uIdx:dIdx, lIdx:rIdx, 1);
        idx_roi = find(roi_rsaf ~= -30);
        roi_rsaf(idx_roi) = roi_rsaf(idx_roi) + abs(max(roi_rsaf(idx_roi))); % normalize
        fig = figure(p_idx*100 + e_idx);
        imagesc(y_roi*1e3, (z_roi-stBFInfo.nRadius)*1e3, roi_rsaf); hold on;
        %     [ContourLine_bg, h] = contour(y_roi*1e3, (z_roi-stBFInfo.nRadius)*1e3, mask_roi, [1 1],'ShowText','off', 'LineColor','r', 'LineWidth', 1);
        caxis([-50 0]); hold off; colormap gray; title(err_list(e_idx).name, 'Interpreter','none');colormap gray; axis tight; axis equal;
        set(gcf, 'Position', [-2052 1639 371 117]);
        savefig(fig, [dir_save_matfig '/psf_apart_rsaf_' num2str(p_idx) 'cm_' err_list(e_idx).name '.fig']);
    end
end

%% 1D psf
% specific depth
d_idx = idx_;

con_ = cell2mat(error_CON(:,:,d_idx));
rsaf_ = cell2mat(error_rSAF(:,:,d_idx));

con_gt_tmp = cell2mat(con_1d_gt(:,:,d_idx));
rsaf_gt_tmp = cell2mat(rsaf_1d_gt(:,:,d_idx));
assert(numel(con_gt_tmp) == numel(rsaf_gt_tmp) , 'Dimension conflicts');

radial_axis = linspace(-0.5*stBFInfo.nFOV, 0.5*stBFInfo.nFOV, numel(con_gt_tmp));

% ground truth
figure(d_idx);
plot(radial_axis, rsaf_gt_tmp, 'LineWidth', 2, 'color', 'k'); hold on; 
figure(d_idx*10);
plot(radial_axis, con_gt_tmp, 'LineWidth', 2, 'color', 'k'); hold on; 
for e_idx = 1:numel(err_list)
    % normalize
    con_norm = con_(:,e_idx);
    con_norm = con_norm(con_norm ~= -30);
    con_norm = con_norm + abs(max(con_norm));
    
    rsaf_norm = rsaf_(:,e_idx);
    rsaf_norm = rsaf_norm(rsaf_norm ~= -30);
    rsaf_norm = rsaf_norm + abs(max(rsaf_norm));
    
    figure(d_idx);
    plot(radial_axis, con_norm, 'LineWidth', 2, 'LineStyle', '-', 'color', [0.15 0.15 0.15]*e_idx);
    
    figure(d_idx*10);
    plot(radial_axis, rsaf_norm, 'LineWidth', 2, 'LineStyle', '-', 'color', [0.15 0.15 0.15]*e_idx);
end

figure(d_idx); hold off;
ylim([-60 0]); xlim([-30 30]); title('CON');
set(gcf, 'Position', [-2935 1200 841 234]);
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

figure(d_idx*10); hold off;
ylim([-60 0]); xlim([-30 30]); title('rSAF');
set(gcf, 'Position', [-2935 806 841 234]);
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
legend('ground truth', '\sigma = 0.1\circ', '\sigma = 0.2\circ', '\sigma = 0.5\circ', ...
    '\sigma = 1\circ', '\sigma = 2\circ', '\sigma = 5\circ');

%%
err_axis = [0.1 0.2 0.5 1 2 5];
for d_idx = (roi_depths-2)
    figure((d_idx+2)*1000); 
    errorbar(err_axis, mean_CON(d_idx, :), std_CON(d_idx, :), 'LineWidth', 2, 'color', 'k'); hold on;
    errorbar(err_axis, mean_rSAF(d_idx, :), std_rSAF(d_idx, :), 'LineWidth', 2, 'color', 'b'); hold off;
    set(gcf, 'Position', [-2914 753 521 426]);
    set(gca, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
end

legend('CON', 'rSAF');