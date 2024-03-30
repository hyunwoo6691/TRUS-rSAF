close all;
%%
max_vals = [max(bp_r5mm(bp_r5mm~=50)) max(bp_r10mm(bp_r10mm~=50)) max(bp_r15mm(bp_r15mm~=50))];

% norm_val = max_vals(1);
norm_val = 1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
depth = 70e-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dor = 2e-3;
%% R = 5mm
bp_tmp = bp_r5mm;

aROI = find(bp_tmp ~= 50);
aOutlier = find(bp_tmp == 50);
mOutput = zeros(size(bp_tmp));
mOutput(aROI) = bp_tmp(aROI);
mOutput_dB = db(mOutput/norm_val);
mOutput_dB(aOutlier) = -150;

figure;
imagesc(y_r5mm*1e3,(z_r5mm-5e-3)*1e3, mOutput_dB); hold on;axis tight; axis equal;
% [ContourLine, h] = contour(axis_y*1e3,(axis_z-nRadius)*1e3,mOutput_dB, [-6 -6],'ShowText','off', 'LineColor','r', 'LineWidth', 1);
xlabel('Elevational [mm]'); ylabel('Axial [mm]');
title(['Elevational aperture size: ' num2str(stTRInfo.nHeight*1e3) 'mm']); hold off;
% caxis([-150 0]); colorbar; colormap default;
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf, 'Position', [378 105 726 535]);


y_axis = y_r5mm;
z_axis = z_r5mm - 5e-3;
% z_axis = z_r5mm;

[z_grid, y_grid] = ndgrid(z_axis, y_axis);

dist_map = sqrt(z_grid.^2 + y_grid.^2);

mask_map_r5mm = logical(zeros(size(dist_map)));
% for d_idx = 1:numel(depth)
%     for k = 1:numel(y_axis)
%         tmp = dist_map(:,k);
%         idx = find(abs(tmp - depth) == min(abs(tmp-depth)));
%         mask_map_r5mm(idx(end), k) = true;
%     end
% end
mask_map_r5mm(dist_map>depth-dor & dist_map<depth+dor) = true;
reject_area = (bp_tmp == 0);

mask_map_r5mm(reject_area) = false;
figure;
imagesc(y_r5mm*1e3, (z_r5mm-5e-3)*1e3,mask_map_r5mm); title('mask map');
axis tight; axis equal;

%% R = 10mm
bp_tmp = bp_r10mm;

aROI = find(bp_tmp ~= 50);
aOutlier = find(bp_tmp == 50);
mOutput = zeros(size(bp_tmp));
mOutput(aROI) = bp_tmp(aROI);
mOutput_dB = db(mOutput/norm_val);
mOutput_dB(aOutlier) = -150;

figure;
imagesc(y_r10mm*1e3,(z_r10mm-10e-3)*1e3, mOutput_dB); hold on;axis tight; axis equal;
% [ContourLine, h] = contour(axis_y*1e3,(axis_z-nRadius)*1e3,mOutput_dB, [-6 -6],'ShowText','off', 'LineColor','r', 'LineWidth', 1);
xlabel('Elevational [mm]'); ylabel('Axial [mm]');
title(['Elevational aperture size: ' num2str(stTRInfo.nHeight*1e3) 'mm']); hold off;
% caxis([-150 0]); colorbar; colormap default;
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf, 'Position', [378 105 726 535]);


y_axis = y_r10mm;
z_axis = z_r10mm - 10e-3;

[z_grid, y_grid] = ndgrid(z_axis, y_axis);

dist_map = sqrt(z_grid.^2 + y_grid.^2);

mask_map_r10mm = logical(zeros(size(dist_map)));
% for d_idx = 1:numel(depth)
%     for k = 1:numel(y_axis)
%         tmp = dist_map(:,k);
%         idx = find(abs(tmp - depth) == min(abs(tmp-depth)));
%         mask_map_r10mm(idx(end), k) = true;
%     end
% end
mask_map_r10mm(dist_map>depth-dor & dist_map<depth+dor) = true;

reject_area = (bp_tmp == 0);
mask_map_r10mm(reject_area) = false;
figure;
imagesc(y_r10mm*1e3, (z_r10mm-10e-3)*1e3,mask_map_r10mm); title('mask map');
axis tight; axis equal;

%% R = 15mm
bp_tmp = bp_r15mm;

aROI = find(bp_tmp ~= 50);
aOutlier = find(bp_tmp == 50);
mOutput = zeros(size(bp_tmp));
mOutput(aROI) = bp_tmp(aROI);
mOutput_dB = db(mOutput/norm_val);
mOutput_dB(aOutlier) = -150;

figure;
imagesc(y_r15mm*1e3,(z_r15mm-15e-3)*1e3, mOutput_dB); hold on;axis tight; axis equal;
% [ContourLine, h] = contour(axis_y*1e3,(axis_z-nRadius)*1e3,mOutput_dB, [-6 -6],'ShowText','off', 'LineColor','r', 'LineWidth', 1);
xlabel('Elevational [mm]'); ylabel('Axial [mm]');
title(['Elevational aperture size: ' num2str(stTRInfo.nHeight*1e3) 'mm']); hold off;
% caxis([-150 0]); colorbar; colormap default;
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf, 'Position', [378 105 726 535]);


y_axis = y_r15mm;
z_axis = z_r15mm - 15e-3;

[z_grid, y_grid] = ndgrid(z_axis, y_axis);

dist_map = sqrt(z_grid.^2 + y_grid.^2);

mask_map_r15mm = logical(zeros(size(dist_map)));
% for d_idx = 1:numel(depth)
%     for k = 1:numel(y_axis)
%         tmp = dist_map(:,k);
%         idx = find(abs(tmp - depth) == min(abs(tmp-depth)));
%         mask_map_r15mm(idx(end), k) = true;
%     end
% end
mask_map_r15mm(dist_map>depth-dor & dist_map<depth+dor) = true;

reject_area = (bp_tmp == 0);
mask_map_r15mm(reject_area) = false;
figure;
imagesc(y_r15mm*1e3, (z_r15mm-15e-3)*1e3,mask_map_r15mm); title('mask map');
axis tight; axis equal;
%%
% bp_1d_r5mm = bp_r5mm(mask_map_r5mm);
% bp_1d_r10mm = bp_r10mm(mask_map_r10mm);
% bp_1d_r15mm = bp_r15mm(mask_map_r15mm);
% 
% axis_r5mm = linspace(-30, 30, numel(bp_1d_r5mm));
% axis_r10mm = linspace(-30, 30, numel(bp_1d_r10mm));
% axis_r15mm = linspace(-30, 30, numel(bp_1d_r15mm));
%%
roi_ref_5 = zeros(size(bp_r5mm));
roi_ref_5(mask_map_r5mm) = bp_r5mm(mask_map_r5mm);
mean_5 = mean(roi_ref_5,1);
mean_5 = mean_5(find(mean_5,1,'first'):find(mean_5,1,'last'));
axis_r5mm = linspace(-28,28, numel(mean_5));


roi_ref_10 = zeros(size(bp_r10mm));
roi_ref_10(mask_map_r10mm) = bp_r10mm(mask_map_r10mm);
mean_10 = mean(roi_ref_10,1);
mean_10 = mean_10(find(mean_10,1,'first'):find(mean_10,1,'last'));
axis_r10mm = linspace(-30, 30, numel(mean_10));

roi_ref_15 = zeros(size(bp_r15mm));
roi_ref_15(mask_map_r15mm) = bp_r15mm(mask_map_r15mm);
mean_15 = mean(roi_ref_15,1);
mean_15 = mean_15(find(mean_15,1,'first'):find(mean_15,1,'last'));
axis_r15mm = linspace(-32, 32, numel(mean_15));

figure;
plot(axis_r5mm, mean_5/max(mean_5), 'LineWidth', 2); hold on;
plot(axis_r10mm, 1*mean_10/max(mean_5), 'LineWidth', 2);
plot(axis_r15mm, 1*mean_15/max(mean_5), 'LineWidth', 2); hold off;
legend('r = 5mm', 'r = 10mm', 'r = 15mm');
xlim([-30 30])
%%
% figure;
% plot(axis_r5mm, bp_1d_r5mm/max(bp_1d_r10mm), 'LineWidth', 2); hold on;
% plot(axis_r10mm, bp_1d_r10mm/max(bp_1d_r10mm), 'LineWidth', 2);
% plot(axis_r15mm, bp_1d_r15mm/max(bp_1d_r10mm), 'LineWidth', 2); hold off;
% legend('r = 5mm', 'r = 10mm', 'r = 15mm');
% figure;
% plot(axis_r5mm, bp_1d_r5mm, 'LineWidth', 2); hold on;
% plot(axis_r10mm, bp_1d_r10mm, 'LineWidth', 2);
% plot(axis_r15mm, bp_1d_r15mm, 'LineWidth', 2); hold off;
% legend('r = 5mm', 'r = 10mm', 'r = 15mm');