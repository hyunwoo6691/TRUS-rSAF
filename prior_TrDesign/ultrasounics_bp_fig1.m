max_vals = [max(bp_h3mm(bp_h3mm~=50)) max(bp_h4mm(bp_h4mm~=50))...
    max(bp_h5mm(bp_h5mm~=50)) max(bp_h6mm(bp_h6mm~=50)) max(bp_h7mm(bp_h7mm~=50))];

norm_val = max_vals(1);
%%
bp_tmp = bp_h5mm;

aROI = find(bp_tmp ~= 50);
aOutlier = find(bp_tmp == 50);
mOutput = zeros(size(bp_tmp));
mOutput(aROI) = bp_tmp(aROI);
mOutput_dB = db(mOutput/norm_val);
mOutput_dB(aOutlier) = -150;

figure;
imagesc(axis_y*1e3,(axis_z-nRadius)*1e3, mOutput_dB); hold on;axis tight; axis equal;
% [ContourLine, h] = contour(axis_y*1e3,(axis_z-nRadius)*1e3,mOutput_dB, [-6 -6],'ShowText','off', 'LineColor','r', 'LineWidth', 1);
xlabel('Elevational [mm]'); ylabel('Axial [mm]');
title(['Elevational aperture size: ' num2str(stTRInfo.nHeight*1e3) 'mm']); hold off;
caxis([-150 0]); colorbar; colormap default;
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf, 'Position', [378 105 726 535]);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
depth = 70e-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y_axis = axis_y;
z_axis = axis_z - nRadius;

[z_grid, y_grid] = ndgrid(z_axis, y_axis);

dist_map = sqrt(z_grid.^2 + y_grid.^2);

mask_map = logical(zeros(size(dist_map)));
for d_idx = 1:numel(depth)
    for k = 1:numel(y_axis)
        tmp = dist_map(:,k);
        idx = find(abs(tmp - depth) == min(abs(tmp-depth)));
        mask_map(idx(end), k) = true;
    end
end

reject_area = (bp_tmp == 50);
mask_map(reject_area) = false;
figure(3);imagesc(mask_map); title('mask map');

bp_1d_h3mm = bp_h3mm(mask_map);
bp_1d_h4mm = bp_h4mm(mask_map);
bp_1d_h5mm = bp_h5mm(mask_map);
bp_1d_h6mm = bp_h6mm(mask_map);
bp_1d_h7mm = bp_h7mm(mask_map);

axis_ = linspace(-30, 30, numel(bp_1d_h3mm));

figure;
plot(axis_, bp_1d_h3mm/max(bp_1d_h7mm), 'LineWidth', 2); hold on;
plot(axis_, bp_1d_h4mm/max(bp_1d_h7mm), 'LineWidth', 2);
plot(axis_, bp_1d_h5mm/max(bp_1d_h7mm), 'LineWidth', 2);
plot(axis_, bp_1d_h6mm/max(bp_1d_h7mm), 'LineWidth', 2);
plot(axis_, bp_1d_h7mm/max(bp_1d_h7mm), 'LineWidth', 2); hold off;
legend('h = 3mm', 'h = 4mm', 'h = 5mm', 'h = 6mm', 'h = 7mm');

figure;
plot(axis_, bp_1d_h3mm, 'LineWidth', 2); hold on;
plot(axis_, bp_1d_h4mm, 'LineWidth', 2);
plot(axis_, bp_1d_h5mm, 'LineWidth', 2);
plot(axis_, bp_1d_h6mm, 'LineWidth', 2);
plot(axis_, bp_1d_h7mm, 'LineWidth', 2); hold off;
legend('h = 3mm', 'h = 4mm', 'h = 5mm', 'h = 6mm', 'h = 7mm');