% close all;
% depth = 20.26e-3;
depth = 17.5e-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y_axis = axis_y;
% z_axis = aZAxis - nRadius;
z_axis = axis_z;

[z_grid, y_grid] = ndgrid(z_axis, y_axis);

dist_map = sqrt(z_grid.^2 + y_grid.^2);

nFOV_Theta = 132.272;
nRadius = 15e-3;
nDth = 72e-3;

pixel_deg = atand(y_grid./(z_grid));
pixel_deg(pixel_deg > 0.5*nFOV_Theta | pixel_deg < -0.5*nFOV_Theta) = 0;
pixel_deg(sqrt(y_grid.^2 + (z_grid).^2) < nRadius | ...
            sqrt(y_grid.^2 + (z_grid).^2) > (nDth+nRadius)) = 0;
%%
% dist_map = dist_map_2362;
% pixel_deg = pixel_deg_2362;
% load('/Users/songhyunwoo/Desktop/y_2362.mat')
% y_axis = axis_y;

mask_map = logical(zeros(size(dist_map)));
for d_idx = 1:numel(depth)
    for k = 1:numel(y_axis)
        tmp = dist_map(:,k);
        r_tmp = depth(d_idx) + nRadius;
        idx = find(abs(tmp - r_tmp) == min(abs(tmp-r_tmp)));
        mask_map(idx(end), k) = true;
    end
end

reject_area = (pixel_deg == 0);
mask_map(reject_area) = false;
figure;imagesc(axis_y*1e3,(axis_z-nRadius)*1e3,mask_map); title('mask map');
axis tight; axis equal; colorbar;
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf, 'Position', [378 105 726 535]);
%%
figure;
imagesc(axis_y*1e3,(axis_z-nRadius)*1e3, mOutput_dB);
hold on;
[ContourLine, h] = contour(axis_y*1e3,(axis_z - nRadius)*1e3, mask_map, [1 1], 'ShowText', 'on', 'LineWidth', 2, 'color','y');
axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]');hold off;
% title(['Elevational aperture size: ' num2str(stTRInfo.nHeight*1e3) 'mm']); 
caxis([-150 0]); colorbar; colormap default;
% xlim([-10 10]);
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf, 'Position', [378 105 726 535]);


bp_1d = mOutput(mask_map);
% db_bp_1d = db(bp_1d/max(mOutput(:)));
db_bp_1d = db(bp_1d/max(bp_1d));

% radial_axis = linspace(-0.5*nFOV_Theta, 0.5*nFOV_Theta, numel(mask_map(mask_map==true)));
y_tmp = y_grid(mask_map);
z_tmp = z_grid(mask_map);
radial_axis = atand(y_tmp./z_tmp);
%%
hold on;
line([0 0], [-10 100],'color','r');
%%

% gl_angle = 40; % 20mm
% gl_angle = 42.72; % 15mm
gl_angle = 40.96; % 10mm
%%%

figure; 
plot(radial_axis,db_bp_1d); hold on;
line([gl_angle gl_angle],[0 -100],'LineStyle','-.','color','k'); hold off;

xlim([0 60]);
ylim([-60 0]);
