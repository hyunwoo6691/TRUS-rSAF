%%
adB_Overlay = [-6 -6];
radius = 10e-3;

xlim_ = [-30 30];
caxis_ = -60;
%%
mOutput = bp_h3mm;

aOutlier = find(mOutput == 0);
mOutput_dB = db(mOutput/max(mOutput(:)));
mOutput_dB(aOutlier) = -150;

figure(1);
imagesc(y_axis*1e3,(z_axis-radius)*1e3, mOutput_dB); hold on;
[ContourLine, h] = contour(y_axis*1e3,(z_axis - radius)*1e3, mOutput_dB, adB_Overlay, 'ShowText', 'off', 'LineWidth', 2, 'color','r');
axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]'); 
title(['Elevational aperture size: 3mm']); hold off;
caxis([caxis_ 0]); colorbar; colormap default;
xlim([-30 30]);
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf, 'Position', [378 105 726 535]);

% h = 4mm
mOutput = bp_h4mm;

aOutlier = find(mOutput == 0);
mOutput_dB = db(mOutput/max(mOutput(:)));
mOutput_dB(aOutlier) = -150;

figure(2);
imagesc(y_axis*1e3,(z_axis-radius)*1e3, mOutput_dB); hold on;
[ContourLine, h] = contour(y_axis*1e3,(z_axis - radius)*1e3, mOutput_dB, adB_Overlay, 'ShowText', 'off', 'LineWidth', 2, 'color','r');
axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]'); 
title(['Elevational aperture size: 4mm']); hold off;
caxis([caxis_ 0]); colorbar; colormap default;
xlim([xlim_(1) xlim_(2)]);
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf, 'Position', [378 105 726 535]);

% h = 5mm
mOutput = bp_h5mm;

aOutlier = find(mOutput == 0);
mOutput_dB = db(mOutput/max(mOutput(:)));
mOutput_dB(aOutlier) = -150;

figure(3);
imagesc(y_axis*1e3,(z_axis-radius)*1e3, mOutput_dB); hold on;
[ContourLine, h] = contour(y_axis*1e3,(z_axis - radius)*1e3, mOutput_dB, adB_Overlay, 'ShowText', 'off', 'LineWidth', 2, 'color','r');
axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]'); 
title(['Elevational aperture size: 5mm']); hold off;
caxis([caxis_ 0]); colorbar; colormap default;
xlim([xlim_(1) xlim_(2)]);
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf, 'Position', [378 105 726 535]);

% h = 6mm
mOutput = bp_h6mm;

aOutlier = find(mOutput == 0);
mOutput_dB = db(mOutput/max(mOutput(:)));
mOutput_dB(aOutlier) = -150;

figure(4);
imagesc(y_axis*1e3,(z_axis-radius)*1e3, mOutput_dB); hold on;
[ContourLine, h] = contour(y_axis*1e3,(z_axis - radius)*1e3, mOutput_dB, adB_Overlay, 'ShowText', 'off', 'LineWidth', 2, 'color','r');
axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]'); 
title(['Elevational aperture size: 6mm']); hold off;
caxis([caxis_ 0]); colorbar; colormap default;
xlim([xlim_(1) xlim_(2)]);
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf, 'Position', [378 105 726 535]);

% h = 7mm
mOutput = bp_h7mm;

aOutlier = find(mOutput == 0);
mOutput_dB = db(mOutput/max(mOutput(:)));
mOutput_dB(aOutlier) = -150;

figure(5);
imagesc(y_axis*1e3,(z_axis-radius)*1e3, mOutput_dB); hold on;
[ContourLine, h] = contour(y_axis*1e3,(z_axis - radius)*1e3, mOutput_dB, adB_Overlay, 'ShowText', 'off', 'LineWidth', 2, 'color','r');
axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]'); 
title(['Elevational aperture size: 7mm']); hold off;
caxis([caxis_ 0]); colorbar; colormap default;
xlim([xlim_(1) xlim_(2)]);
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf, 'Position', [378 105 726 535]);