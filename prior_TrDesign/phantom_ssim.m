clc; clear; close all;

element_idx = 63;
mode = {'Conv', 'SA'};

dir_ = uigetdir('./','');
%% load ground truth
for k = 1:2
    load([dir_ '/gt/img_' mode{k} '_' num2str(element_idx) '.mat']);
    if(k == 1)
        con_gt = output.img_mag;
    else
        rsaf_gt = output.img_mag;
    end
end
y_axis = output.y_axis;
z_axis = output.z_axis;

%% load error
err_ = {'0_1', '0_2', '0_5', '1', '2', '5'};
con_err = []; rsaf_err = [];
for e_idx = 1:numel(err_)
    for k = 1:2
        load([dir_ '/err_' err_{e_idx} '/img_' mode{k} '_' num2str(element_idx) '.mat']);
        if(k == 1)
            con_err = cat(3, con_err, output.img_mag);
        else
            rsaf_err = cat(3, rsaf_err, output.img_mag);
        end
    end
end
%% select roi 
figure(1);
imagesc(mag_to_db(rsaf_gt));
colormap gray; caxis([-55 0]);

roi_tmp = drawrectangle();
roi_coord = round(roi_tmp.Position);

%% calculate SSIM
ssim_con = zeros(1, numel(err_));
ssim_rsaf = zeros(1, numel(err_));

edge_sim_con = zeros(1, numel(err_));
edge_sim_rsaf = zeros(1, numel(err_));

for e_idx = 1:numel(err_)
    % ground truth
    roi_con_gt = con_gt(roi_coord(2):roi_coord(2)+roi_coord(4), ...
        roi_coord(1):roi_coord(1)+roi_coord(3));
    
    roi_rsaf_gt = rsaf_gt(roi_coord(2):roi_coord(2)+roi_coord(4), ...
        roi_coord(1):roi_coord(1)+roi_coord(3));
    
    % error
    roi_con_err = con_err(roi_coord(2):roi_coord(2)+roi_coord(4), ...
        roi_coord(1):roi_coord(1)+roi_coord(3),e_idx);
    
    roi_rsaf_err = rsaf_err(roi_coord(2):roi_coord(2)+roi_coord(4), ...
        roi_coord(1):roi_coord(1)+roi_coord(3),e_idx);
    
    
    % ssim
    [ssimval, ssimmap] = ssim(mag_to_db(roi_con_gt), mag_to_db(roi_con_err));
    ssim_con(e_idx) = ssimval;
    [ssimval, ssimmap] = ssim(mag_to_db(roi_rsaf_gt), mag_to_db(roi_rsaf_err));
    ssim_rsaf(e_idx) = ssimval;
    
    % edge
    [edge_rsaf_gt, thresh_rsaf] = edge(roi_rsaf_gt, 'Canny');
    [edge_con_gt, thresh_con] = edge(roi_con_gt, 'Canny', thresh_rsaf);
    edge_con_err = edge(roi_con_err, 'Canny');
    edge_rsaf_err = edge(roi_rsaf_err, 'Canny');
    
    xcorr_map = normxcorr2(edge_con_gt, edge_con_err);
    center_coord = ceil(size(xcorr_map)/2);
    edge_sim_con(e_idx) = xcorr_map(center_coord(1), center_coord(2));
    xcorr_map = normxcorr2(edge_rsaf_gt, edge_rsaf_err);
    center_coord = ceil(size(xcorr_map)/2);
    edge_sim_rsaf(e_idx) = xcorr_map(center_coord(1), center_coord(2));
end

%% plot
err_axis = [0.1 0.2 0.5 1 2 5];
figure(2);
plot(err_axis, ssim_con, 'Marker', 's', 'LineWidth', 2.5, 'color', 'k'); hold on;
plot(err_axis, ssim_rsaf, 'Marker', 'o', 'LineWidth', 2.5, 'color', 'b'); hold off;

grid on; xlabel('\sigma'); ylabel('SSIM');
legend('CON', 'rSAF');

set(gcf, 'Position', [159 300 573 266]);
set(gca, 'FontSize', 14, 'FontWeight', 'bold');

figure(3);imagesc(y_axis(roi_coord(1):roi_coord(1)+roi_coord(3)), ...
                    z_axis(roi_coord(2):roi_coord(2)+roi_coord(4)), ...
                    roi_con_gt); axis tight; axis equal;
set(gcf, 'Position', [159 622 655 254]); title('con');
                
figure(4);imagesc(y_axis(roi_coord(1):roi_coord(1)+roi_coord(3)), ...
                    z_axis(roi_coord(2):roi_coord(2)+roi_coord(4)), ...
                    roi_rsaf_gt); axis tight; axis equal;
set(gcf, 'Position', [800 622 655 254]); title('rsaf')

figure(5);
plot(err_axis, edge_sim_con, 'Marker', 's', 'LineWidth', 2.5, 'color', 'k'); hold on;
plot(err_axis, edge_sim_rsaf, 'Marker', '^', 'LineWidth', 2.5, 'color', 'b'); hold off;

grid on; xlabel('\sigma'); ylabel('edge similarity');
legend('CON', 'rSAF');

set(gcf, 'Position', [800 300 573 266]);
set(gca, 'FontSize', 14, 'FontWeight', 'bold');










