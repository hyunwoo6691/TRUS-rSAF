clc; clear; close all;

addpath('functions');
%%
base_angle = 0.4724;
ratio = [1 1.125 1.25 1.375 1.5 1.625 1.75 1.875 2 2.125 2.25 2.375 2.5 2.625 2.75 2.875 3 3.0205];
d_theta_ = fliplr(base_angle./ratio);
%%
dir_ = uigetdir('./data','');

case_list = dir(dir_);
flag = 0;
if(strcmp(case_list(3).name,'.DS_Store')), flag = 1; end
case_list = case_list(3+flag:end);

%%
radius = 15e-3;
depth = 10e-3;
dor = 5e-3; % depth of ROI
%% reference roi
disp('load reference data');
dir_ref = [dir_ '/' case_list(1).name];
[dsc_ref, y_ref, z_ref] = load_bfedData(dir_ref);
disp('complete');

%% load image of all cases
dsc_obj = {};
mask_map_obj = {};
for case_idx = 1:numel(case_list)
    disp(['>>> ' case_list(case_idx).name]);
    dir_obj = [dir_ '/' case_list(case_idx).name];
    
    [dsc_tmp, y_tmp, z_tmp] = load_bfedData(dir_obj);
    mask_map_tmp = getROIMask(y_tmp, z_tmp-radius, depth, dor, find(dsc_tmp == 0),0);
    
    dsc_obj = cat(3, dsc_obj, dsc_tmp);
    mask_map_obj = cat(3, mask_map_obj, mask_map_tmp);
end

%% select ROI
dsc_tmp = dsc_obj{end};
figure(1231);
imagesc(db(dsc_tmp/max(dsc_tmp(:)))); 
axis tight; axis equal; caxis([-80 0]); colormap gray;
set(gcf,'Position',[1196 81 1566 1089]);
roi = drawrectangle();

%% set ROI for reference image
roi_ref = dsc_ref(round(roi.Position(2)):(round(roi.Position(2))+round(roi.Position(4))), ...
                    round(roi.Position(1)):(round(roi.Position(1))+round(roi.Position(3))));
                
figure; imagesc(db(roi_ref/max(dsc_ref(:))));
% mask_map_ref = getROIMask(y_ref, z_ref-radius, depth, dor, find(dsc_ref == 0),1);
% roi_ref = zeros(size(dsc_ref));
% roi_ref(mask_map_ref) = dsc_ref(mask_map_ref);

%% xcorr 2
xcorr_maps = {};
for case_idx = 1:numel(case_list)
    dsc_tmp = dsc_obj{case_idx};
%     mask_map_tmp = mask_map_obj{case_idx};
%     roi_tmp = zeros(size(dsc_tmp));
%     roi_tmp(mask_map_tmp) = dsc_tmp(mask_map_tmp);

    roi_tmp = dsc_tmp(round(roi.Position(2)):(round(roi.Position(2))+round(roi.Position(4))), ...
                        round(roi.Position(1)):(round(roi.Position(1))+round(roi.Position(3))));

    % 2d xcorr
     xcorr_map = normxcorr2(roi_ref/max(roi_ref(:)), roi_tmp/max(roi_ref(:)));
     xcorr_maps = cat(3, xcorr_maps, xcorr_map);
end
%%
xcorr_vals = zeros(1,numel(case_list));
for case_idx = 1:numel(case_list)
    xcorr_map_tmp = xcorr_maps{case_idx};
    [rows_, cols_] = size(xcorr_map_tmp);
    xcorr_vals(case_idx) = xcorr_map_tmp(floor(0.5*(rows_+1)), floor(0.5*(cols_+1)));
end

%%
figure(33);
plot(d_theta_, xcorr_vals,'LineWidth',2,'Marker','^','color','k');

%%
indices = [1 2 3 10 11 12];
for k = indices
    dsc_tmp = dsc_obj{k};
    figure(k);
    imagesc(db(dsc_tmp/max(dsc_tmp(:))));
    title(num2str(d_theta_(k))); 
    colormap gray;
    caxis([-80 0]);
    axis tight; axis equal;
    set(gcf,'Position',[1196 81 1566 1089]);
end
%%
% for k = 1:numel(case_list)
%     dsc_tmp = dsc_obj{k};
%     figure(k);
%     imagesc(db(dsc_tmp/max(dsc_tmp(:))));
%     title(num2str(d_theta_(k))); 
%     colormap gray;
%     caxis([-80 0]);
%     axis tight; axis equal;
%     set(gcf,'Position',[1196 81 1566 1089]);
%     waitforbuttonpress;
% end