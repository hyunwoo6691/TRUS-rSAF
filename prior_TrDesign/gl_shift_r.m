clc; clear; close all;

addpath('functions');
%%
dir_ = uigetdir('./data','');

case_list = dir(dir_);
flag = 0;
if(strcmp(case_list(3).name,'.DS_Store')), flag = 1; end
case_list = case_list(3+flag:end);

%%
radius = [5e-3 10e-3 15e-3];
depth = 70e-3;
dor = 1e-3; % depth of ROI

%% load image of all cases
dsc_obj = {};
mask_map_obj = {};
for case_idx = 1:numel(case_list)
    disp(['>>> ' case_list(case_idx).name]);
    dir_obj = [dir_ '/' case_list(case_idx).name];
    
    [dsc_tmp, y_tmp, z_tmp] = load_bfedData(dir_obj);
    mask_map_tmp = getROIMask(y_tmp, z_tmp, depth+radius(case_idx), dor, find(dsc_tmp == 0),0);
    
    dsc_obj = cat(3, dsc_obj, dsc_tmp);
    mask_map_obj = cat(3, mask_map_obj, mask_map_tmp);
end
disp('    done');

%% get roi for obj images
roi_obj = {};
for case_idx = 1:numel(dsc_obj)
    dsc_tmp = dsc_obj{case_idx};
    roi_tmp = zeros(size(mask_map_obj{case_idx}));
    roi_tmp(mask_map_obj{case_idx}) = dsc_tmp(mask_map_obj{case_idx});
    
    % Nan to Zero
    roi_tmp(isnan(roi_tmp)) = 0;
    
    roi_obj = cat(3, roi_obj, roi_tmp);
end

%% mean projection
% fov_ = [55 65 63];
fov_ = [63 63 63];
figure(123);
MeanPrj_others = {};
for case_idx = 1:numel(roi_obj)
    roi_tmp = roi_obj{case_idx};
    
    MeanPrj_tmp = mean(roi_tmp,1);
    MeanPrj_tmp = MeanPrj_tmp(find(MeanPrj_tmp,1,'first'):find(MeanPrj_tmp,1,'last'));
    axis_ = linspace(-fov_(case_idx), fov_(case_idx), numel(MeanPrj_tmp));
    plot(axis_,db(MeanPrj_tmp/max(MeanPrj_tmp)),'LineWidth',2); hold on;
    
    MeanPrj_others = cat(3, MeanPrj_others,MeanPrj_tmp);
end
hold off;
xlim([0 66]);
ylim([-60 0]);
legend('r=5mm','r=10mm','r=15mm');
