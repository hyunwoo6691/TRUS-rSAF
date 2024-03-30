clc; clear; close all;

addpath('functions');
%%
dir_ = uigetdir('./data','');

case_list = dir(dir_);
flag = 0;
if(strcmp(case_list(3).name,'.DS_Store')), flag = 1; end
case_list = case_list(3+flag:end);

%%
radius = 15e-3;
depth = 10e-3;
dor = 6e-3; % depth of ROI
%% load reference image
disp('load reference data');
dir_ref = [dir_ '/' case_list(1).name];
[dsc_ref, y_ref, z_ref] = load_bfedData(dir_ref);
disp('complete');

reject_area = find(dsc_ref == 0);
%% set roi for reference image
mask_map_ref = getROIMask(y_ref, z_ref, depth+radius, dor, find(dsc_ref == 0), 1);
roi_ref = zeros(size(dsc_ref));
roi_ref(mask_map_ref) = dsc_ref(mask_map_ref);
% Nan to Zero
roi_ref(isnan(roi_ref)) = 0;

figure(120000); 
imagesc(y_ref*1e3, (z_ref-radius)*1e3, db(roi_ref/max(dsc_ref(:))));
axis tight; axis equal; 
% caxis([-80 0]); 
colormap gray;
set(gcf,'Position',[204 678 894 467]);
%% load image of all cases
dsc_obj = {};
mask_map_obj = {};
for case_idx = 2:numel(case_list)
    disp(['>>> ' case_list(case_idx).name]);
    dir_obj = [dir_ '/' case_list(case_idx).name];
    
    [dsc_tmp, y_tmp, z_tmp] = load_bfedData(dir_obj);
    mask_map_tmp = getROIMask(y_tmp, z_tmp, depth+radius, dor, find(dsc_tmp == 0),0);
    
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
MeanPrj_ref = mean(roi_ref,1);
MeanPrj_ref = MeanPrj_ref(find(MeanPrj_ref,1,'first'):find(MeanPrj_ref,1,'last'));
axis_ = linspace(-60, 60, numel(MeanPrj_ref));
figure(123);
plot(axis_,db(MeanPrj_ref/max(MeanPrj_ref)),'LineWidth',2); hold on;

MeanPrj_others = {};
for case_idx = 1:numel(roi_obj)
    roi_tmp = roi_obj{case_idx};
    
    MeanPrj_tmp = mean(roi_tmp,1);
    MeanPrj_tmp = MeanPrj_tmp(find(MeanPrj_tmp,1,'first'):find(MeanPrj_tmp,1,'last'));
    axis_ = linspace(-60, 60, numel(MeanPrj_tmp));
    plot(axis_,db(MeanPrj_tmp/max(MeanPrj_tmp)),'LineWidth',2);
    
    MeanPrj_others = cat(3, MeanPrj_others,MeanPrj_tmp);
end
hold off;
xlim([0 66]);
ylim([-60 0]);
legend('1564','2362','2699','3149','3779','4724');

%% difference
samples_ = zeros(1,numel(roi_obj));
for k = 1:numel(roi_obj)
    samples_(k) = numel(MeanPrj_others{k});
end
num_ = min(samples_);

MeanPrj_ref_norm = db(MeanPrj_ref/max(MeanPrj_ref));
MeanPrj_ref_norm = MeanPrj_ref_norm(1:num_);

figure(12314);
for case_idx = 1:numel(roi_obj)
    MeanPrj_tmp = MeanPrj_others{case_idx};
    MeanPrj_tmp = db(MeanPrj_tmp/max(MeanPrj_tmp));
    
    diff_tmp = MeanPrj_tmp(1:num_)-MeanPrj_ref_norm;
    
    axis_ = linspace(-60, 60, numel(diff_tmp));
    plot(axis_,(diff_tmp),'LineWidth',2); hold on;   
end
hold off;
xlim([0 60]);
legend('2362','2699','3149','3779','4724');