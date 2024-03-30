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
%% reference roi
disp('load reference data');
dir_ref = [dir_ '/' case_list(1).name];
[dsc_ref, y_ref, z_ref] = load_bfedData(dir_ref);
disp('complete');

%% load image of all cases
dsc_obj = {};
for case_idx = 1:numel(case_list)
    disp(['>>> ' case_list(case_idx).name]);
    dir_obj = [dir_ '/' case_list(case_idx).name];
    
    [dsc_tmp, y_tmp, z_tmp] = load_bfedData(dir_obj);
    dsc_obj = cat(3, dsc_obj, dsc_tmp);
end

%%
dsc_obj_ = zeros(809,1591);

dsc_tmp = dsc_obj{1};
dsc_obj_(:,:,1) = dsc_tmp;

dsc_tmp = dsc_obj{2};
dsc_obj_(:,:,2) = dsc_tmp;

dsc_tmp = dsc_obj{3};
dsc_obj_(:,1:1590,3) = dsc_tmp;

dsc_tmp = dsc_obj{4};
dsc_obj_(:,2:1590,4) = dsc_tmp;

dsc_tmp = dsc_obj{5};
dsc_obj_(:,2:1590,5) = dsc_tmp;

dsc_tmp = dsc_obj{6};
dsc_obj_(:,:,6) = dsc_tmp;


%% diff
aOutlier = find(dsc_ref == 0);

dsc_ref_ = dsc_obj_(:,:,1)/max(max(dsc_obj_(:,:,1)));
for case_idx = 1:numel(case_list)
    dsc_tmp = dsc_obj_(:,:,case_idx)/max(max(dsc_obj_(:,:,case_idx)));
%     dsc_tmp = dsc_obj_(:,:,case_idx)/max(max(dsc_obj_(:,:,1)));
    img_tmp = dsc_tmp - dsc_ref_;
    img_tmp(aOutlier) = 1e-1;
%     img_tmp(isnan(img_tmp)) = 1e-1;
    figure(case_idx);
%     imagesc(db(img_tmp/max(img_tmp(:))));
    imagesc(y_ref*1e3, (z_ref-radius)*1e3,db(img_tmp));
    axis tight; axis equal;
    title(case_list(case_idx).name,'Interpreter','none');
    xlim([0 y_ref(end)]*1e3);
    colormap gray;
    caxis([-100 0]);
    
    disp(max(max(dsc_obj_(:,:,case_idx))));
end