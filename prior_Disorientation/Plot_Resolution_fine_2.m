% performance assessment - different focal point
clc; clear; close all;

%% set path
dir_ = uigetdir('./01.Data', '');

dir_data = [dir_ '/Data'];
folder_list = dir(dir_data);
flag = 0;
if(strcmp(folder_list(3).name, '.DS_Store'))
    flag = 1;
end
folder_list = folder_list(3+flag:end);

dir_save = [dir_ '/Figures'];
mkdir(dir_save);

%% get the worst case error at each focal point data

resol_err = {};
for f_idx = 1:numel(folder_list)
    folder_name = folder_list(f_idx).name;
    disp(['>>> Folder name: ' folder_name ]);
    
    dir_tmp = [dir_data '/' folder_name '/Resolutions'];
    error_list = dir(dir_tmp);
    flag = 0;
    if(strcmp(error_list(3).name, '.DS_Store'))
        flag = 1;
    end
    load([dir_tmp '/' error_list(3+flag).name '/Sample001.mat']);
    resol_gt = mLateralResol;
    clear mLateralResol;
    
    error_list = error_list(4+flag:end);
    
    err_rate = zeros(2, size(resol_gt,2), numel(error_list)); % # of beamforming * depth * error case
    
    for e_idx = 1:numel(error_list)
        resol_tmp = zeros(size(resol_gt,1), size(resol_gt,2), 100);
        
        error_tmp = error_list(e_idx).name;
        dir_tmpp = [dir_tmp '/' error_tmp];
        sample_list = dir(dir_tmpp);
        flag = 0;
        if(strcmp(sample_list(3).name, '.DS_Store'))
            flag = 1;
        end
        sample_list = sample_list(3+flag:end);
        
        for s_idx = 1:numel(sample_list)
            sample_tmp = sample_list(s_idx).name;
            load([dir_tmpp '/' sample_tmp]);
            resol_tmp(:,:,s_idx) = mLateralResol;
        end
        
        err_rate_tmp = zeros(2, size(resol_gt,2));
        for bf_idx = 1:2
            gt_ = repmat(resol_gt(bf_idx, :), 100, 1);
            err_ = squeeze(resol_tmp(bf_idx, :, :))';
            
            err_rate_tmpp = abs(err_ - gt_) ./gt_ * 100;
            err_rate_tmp(bf_idx, :) = max(err_rate_tmpp);
%             err_rate_tmp(bf_idx, :) = mean(err_rate_tmpp);
        end
        
        err_rate(:,:,e_idx) = err_rate_tmp;
        
    end
    
    resol_err = cat(3, resol_err, err_rate);
    
end

%% set axis
error_axis = 0.1:0.1:1;
focal_point_axis = 8:1:25;
depth_axis = 1:1:10;

%% plot - for each error case

for bf_idx = 1:2
    
    switch bf_idx
        case 1
            bf_mode = 'Conv';
        case 2
            bf_mode = 'SA';
    end
    
    for e_idx = 1:numel(error_list)
        err_rate_map = zeros(numel(folder_list), size(resol_gt,2));
        for f_idx = 1:numel(folder_list)
            tmp = resol_err{f_idx};
            tmpp = tmp(bf_idx, :, e_idx); % all depth
            err_rate_map(f_idx, :) = tmpp;
        end
        fig_ = figure(e_idx);
        imagesc(depth_axis, focal_point_axis,err_rate_map); colormap jet; colorbar; caxis([0 200]);
        xlabel('depth [cm]'); ylabel('focal point [mm]'); title(['E_max : ' num2str(e_idx/10) 'deg'], 'Interpreter', 'None');
        saveas(fig_, [dir_save '/[' bf_mode ']E_max_' num2str(e_idx/10) '_deg.jpg']);
    end
end

%% plot - for each depth

for bf_idx = 1:2
    
    switch bf_idx
        case 1
            bf_mode = 'Conv';
        case 2
            bf_mode = 'SA';
    end
    
    for d_idx = 1:numel(depth_axis)
        err_rate_map = zeros(numel(folder_list), numel(error_axis));
        for f_idx = 1:numel(folder_list)
            tmp = resol_err{f_idx};
            tmpp = tmp(bf_idx, d_idx, :); % all error case
            err_rate_map(f_idx, :) = tmpp;
        end
        fig_ = figure(d_idx);
        imagesc(error_axis, focal_point_axis,err_rate_map); colormap jet; colorbar; caxis([0 200]);
        xlabel('E_max [deg]','Interpreter', 'None'); ylabel('focal point [mm]'); title(['Depth : ' num2str(d_idx) 'cm']);
        saveas(fig_, [dir_save '/[' bf_mode ']Depth_' num2str(d_idx) 'cm.jpg']);
    end
end

close all;