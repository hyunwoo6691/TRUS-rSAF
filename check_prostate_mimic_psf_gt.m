clc; clear; close all;
addpath('functions');
%% set directory
dir_ = uigetdir('./data','');

case_list = dir(dir_);
flag = 0;
if(strcmp('.DS_Store',case_list(3).name)), flag = 1; end
case_list = case_list(3+flag:end);

%%
roi_gt_all = {};
axis_ = {};
radii_ = [];

for caseIdx = 1:numel(case_list)
    case_name = case_list(caseIdx).name;
    dir_tmp = [dir_ '/' case_name];
    
    %% load parameter
    load([dir_tmp '/Parameters.mat']);
    acoustic_ = stParam.stRFInfo;
    bf_ = stParam.stBFInfo;
    trans_ = stParam.stTRInfo;
    
    imgpos_y = stParam.mImgY;
    imgpos_z = stParam.mImgZ;
    
    clear stParam;
    
    % mid processing paramete
    mid_.nTGC_Atten = 0.5;                % [dB]
    
    mid_.nDCRType = 'high';
    mid_.nDCRTap = 128;                   % BPF tap #
    mid_.nDCRFcut = 1e6;
    
    % dsc parameter
    scanline_theta = linspace(-0.5*bf_.nFOV, 0.5*bf_.nFOV, bf_.nScline); % Ground truth transmitted angle
    depth_ = linspace(bf_.nRadius, bf_.nRadius+bf_.nDth, bf_.nDthSpl);
    
    da = abs(scanline_theta(1)-scanline_theta(2));
    dr = abs(depth_(1)-depth_(2));
    view_depth = bf_.nDth + (bf_.nRadius*(1-cosd(0.5*bf_.nFOV)));
    view_width = 2* (bf_.nRadius+bf_.nDth)*sind(0.5*bf_.nFOV);
    
    dz = 1e-4;
    dy = 1e-4;
    height = round(view_depth / dz);
    width = round(view_width / dy);
    
    %% set data directory
    dir_data = [dir_tmp '/errors_bf'];
    
    total_case = dir(dir_data);
    flag = 0;
    if(strcmp(total_case(3).name,'.DS_Store')), flag = 1; end
    total_case = total_case(3+flag:end);
    
    dir_gt = [dir_data '/' total_case(1).name];
    
    %% GROUND TRUTH : show image
    dir_gt_ = [dir_gt '/Sample001/Element_64'];
    bf_cases = dir(dir_gt_);
    
    flag = 0;
    if(strcmp('.DS_Store',bf_cases(3).name)), flag = 1; end
    bf_cases = bf_cases(3+flag:end);
    
    for bfIdx = 1:numel(bf_cases)
        
        disp('>>> load ground truth data');
        load([dir_gt_ '/' bf_cases(bfIdx).name]);
        env_data = mid_proc(stSaveInfo.mBFedData, mid_, acoustic_, bf_);
        [axis_y, axis_z, dsc_gt, reject_idx] = dsc(env_data, dr, da, bf_, height, width, dz, dy);
        
        figure(1);
        imagesc(axis_y*1e3, (axis_z-bf_.nRadius)*1e3, db(dsc_gt/max(dsc_gt(:))));
        axis equal; axis tight; colormap gray; title('GROUND TRUTH');
        caxis([-50 0]);
        set(gcf, 'Position', [561 431 936 516]);
        disp('done'); pause(1);
        
        %% SET ROI POSITION
        depth_pos = [12.5 27.5 42.5 57.5]*1e-3; % point target
        
        % depth_pos = [20 35 50]*1e-3;
        
        dor = 1e-3; % point target
        % dor = 6e-3; % mass target
        %% GROUND TRUTH : extract 1d psf (mean projection)
        disp('>>> extract ground truth data');
        roi_gt = [];
        for d_idx = 1:numel(depth_pos)
            depth_tmp = depth_pos(d_idx);
            
            mask_map = getROIMask(axis_y, axis_z, depth_tmp+bf_.nRadius, bf_.nRadius, dor, reject_idx,1);
            
            roi_tmp = zeros(size(mask_map));
            roi_tmp(mask_map) = dsc_gt(mask_map);
            
            % mean projection
            mean_prj = mean(roi_tmp,1);
            
            roi_gt = cat(1, roi_gt, mean_prj);
        end
        
        axis_ = cat(3, axis_, axis_y);
        roi_gt_all = cat(3, roi_gt_all, roi_gt);
        
        radii_ = cat(2, radii_, bf_.nRadius);
    end
    
    %% GROUND TRUTH : plot 1d psf
    for d_idx = 1:numel(depth_pos)
        roi_tmpp = roi_gt(d_idx,:);
        roi_tmp = roi_tmpp / max(roi_tmpp);
        
        figure(d_idx);
        plot(axis_y*1e3, roi_tmp - mean(roi_tmp), 'LineWidth',2,'color','k');
        
        title(['GROUND TRUTH - Depth: ' num2str(depth_pos(d_idx)*1e3) 'mm']);
        xlim([-15 15]); ylim([0 1]);
        set(gcf,'Position',[150 1400-300*d_idx 400 160]);
    end
    
end
%%
color_ = {'b','k','r'};

for dIdx = 1:numel(depth_pos)
    for cIdx = 1:numel(roi_gt_all)
        roi_gt_tmp = roi_gt_all{cIdx};
        roi_gt_tmpp = roi_gt_tmp(dIdx,:)/max(roi_gt_tmp(dIdx,:));
        figure(dIdx);
%         plot(atand(axis_{cIdx}/(depth_pos(dIdx))), roi_gt_tmpp, 'LineWidth',2,'color',color_{cIdx}); hold on;
        plot(axis_{cIdx}*1e3, roi_gt_tmpp, 'LineWidth',2,'color',color_{cIdx}); hold on;
    end
end
%%
for dIdx = 1:numel(depth_pos)
    figure(dIdx);
    xlim([-20 20]); ylim([0 1]); hold off;
end




