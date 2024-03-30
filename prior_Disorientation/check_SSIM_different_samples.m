clear; close all; clc;

%% set directory
dir_ = uigetdir('', '');

dir_data = [dir_ '/Errors'];
dir_save = [dir_ '/Saved']; mkdir(dir_save);

dir_gt = [dir_data '/error_0/Sample001/Element_64'];

gt_list = dir(dir_gt);

%% load parameter
load([dir_ '/Parameters.mat']);
stRFInfo = stParam.stRFInfo;
stBFInfo = stParam.stBFInfo;
stTRInfo = stParam.stTRInfo;
stTxInfo = stParam.stTxInfo;
mImgY = stParam.mImgY;
mImgZ = stParam.mImgZ;
aYAxis_DSC = stParam.aYAxis_DSC;
aZAxis_DSC = stParam.aZAxis_DSC;

clear stParam

%% load ground truth
name_rsaf = gt_list(end).name;
name_con = gt_list(end-1).name;

load([dir_gt '/' name_rsaf]);
img_rsaf_gt = stSaveInfo.mImg_DSC;

load([dir_gt '/' name_con]);
img_con_gt = stSaveInfo.mImg_DSC;

clear stSaveInfo;
%% Image
% deg = 3.75;
% line_one = [-stBFInfo.nRadius*sind(deg), -stBFInfo.nRadius*(1-cosd(deg));...
%     -(stBFInfo.nRadius+stBFInfo.nDth)*sind(deg), (stBFInfo.nRadius + stBFInfo.nDth)*cosd(deg)-stBFInfo.nRadius]*1e3; % y, z

% line_two = [stBFInfo.nRadius*sind(deg), -stBFInfo.nRadius*(1-cosd(deg));...
%     (stBFInfo.nRadius+stBFInfo.nDth)*sind(deg), (stBFInfo.nRadius + stBFInfo.nDth)*cosd(deg)-stBFInfo.nRadius]*1e3; % y, z

figure(10);
subplot(1,2,1);imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, img_con_gt);
colormap gray; colorbar; xlabel('Elevational'); ylabel('Axial'); title('CON'); caxis([-50 0]);
subplot(1,2,2);imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, img_rsaf_gt);
colormap gray; colorbar; xlabel('Elevational'); ylabel('Axial'); title('rSAF'); caxis([-50 0]);
set(gcf, 'Position', [-2933 1039 1533 556]);

figure(20);
imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, img_con_gt); hold on;
% line([line_one(1,1) line_one(2,1)], [line_one(1,2) line_one(2,2)], 'color', 'r', 'LineWidth', 1);
% line([line_two(1,1) line_two(2,1)], [line_two(1,2) line_two(2,2)], 'color', 'r', 'LineWidth', 1); hold off;
colormap gray; colorbar; xlabel('Elevational'); ylabel('Axial'); caxis([-50 0]);
set(gcf, 'Position', [-2287 856 788 640]);
%% load data with error
err_list = dir(dir_data);
flag = 0; if(strcmp(err_list(3).name, '.DS_Store')), flag = 1; end
err_list = err_list(4+flag:end);

disp('>>> Data loading....');
CON = {}; rSAF = {};
for e_idx = 1:numel(err_list)
    err_tmp = err_list(e_idx).name;
    
    disp(['    Error: ' err_tmp]);
    
    dir_tmp = [dir_data '/' err_tmp];
    
    sample_list = dir(dir_tmp);
    flag = 0; if(strcmp(sample_list(3).name, '.DS_Store')), flag = 1; end
    sample_list = sample_list(3+flag:end);
    
    con_ = []; rsaf_ = [];
    for s_idx = 1:numel(sample_list)
        sample_tmp = sample_list(s_idx).name;
        dir_tmpp = [dir_tmp '/' sample_tmp '/Element_64'];
        
        bf_cases = dir(dir_tmpp);
        con_tmp = bf_cases(end-1).name;
        rsaf_tmp = bf_cases(end).name;
        
        load([dir_tmpp '/' con_tmp]);
        con_ = cat(3, con_, stSaveInfo.mImg_DSC);
        load([dir_tmpp '/' rsaf_tmp]);
        rsaf_ = cat(3, rsaf_, stSaveInfo.mImg_DSC);
    end
    
    CON = cat(3, CON, con_);
    rSAF = cat(3, rSAF, rsaf_);
end


%% set number of samples
N = [1, 5, 10, 15, 20, 30, 50, 80, 100];

%% compute SSIM for entire image
ssim_CON_single = zeros(numel(N), numel(err_list));
ssim_rSAF_single = zeros(numel(N), numel(err_list));

for n = 1:numel(N)
    N_tmp = N(n);
%     samples_idx = randperm(numel(sample_list), N_tmp);
    samples_idx = 1:N_tmp;
    disp(['>>> N = ' num2str(N_tmp) ]);
    
    for e_idx = 1:numel(err_list)
        disp(['    Error: ' err_list(e_idx).name ]);
        
        rSAF_err = cell2mat(rSAF(e_idx));
        con_err = cell2mat(CON(e_idx));
        
        averaged_frame_con = mean(con_err(:,:,samples_idx), 3);
        averaged_frame_rsaf = mean(rSAF_err(:,:,samples_idx), 3);
        
        %%%%%% normalize the image %%%%%%
        % CON
        img_con_gt_norm = img_con_gt;
        idx_ = find(img_con_gt ~= -30);
        img_con_gt_norm(idx_) = img_con_gt_norm(idx_) + abs(max(img_con_gt_norm(idx_)));
        
        img_con_err_norm = averaged_frame_con;
        idx_ = find(img_con_err_norm ~= -30);
        img_con_err_norm(idx_) = img_con_err_norm(idx_) + abs(max(img_con_err_norm(idx_)));
        
        % rSAF
        img_rsaf_gt_norm = img_rsaf_gt;
        idx_ = find(img_rsaf_gt_norm ~= -30);
        img_rsaf_gt_norm(idx_) = img_rsaf_gt_norm(idx_) + abs(max(img_rsaf_gt_norm(idx_)));
        
        img_rsaf_err_norm = averaged_frame_rsaf;
        idx_ = find(img_rsaf_err_norm ~= -30);
        img_rsaf_err_norm(idx_) = img_rsaf_err_norm(idx_) + abs(max(img_rsaf_err_norm(idx_)));
        
        %%%%%% calculate SSIM %%%%%%
        [ssimval, ssimmap] = ssim(img_con_gt_norm(img_con_gt_norm ~= -30), img_con_err_norm(img_con_err_norm ~= -30));
        ssim_CON_single(n, e_idx) = ssimval;
        
        [ssimval, ssimmap] = ssim(img_rsaf_gt_norm(img_rsaf_gt_norm ~= -30), img_rsaf_err_norm(img_rsaf_err_norm ~= -30));
        ssim_rSAF_single(n, e_idx) = ssimval;
    end
end


%% plot the result
markers = {'*', '+', 'o', 'x', 'square', 'pentagram'};
assert(numel(markers)==numel(err_list), ...
        'Number of markers does not match with err_list');

figure(100);
for e_idx = 1:numel(err_list)
    subplot(1,2,1) % CON
    plot(N, ssim_CON_single(:, e_idx), 'LineWidth', 3, 'Marker', markers{e_idx}); hold on;
    xlabel('N'); ylabel('SSIM value'); title('CON'); grid on;
    
    subplot(1,2,2) % rSAF
    plot(N, ssim_rSAF_single(:, e_idx), 'LineWidth', 3, 'Marker', markers{e_idx}); hold on;
    xlabel('N'); ylabel('SSIM value'); title('rSAF'); grid on;
end
% set ylim
ylim_min = round(min(min(ssim_CON_single(:)), min(ssim_rSAF_single(:)))*10)/10 - 0.05;

set(gcf, 'Position', [-3154 743 787 948]);

subplot(1,2,1); hold off; ylim([ylim_min 1]);

subplot(1,2,2); hold off; ylim([ylim_min 1]);
            
% 
% legend('\epsilon = 0.1 \circ', '\epsilon = 0.2 \circ', '\epsilon = 0.5 \circ', ...
%         '\epsilon = 1 \circ', '\epsilon = 2 \circ', '\epsilon = 5 \circ', ...
%         'Location', 'northeastoutside');



