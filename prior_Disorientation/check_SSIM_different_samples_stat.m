clear; close all; clc;

%%
load_ = 0;

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
deg = 3.75;
line_one = [-stBFInfo.nRadius*sind(deg), -stBFInfo.nRadius*(1-cosd(deg));...
    -(stBFInfo.nRadius+stBFInfo.nDth)*sind(deg), (stBFInfo.nRadius + stBFInfo.nDth)*cosd(deg)-stBFInfo.nRadius]*1e3; % y, z

line_two = [stBFInfo.nRadius*sind(deg), -stBFInfo.nRadius*(1-cosd(deg));...
    (stBFInfo.nRadius+stBFInfo.nDth)*sind(deg), (stBFInfo.nRadius + stBFInfo.nDth)*cosd(deg)-stBFInfo.nRadius]*1e3; % y, z

figure(1);
subplot(1,2,1);imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, img_con_gt);
colormap gray; colorbar; xlabel('Elevational'); ylabel('Axial'); title('CON'); caxis([-50 0]);
subplot(1,2,2);imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, img_rsaf_gt);
colormap gray; colorbar; xlabel('Elevational'); ylabel('Axial'); title('rSAF'); caxis([-50 0]);
set(gcf, 'Position', [-2933 1039 1533 556]);

figure(2);
imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, img_rsaf_gt); hold on;
% line([line_one(1,1) line_one(2,1)], [line_one(1,2) line_one(2,2)], 'color', 'r', 'LineWidth', 1);
% line([line_two(1,1) line_two(2,1)], [line_two(1,2) line_two(2,2)], 'color', 'r', 'LineWidth', 1); hold off;
colormap gray; colorbar; xlabel('Elevational'); ylabel('Axial'); caxis([-50 0]);
set(gcf, 'Position', [-2287 856 788 640]);


%% set number of samples
N = [1, 5, 10, 15, 20, 30, 50, 80, 100];
N = 1;
num_iter = 100;

%%
if(load_)
    load([dir_save '/ssim_CON.mat']);
    load([dir_save '/ssim_rSAF.mat']);
else
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
    
    % compute SSIM for entire image
    % calculate SSIM
    ssim_CON = zeros(numel(N), numel(err_list), num_iter);
    ssim_rSAF = zeros(numel(N), numel(err_list), num_iter);
    
    for i_idx = 1:num_iter
        disp(['Iteration: ' num2str(i_idx)]);tic;
        for n = 1:numel(N)
            N_tmp = N(n);
            samples_idx = randperm(numel(sample_list), N_tmp);
            %         disp(['>>> N = ' num2str(N_tmp) ]);
            
            for e_idx = 1:numel(err_list)
                %             disp(['    Error: ' err_list(e_idx).name ]);
                
                rSAF_err = cell2mat(rSAF(e_idx));
                con_err = cell2mat(CON(e_idx));
                
                averaged_frame_con = mean(con_err(:,:,samples_idx), 3);
                averaged_frame_rsaf = mean(rSAF_err(:,:,samples_idx), 3);
                
                %%%%%% normalize the image %%%%%%
                % CON
                img_con_gt_norm = img_con_gt;
                idx_ = find(img_con_gt ~= -30);
                img_con_gt_norm(idx_) = img_con_gt_norm(idx_) + abs(max(img_con_gt_norm(idx_))); % normalize
                
                img_con_err_norm = averaged_frame_con;
                idx_ = find(img_con_err_norm ~= -30);
                img_con_err_norm(idx_) = img_con_err_norm(idx_) + abs(max(img_con_err_norm(idx_))); % normalize
                
                % rSAF
                img_rsaf_gt_norm = img_rsaf_gt;
                idx_ = find(img_rsaf_gt_norm ~= -30);
                img_rsaf_gt_norm(idx_) = img_rsaf_gt_norm(idx_) + abs(max(img_rsaf_gt_norm(idx_))); % normalize
                
                img_rsaf_err_norm = averaged_frame_rsaf;
                idx_ = find(img_rsaf_err_norm ~= -30);
                img_rsaf_err_norm(idx_) = img_rsaf_err_norm(idx_) + abs(max(img_rsaf_err_norm(idx_))); % normalize
                
                %%%%%% calculate SSIM %%%%%%
                [ssimval, ssimmap] = ssim(img_con_gt_norm(img_con_gt_norm ~= -30), img_con_err_norm(img_con_err_norm ~= -30)); % only take the imaging region
                ssim_CON(n, e_idx, i_idx) = ssimval;
                
                [ssimval, ssimmap] = ssim(img_rsaf_gt_norm(img_rsaf_gt_norm ~= -30), img_rsaf_err_norm(img_rsaf_err_norm ~= -30)); % only take the imaging region
                ssim_rSAF(n, e_idx, i_idx) = ssimval;
            end
        end
        toc;
    end
    
    % save ssim_CON, ssim_rSAF
%     save([dir_save '/ssim_CON'],'ssim_CON', '-v7.3');
%     save([dir_save '/ssim_rSAF'], 'ssim_rSAF', '-v7.3');
end
%% plot the result (N as x-axis)
markers = {'*', '+', 'o', 'x', 'square', 'pentagram'};
assert(numel(markers)==numel(err_list), ...
    'Number of markers does not match with err_list');


color_e = 0.12;


con_mean_N = zeros(numel(N), numel(err_list)); con_std_N = zeros(numel(N), numel(err_list));
rsaf_mean_N = zeros(numel(N), numel(err_list)); rsaf_std_N = zeros(numel(N), numel(err_list));
figure(1);
for e_idx = 1:numel(err_list)
    con_tmp = squeeze(ssim_CON(:,e_idx,:));
    mean_ = mean(con_tmp, 2);
    std_ = std(con_tmp, 0, 2);
    con_mean_N(:,e_idx) = mean_;
    con_std_N(:,e_idx) = std_;
    
    subplot(1,2,1) % CON
    errorbar(N, mean_, std_, 'LineWidth', 2, 'Marker', markers{e_idx}, 'color', [0 0 0]+color_e*e_idx); hold on;
    xlabel('N'); ylabel('SSIM value'); title('CON'); grid on;
    set(gca,'Fontsize',14,'FontWeight','bold');
    
    rsaf_tmp = squeeze(ssim_rSAF(:,e_idx,:));
    mean_ = mean(rsaf_tmp, 2);
    std_ = std(rsaf_tmp, 0, 2);
    rsaf_mean_N(:,e_idx) = mean_;
    rsaf_std_N(:,e_idx) = std_;
    
    subplot(1,2,2) % rSAF
    errorbar(N, mean_, std_, 'LineWidth', 2, 'Marker', markers{e_idx}, 'color', [0 0 0]+color_e*e_idx); hold on;
    xlabel('N'); ylabel('SSIM value'); title('rSAF'); grid on;
    set(gca,'Fontsize',14,'FontWeight','bold');
end
% set ylim
ylim_min = round(min(min(ssim_CON(:)), min(ssim_rSAF(:)))*10)/10 - 0.05;

set(gcf, 'Position', [-3154 743 787 948]);

subplot(1,2,1); hold off; ylim([ylim_min 1]);

subplot(1,2,2); hold off; ylim([ylim_min 1]);
% legend('\sigma = 0.1 \circ', '\sigma = 0.2 \circ', '\sigma = 0.5 \circ', ...
%         '\sigma = 1 \circ', '\sigma = 2 \circ', '\sigma = 5 \circ', ...
%         'Location', 'northeastoutside');

% standard deviation
err_axis = [0.1 0.2 0.5 1 2 5];
figure(2);

subplot(1,2,1); % CON
imagesc(con_std_N); colormap default; caxis([min(min(con_std_N(:)), min(rsaf_std_N(:))) max(max(con_std_N(:)), max(rsaf_std_N(:)))]);
xlabel('\epsilon', 'Interpreter', 'tex'); ylabel('Number of frames for averaging'); title('CON');
set(gca,'Fontsize',20,'FontWeight','bold');

subplot(1,2,2); % rSAF
imagesc(rsaf_std_N); colormap default; caxis([min(min(con_std_N(:)), min(rsaf_std_N(:))) max(max(con_std_N(:)), max(rsaf_std_N(:)))]);
xlabel('\epsilon', 'Interpreter', 'tex'); ylabel('Number of frames for averaging'); title('rSAF');
set(gca,'Fontsize',20,'FontWeight','bold');

set(gcf,'Position',[-2126 887 1027 478]);

%%
errors_ = [0.1, 0.2, 0.5, 1, 2, 5];
figure(4); legend_ = {};
for e_idx = 1:numel(err_list)
    con_tmp = squeeze(ssim_CON(:,e_idx,:));
    mean_ = mean(con_tmp, 2);
    std_ = std(con_tmp, 0, 2);
    con_mean_N(:,e_idx) = mean_;
    con_std_N(:,e_idx) = std_;
    
    % CON
    errorbar(N, mean_, std_, 'LineWidth', 2, 'Marker', markers{e_idx}, 'color', [0 0 0]+color_e*e_idx, 'LineStyle', '--'); hold on;
    xlabel('N'); ylabel('SSIM value'); grid on;
    set(gca,'Fontsize',14,'FontWeight','bold');
    
    legend_{2*(e_idx-1)+1} = ['[CON] \sigma = ' num2str(errors_(e_idx)) '\circ'];
    
    rsaf_tmp = squeeze(ssim_rSAF(:,e_idx,:));
    mean_ = mean(rsaf_tmp, 2);
    std_ = std(rsaf_tmp, 0, 2);
    rsaf_mean_N(:,e_idx) = mean_;
    rsaf_std_N(:,e_idx) = std_;
    
    % rSAF
    errorbar(N, mean_, std_, 'LineWidth', 2, 'Marker', markers{e_idx}, 'color', [0 0 0]+color_e*e_idx);
    xlabel('N'); ylabel('SSIM value');grid on;
    set(gca,'Fontsize',14,'FontWeight','bold');
    
    legend_{2*e_idx} = ['[rSAF] \sigma = ' num2str(errors_(e_idx)) '\circ'];
end
line([0 N(end)], [0.97 0.97], 'color', 'r', 'LineWidth', 2);
% set ylim
ylim_min = round(min(min(ssim_CON(:)), min(ssim_rSAF(:)))*10)/10 - 0.05;
hold off; ylim([ylim_min 1]);
legend(legend_, 'Interpreter', 'tex', 'location', 'northeastoutside');
set(gcf, 'Position', [-3154 743 787 948]);
% ylim([0.9 1]);
%%
figure(6); legend_ = {};

err_ = 5;

e_idx = find(err_axis == err_);
con_tmp = squeeze(ssim_CON(:,e_idx,:));
mean_ = mean(con_tmp, 2);
std_ = std(con_tmp, 0, 2);
con_mean_N(:,e_idx) = mean_;
con_std_N(:,e_idx) = std_;

N_p = linspace(N(1), N(end), 1000);
% CON
errorbar(N, mean_, std_, 'LineWidth', 2, 'Marker', markers{e_idx}, 'LineStyle', '--'); hold on; grid on;

for e_idxx = 1:numel(err_list)
    rsaf_tmp = squeeze(ssim_rSAF(:,e_idxx,:));
    mean_ = mean(rsaf_tmp, 2);
    std_ = std(rsaf_tmp, 0, 2);
    rsaf_mean_N(:,e_idx) = mean_;
    rsaf_std_N(:,e_idx) = std_;
    
    mean_interp = interpn(N, mean_, N_p, 'pchip');
    % rSAF
    plot(N_p, mean_interp);
    
end
xlabel('N'); ylabel('SSIM value'); title('rSAF'); grid on;
set(gca,'Fontsize',14,'FontWeight','bold');

line([0 N(end)], [0.97 0.97], 'color', 'r', 'LineWidth', 2);
% set ylim
hold off;
% ylim([0.8 0.97]);
% legend(legend_, 'Interpreter', 'tex', 'location', 'northeastoutside');
set(gcf, 'Position', [-3154 743 787 948]);

%% frame rate map
CON_0_1 = [5, 10, 15, 20, 30, 50, 80, 100];
rSAF_0_2 = [2.09, 2.486, 2.586, 2.685, 2.784, 2.853, 2.883, 2.982];

CON_0_1_p = [1, 5, 10, 15, 20, 30];
rSAF_0_5 = [5.162, 22.21, 43.22, 51.34, 61.15, 77.6];

CON_0_2 = [15, 20, 30, 50, 80, 100];
rSAF_0_2_p = [1.099, 1.178, 1.297, 1.396, 1.455, 1.495];

CON_0_2_p = [1, 5, 10, 15, 20, 30, 50, 80, 100];
rSAF_0_5_p = [2.685, 5.955, 8.135, 9.225, 9.622, 10.41, 11.21, 11.5, 11.7];

CON_0_5 = [5, 10, 15, 20, 30, 50, 80, 100];
rSAF_0_5_pp = [2.189, 2.784, 2.982, 3.081, 3.279, 3.37, 3.378, 3.477];

CON_0_5_p = [1, 1.892, 2.883, 3.577];
rSAF_1 = [4.468, 7.342, 18.05, 59.96];

CON_1 = [5, 10, 15, 20, 30, 50, 80, 100];
rSAF_1_p = [3.279, 4.072, 4.468, 4.766, 5.162, 5.658, 5.955, 6.054];

CON_2 = [1, 5, 10, 15, 20, 30, 50, 80, 100];
rSAF_2 = [1.694, 6.946, 13.19, 18.34, 23.59, 30, 41.14, 50.55, 52.33];

CON_5 = [1, 2, 3, 4];
rSAF_5 = [3.477, 7.342, 18.24, 58.48];

figure(7);
plot(CON_0_1, rSAF_0_2, 'LineWidth', 2); hold on;
plot(CON_0_1_p, rSAF_0_5, 'LineWidth', 2);
plot(CON_0_2, rSAF_0_2_p, 'LineWidth', 2);
plot(CON_0_2_p, rSAF_0_5_p, 'LineWidth', 2);

plot(CON_0_5, rSAF_0_5_pp, 'LineWidth', 2);
plot(CON_0_5_p, rSAF_1, 'LineWidth', 2);
plot(CON_1, rSAF_1_p, 'LineWidth', 2);
plot(CON_2, rSAF_2, 'LineWidth', 2);
plot(CON_5, rSAF_5, 'LineWidth', 2);

line([0 100], [0 100], 'LineStyle', '--', 'color', 'k', 'LineWidth', 2); hold off;
grid on;
xlabel('CON'); ylabel('rSAF');
legend('\sigma_{CON} = 0.1\circ | \sigma_{rSAF} = 0.2\circ', ...
    '\sigma_{CON} = 0.1\circ | \sigma_{rSAF} = 0.5\circ', ...
    '\sigma_{CON} = 0.2\circ | \sigma_{rSAF} = 0.2\circ', ...
    '\sigma_{CON} = 0.2\circ | \sigma_{rSAF} = 0.5\circ', ...
    '\sigma_{CON} = 0.5\circ | \sigma_{rSAF} = 0.5\circ', ...
    '\sigma_{CON} = 0.5\circ | \sigma_{rSAF} = 1\circ', ...
    '\sigma_{CON} = 1\circ | \sigma_{rSAF} = 1\circ', ...
    '\sigma_{CON} = 2\circ | \sigma_{rSAF} = 2\circ', ...
    '\sigma_{CON} = 5\circ | \sigma_{rSAF} = 5\circ', ...
    'location', 'northeastoutside');
set(gca, 'FontSize', 14);

% for individual plot
figure(8);
plot(CON_0_1, rSAF_0_2, 'color', 'k', 'LineWidth', 2); hold on;
plot(CON_0_1_p, rSAF_0_5, 'color', 'k', 'LineWidth', 2);
plot(CON_0_2, rSAF_0_2_p, 'color', 'k', 'LineWidth', 2);
plot(CON_0_2_p, rSAF_0_5_p, 'color', 'k', 'LineWidth', 2);

line([0 100], [0 100], 'LineStyle', '--', 'LineWidth', 2); hold off;
grid on;
xlabel('CON'); ylabel('rSAF');
% legend('\sigma_{CON} = 0.1\circ | \sigma_{rSAF} = 0.2\circ', ...
%     '\sigma_{CON} = 0.1\circ | \sigma_{rSAF} = 0.5\circ', ...
%     '\sigma_{CON} = 0.2\circ | \sigma_{rSAF} = 0.2\circ', ...
%     '\sigma_{CON} = 0.2\circ | \sigma_{rSAF} = 0.5\circ', ...
%     'location', 'northeastoutside');
set(gca, 'FontSize', 14);


