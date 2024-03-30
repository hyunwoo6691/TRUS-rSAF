clc; close all; clear;
%%
%%%%%%%%%%%%%%% dB_ select %%%%%%%%%%%%%%%%%%%%%%%
dB_ = -6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
dir_master = uigetdir('./01.Data', '');
dir_master = [dir_master '/Data'];
folder_list = dir(dir_master);
flag = 0;
if(strcmp(folder_list(3).name, '.DS_Store'))
    flag = 1;
end
folder_list = folder_list(3+flag:end);

for f_idx = 1:numel(folder_list)
    %% Load data
    %     dir_ = uigetdir('./01.Data', '');
    folder_name = folder_list(f_idx).name;
    dir_ = [dir_master '/' folder_name];
    
    disp(['>>> Folder name: ' folder_name]);
    
    dir_save = [dir_ '/Figures_' num2str(dB_) 'dB'];
    mkdir(dir_save);
    
    % load ground truth
    dir_gt = [dir_ '/Resolutions_' num2str(dB_) 'dB/error_0'];
    load([dir_gt '/Sample001.mat']);
    resol_gt = mLateralResol;
    
    % load error
    error_list = dir([dir_ '/Resolutions_' num2str(dB_) 'dB']);
    flag = 0;
    if(strcmp(error_list(3).name, '.DS_Store'))
        flag = 1;
    end
    error_list = error_list(4+flag:end);
    
    resol_err = {};
    for e_idx = 1:numel(error_list)
        resol_tmp = zeros(size(resol_gt,1), size(resol_gt,2), 100);
        error_tmp = error_list(e_idx).name;
        dir_tmp = [dir_ '/Resolutions_' num2str(dB_) 'dB/' error_tmp];
        sample_list = dir(dir_tmp);
        flag = 0;
        if(strcmp(sample_list(3).name, '.DS_Store'))
            flag = 1;
        end
        sample_list = sample_list(3+flag:end);
        
        for s_idx = 1:numel(sample_list)
            sample_tmp = sample_list(s_idx).name;
            load([dir_tmp '/' sample_tmp]);
            resol_tmp(:,:,s_idx) = mLateralResol;
        end
        
        resol_err = cat(3, resol_err, resol_tmp);
    end
    clear mLateralResol;
    clear resol_tmp;
    
    %% plot1 - every error, along depth
%     
    z_axis = (1:size(resol_gt,2)) * 1e-2;
%     for bf_idx = 1:2
%         for e_idx = 1:numel(error_list)
%             err_tmp = resol_err{e_idx};
%             err_tmpp = squeeze(err_tmp(bf_idx, :, :))';
%             
%             ground_truth = resol_gt(bf_idx, :);
%             
%             fig_ = figure(e_idx);
%             plot(ground_truth,'LineWidth', 2, 'Marker', 'o', 'LineStyle','--'); hold on;
%             boxplot(err_tmpp); hold off;
%             xlabel('Depth [cm]'); ylabel('Elevational resolution [mm]'); grid minor;
%             ylim([0 50]);
%             leg = legend('Ground truth','Location','northwest');
%             set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
%             set(gcf,'Position',[126 87 800 300]);
%             
%             if(bf_idx == 1)
%                 saveas(fig_, [dir_save '/[Conv]' error_list(e_idx).name '.jpg']);
%             else
%                 saveas(fig_, [dir_save '/[SA]' error_list(e_idx).name '.jpg']);
%             end
%         end
%     end
%     clear err_tmpp;
%     close all;
    
    %% Pick worst & mean error, standard deviation among 100 samples
    gt_tmp = repmat(resol_gt, 1, 1, 100);
    
    worst_case = zeros(numel(error_list), size(resol_gt,2), 2); % row : error_case    col: depth
    mean_case = zeros(numel(error_list), size(resol_gt,2), 2);
    std_ = zeros(numel(error_list), size(resol_gt,2), 2);
    
    for bf_idx = 1:2
        gt_tmpp = squeeze(gt_tmp(bf_idx, :, :))';
        for e_idx = 1:numel(error_list)
            err_tmp = resol_err{e_idx};
            err_tmpp  = squeeze(err_tmp(bf_idx, :, :))';
            
            error_ = abs(gt_tmpp - err_tmpp);
%             error_ = gt_tmpp - err_tmpp;
            
            worst_case(e_idx, :, bf_idx) = max(error_)./resol_gt(bf_idx,:) * 100; % percent error
            mean_case(e_idx, :, bf_idx) = mean(error_)./resol_gt(bf_idx,:) * 100; % percent error
            std_(e_idx, :, bf_idx) = std(error_)./resol_gt(bf_idx,:)*100;
        end
    end
    
    % for every depth
    err_axis = (1:10) / 10;
    for d_idx = 1:size(resol_gt,2)
        
        % worst case
%         conv_tmp = worst_case(:, d_idx, 1);
%         sa_tmp = worst_case(:, d_idx, 2);
%         
%         fig_ = figure(d_idx);
%         plot(err_axis, conv_tmp, 'LineWidth', 4, 'Marker', 'o', 'MarkerSize', 20, 'LineStyle','--'); hold on;
%         plot(err_axis, sa_tmp, 'LineWidth', 4, 'Marker', '*', 'MarkerSize', 20, 'LineStyle','--'); hold off;
%         
%         xlabel('E_max', 'Interpreter', 'None'); ylabel('Worst case error [%]'); grid minor;
%         %     ylim([0 200]);
%         leg = legend('Conventional', 'Synthetic Aperture','Location','northwest');
%         set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
%         set(gcf,'Position',[126 87 800 300]);
%         
%         saveas(fig_, [dir_save '/[Compare_worst]depth_' num2str(d_idx) 'cm.jpg']);
        
        % mean
        conv_tmp = mean_case(:, d_idx, 1);
        sa_tmp = mean_case(:, d_idx, 2);
        std_tmp_conv = std_(:, d_idx, 1);
        std_tmp_sa = std_(:, d_idx, 2);
        
        fig_ = figure(d_idx);
%         plot(err_axis, conv_tmp, 'LineWidth', 4, 'Marker', 'o', 'MarkerSize', 20, 'LineStyle','--'); hold on;
%         plot(err_axis, sa_tmp, 'LineWidth', 4, 'Marker', '*', 'MarkerSize', 20, 'LineStyle','--'); 
%         scatter(err_axis, std_tmp_conv, 'Marker', '+', 'SizeData', 50, 'LineWidth', 3);
%         scatter(err_axis, std_tmp_sa, 'Marker', 'x', 'SizeData', 50, 'LineWidth', 3); hold off;
        errorbar(err_axis, conv_tmp, std_tmp_conv, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 10, 'LineStyle','--'); hold on;
        errorbar(err_axis, sa_tmp, std_tmp_sa, 'LineWidth', 2, 'Marker', '*', 'MarkerSize', 10, 'LineStyle','--'); hold off;

        xlabel('E_max', 'Interpreter', 'None'); ylabel('Mean error [%]'); grid minor;
%         ylim([0 50]);
        leg = legend('Conventional', 'Synthetic Aperture', 'std_conv', 'std_sa', 'Location','northwest', 'Interpreter', 'None');
        set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
        set(gcf,'Position',[126 87 800 300]);
        
        saveas(fig_, [dir_save '/[Errorbar]depth_' num2str(d_idx) 'cm.jpg']);
        
    end
    
    % standard deviation
    % text file
%     fileID = fopen([dir_ '/std.txt'], 'w');
%     fprintf(fileID, 'Standard deviation along depth \n');
%     fprintf(fileID, 'Conventional\n');
%     for d_idx = 1:size(resol_gt,2)
%         fprintf(fileID, [num2str(d_idx) ' cm \n']);
%         fprintf(fileID, [num2str(std_(:,d_idx,1)) '\n']);
%     end
%     fprintf('Synthetic aperture \n');
%     for d_idx = 1:size(resol_gt,2)
%         fprintf(fileID, [num2str(d_idx) ' cm \n']);
%         fprintf(fileID, [num2str(std_(:,d_idx,2)) '\n']);
%     end
    % figure
%     fig_ = figure(999999);
%     subplot(1,2,1); % conventional
%     imagesc(z_axis*1e2, err_axis, std_(:,:,1)); xlabel('depth [cm]'); ylabel('error [deg]'); title('Conventional');
%     colorbar; caxis([0 max(max(std_(:,:,1)))]);
%     subplot(1,2,2); % synthetic aperture
%     imagesc(z_axis*1e2, err_axis, std_(:,:,2)); xlabel('depth [cm]'); ylabel('error [deg]'); title('Synthetic aperture');
%     colorbar; caxis([0 max(max(std_(:,:,1)))]);
%     set(gcf, 'Position', [560 482 1069 466]);
%     
%     saveas(fig_, [dir_save '/Standard_deviation.jpg']);
end

close all;
