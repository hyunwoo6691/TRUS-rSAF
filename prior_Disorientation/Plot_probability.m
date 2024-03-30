clc; clear; close all;

%%
dir_ = uigetdir('03.Beamfield', '');

dir_save = [dir_ '/Figures'];
mkdir(dir_save);

dir_ = [dir_ '/Probability_gauss'];

err_list = dir(dir_);
flag = 0;
if(strcmp(err_list(3).name, '.DS_Store')), flag = 1; end
err_list = err_list(3+flag:end);

fov = 60;
angle_axis = linspace(-0.5*fov, 0.5*fov, 301);
err_axis = (1:10)*0.1;

depth_list = 1:10;
%% At specific depth
depth_ = 10;

for e_idx = 1:numel(err_list)
    err_tmp = err_list(e_idx).name;
    dir_tmp = [dir_ '/' err_tmp];
    
    load([dir_tmp '/Conv.mat']);
    prob_conv = prob_;
    
    load([dir_tmp '/SA.mat']);
    prob_sa = prob_;
    
    clear prob_;
    
    prob_conv_tmp = prob_conv(depth_,:) / 100;
    prob_sa_tmp = prob_sa(depth_, :) / 100;
    
    prob_conv_d = zeros(1, 301);
    prob_sa_d = zeros(1, 301);
    prob_conv_d(1:300) = prob_conv_tmp;
    prob_sa_d(1:298) = prob_sa_tmp(3:300);
    
    figure(1);sgtitle(['Depth : ' num2str(depth_) 'cm']);
    subplot(2, ceil(numel(err_list)/2), e_idx);
    plot(angle_axis, prob_conv_d, 'LineWidth', 2); hold on;
    plot(angle_axis, prob_sa_d, 'LineWidth', 2); hold off;
    grid on; ylim([0 1]); xlim([0.1*angle_axis(1) 0.1*angle_axis(end)]);
    xlabel('Angle [deg]'); ylabel('Probability'); title(['Error : ' num2str(e_idx/10) 'deg']);
    legend('CON', 'eSAF');
    
end
set(gcf, 'Position', [25 385 1641 563]);

%% Save every figure
depth_roi = [3 6 9];
for e_idx = 1:numel(err_list)
    % for e_idx = (numel(err_list)-1):numel(err_list)
    %     for d_idx = 1:numel(depth_list)
    for d_idx = 1:numel(depth_roi)
        depth = depth_roi(d_idx);
        err_tmp = err_list(e_idx).name;
        dir_tmp = [dir_ '/' err_tmp];
        
        load([dir_tmp '/Conv.mat']);
        prob_conv = prob_;
        
        load([dir_tmp '/SA.mat']);
        prob_sa = prob_;
        
        clear prob_;
        
        prob_conv_d = zeros(1, 301);
        prob_sa_d = zeros(1, 301);
        
        prob_conv_tmp = prob_conv(depth,:) / 100;
        shift_idx = find(prob_conv_tmp == max(prob_conv_tmp))-150;
        shift_idx = shift_idx(end);
        shift_idx = 1;
        if(shift_idx > 0)
            prob_conv_d(1:301-shift_idx) = prob_conv_tmp(shift_idx:end);
        elseif(shift_idx < 0)
            prob_conv_d(1+abs(shift_idx):301) = prob_conv_tmp(1+abs(shift_idx):end);
        else
            prob_conv_d(2:301) = prob_conv_tmp;
        end
        
        prob_sa_tmp = prob_sa(depth, :) / 100;
        shift_idx = find(prob_sa_tmp == max(prob_sa_tmp))-150;
        %         shift_idx = 1;
        if(shift_idx > 0)
            prob_sa_d(1:301-shift_idx) = prob_sa_tmp(shift_idx:end);
        elseif(shift_idx < 0)
            prob_sa_d(1+abs(shift_idx):301) = prob_sa_tmp(1+abs(shift_idx):end);
        else
            prob_sa_d(2:301) = prob_sa_tmp;
        end
        
        %         prob_conv_d(1:300) = prob_conv_tmp(1:end);
        %         prob_sa_d(1:298) = prob_sa_tmp(3:end);
        
        %         fig_ = figure(e_idx+10); %sgtitle(err_tmp, 'Interpreter', 'None');
        %         subplot(1, numel(depth_roi), d_idx);
        figure(depth);
        line([0 0], [0 1], 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--'); hold on;
        plot(angle_axis, prob_conv_d, 'LineWidth', 4); %hold on;
        plot(angle_axis, prob_sa_d, 'LineWidth', 4); hold off;
        %         grid on;
        ylim([0 1]); xlim([-5 5]);
        xlabel('Angle [deg]'); ylabel('Probability'); %title(['Depth : ' num2str(d_idx*10) 'mm']);
        %         legend('CON', 'eSAF');
        set(gca, 'FontSize', 14, 'FontName','Times New Roman', 'FontWeight','bold');
        set(gcf, 'Position', [448 643 425 183]);
    end
    %     set(gcf, 'Position', [31 514 1641 363]);
    %%
    %     saveas(fig_, [dir_save '/[Prob]' err_tmp '.jpg']);
end
%% Overlap epsilon_max at each method
method = {'Conv', 'SA'};
depth_roi = [3 6 9];
err_roi = {'0.1', '0.2', '0.5', '1', '2', '5'};
close all;

for m_idx = 1:numel(method)
    method_tmp = method{m_idx};
    for e_idx = 1:numel(err_roi)
        err_tmp = err_roi{e_idx};
        err_tmp = ['error_' err_tmp];
        for d_idx = 1:numel(depth_roi)
            depth = depth_roi(d_idx);
            
            %             err_tmp = err_list(e_idx).name;
            dir_tmp = [dir_ '/' err_tmp];
            
            load([dir_tmp '/' method_tmp '.mat']);
            
            prob_tmp = prob_(depth,:) / 100;
            
            prob_d = zeros(1, 301);
            
            if(strcmp(method_tmp, 'Conv'))
                shift_idx = 1;
                prob_d(1:301-shift_idx) = prob_tmp(shift_idx:end);
                if(d_idx==3 && m_idx == 1)
                    prob_tmp = prob_(7,:) / 100;
                    prob_d(1:301-shift_idx) = prob_tmp(shift_idx:end);
                end
            else
                shift_idx = find(prob_tmp == max(prob_tmp))-150;
                if(e_idx == 5 || e_idx == 6), shift_idx = 1; end
                
                if(shift_idx > 0)
                    prob_d(1:301-shift_idx) = prob_tmp(shift_idx:end);
                elseif(shift_idx < 0)
                    %                     prob_d(1+abs(shift_idx):301) = prob_tmp(1+abs(shift_idx):end);
                else
                    prob_d(2:301) = prob_tmp;
                end
            end
            
            figure_idx = 10*(m_idx-1) + depth;
            figure(figure_idx); hold on;
            line([0 0], [0 1], 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
            plot(angle_axis, prob_d, 'LineWidth', 3);
            ylim([0 1]); xlim([-5 5]);
            set(gca, 'FontSize', 14, 'FontName','Times New Roman', 'FontWeight','bold');
            set(gcf, 'Position', [448 643 425 183]);
        end
        %          hold off;
    end
end

%% for TMI
method = {'Conv', 'SA'};
depth_roi = [3 6 9];
err_roi = {'0.1', '0.2', '0.5', '1', '2', '5'};
close all;


for m_idx = 1:numel(method)
    method_tmp = method{m_idx};
    for e_idx = 1:numel(err_roi)
        err_tmp = err_roi{e_idx};
        err_tmp = ['error_' err_tmp];
        for d_idx = 1:numel(depth_roi)
            depth = depth_roi(d_idx);
            
            %             err_tmp = err_list(e_idx).name;
            dir_tmp = [dir_ '/' err_tmp];
            
            load([dir_tmp '/' method_tmp '.mat']);
            
            prob_tmp = prob_(depth,:) / 100;
            
            prob_d = zeros(1, 301);
            
            shift_idx = find(prob_tmp == max(prob_tmp))-150;
            if(e_idx == 5 || e_idx == 6), shift_idx = 1; end
            
            if(shift_idx > 0)
                prob_d(1:301-shift_idx) = prob_tmp(shift_idx:end);
            elseif(shift_idx < 0)
                %                     prob_d(1+abs(shift_idx):301) = prob_tmp(1+abs(shift_idx):end);
            else
                prob_d(2:301) = prob_tmp;
            end
            
            % normalize
            prob_d = prob_d;% / max(prob_d);
            
            figure_idx = 10*(m_idx-1) + depth;
            figure(figure_idx); hold on;
            plot(angle_axis, prob_d, 'LineStyle', '-', 'LineWidth', 2,'color', [0.15 0.15 0.15]*e_idx);
            ylim([0 1]); xlim([-5 5]);
            set(gca, 'FontSize', 14, 'FontName','Times New Roman', 'FontWeight','bold');
            set(gcf, 'Position', [448 643 425 183]);
        end
        %          hold off;
    end
end
%% Plot the width (50%) of probabilistic curve
width_conv = zeros(numel(depth_list), numel(err_list));
width_sa = zeros(numel(depth_list), numel(err_list));

angle_axiss = linspace(-0.5*fov, 0.5*fov, 300);
interpn_angle = linspace(-0.5*fov, 0.5*fov, 10000);

for e_idx = 1:numel(err_list)
    % for e_idx = numel(err_list)-1:numel(err_list)
    for d_idx = 1:numel(depth_list)
        err_tmp = err_list(e_idx).name;
        dir_tmp = [dir_ '/' err_tmp];
        
        load([dir_tmp '/Conv.mat']);
        prob_conv = prob_;
        
        load([dir_tmp '/SA.mat']);
        prob_sa = prob_;
        
        clear prob_;
        
        % CON
        prob_conv_d = prob_conv(d_idx,:) / 100;
        conv_interpn = interpn(angle_axiss, prob_conv_d, interpn_angle, 'linear');
        Peak_val = max(conv_interpn);
        Peak_idx = find(conv_interpn == max(conv_interpn));
        left_ = conv_interpn(1:Peak_idx-1);
        left_angle = interpn_angle(1:Peak_idx-1);
        right_ = conv_interpn(Peak_idx+1:end);
        right_angle = interpn_angle(Peak_idx+1:end);
        
        idx_left = find(left_>0.5*Peak_val, 1, 'first');
        left_half = left_(idx_left);
        left_half_angle = left_angle(idx_left);
        
        idx_right = find(right_>0.5*Peak_val, 1, 'last');
        right_half = right_(idx_right);
        right_half_angle = right_angle(idx_right);
        
        width = right_half_angle - left_half_angle;
        width_conv(d_idx, e_idx) = width;
        
        % eSAF
        prob_sa_d = prob_sa(d_idx, :) / 100;
        sa_interpn = interpn(angle_axiss, prob_sa_d, interpn_angle, 'linear');
        Peak_val = max(sa_interpn);
        Peak_idx = find(sa_interpn == max(sa_interpn));
        left_ = sa_interpn(1:Peak_idx-1);
        left_angle = interpn_angle(1:Peak_idx-1);
        right_ = sa_interpn(Peak_idx+1:end);
        right_angle = interpn_angle(Peak_idx+1:end);
        
        idx_left = find(left_>0.5*Peak_val, 1, 'first');
        left_half = left_(idx_left);
        left_half_angle = left_angle(idx_left);
        
        idx_right = find(right_>0.5*Peak_val, 1, 'last');
        right_half = right_(idx_right);
        right_half_angle = right_angle(idx_right);
        
        width = right_half_angle - left_half_angle;
        width_sa(d_idx, e_idx) = width;
    end
    fig_ = figure(e_idx*100);
    plot(depth_list(3:end)*10, width_conv(3:end, e_idx), 'LineWidth', 4, 'Marker', 'o'); hold on;
    plot(depth_list(3:end)*10, width_sa(3:end, e_idx), 'LineWidth', 4, 'Marker', 'd'); hold off;
    grid on; legend('CON', 'eSAF'); xlabel('Depth [mm]'), ylabel('Angle range [deg]');
    title([err_tmp 'deg'], 'Interpreter', 'None');
    %     ylim([0 10]);
    
    set(gcf, 'Position', [560 762 560 186]);
    saveas(fig_, [dir_save '/[Acceptance_range]' err_tmp '.jpg']);
end

% overall plot
err_axis = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 5];
err_axiss = [0.1 0.2 0.5 1 2 5];
% err_axiss = err_axis;

mask_ = zeros(1, numel(err_axis));
for k = 1:numel(err_axiss)
    idx_tmp = find(err_axiss(k) == err_axis);
    mask_(idx_tmp) = 1;
end
mask_ = logical(mask_);

mean_width_conv = mean(width_conv(3:end,:));
mean_width_sa = mean(width_sa(3:end,:));

%%
fig_ = figure(99999);
plot(err_axiss, (mean_width_conv), 'LineWidth', 2.5, 'Marker', '*', 'MarkerSize', 10, 'LineStyle', '--'); hold on;
plot(err_axiss, (mean_width_sa), 'LineWidth', 2.5, 'Marker', 'o', 'MarkerSize', 10, 'LineStyle', '--'); hold off;
grid on; legend('CON', 'eSAF', 'location', 'northeast');
xlabel('\epsilon _{max}'); ylabel('Robustness score');
set(gcf, 'Position', [414 550 697 266]);
set(gca, 'FontSize', 14, 'FontName','Times New Roman', 'FontWeight','bold');
saveas(fig_, [dir_save '/Robustness_score.jpg']);
