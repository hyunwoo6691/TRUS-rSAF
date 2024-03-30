clc; clear; close all;

%%
errors_deg = [0.1 0.2 0.5 1 2 5];
errors_rad = errors_deg * pi/180;
sigma_ = errors_rad;

r = 5e-3;

fc = 6.5e6;
c = 1540;
lambda_ = c/fc;
k_ = 2*pi/lambda_;

height = 5e-3;

z_f_ = [3 6 9]*1e-2+1e-3;
y_f = 0e-3;

num_syn = 4;
delta_alpha = (60/127)*pi/180
alpha_max = delta_alpha*(num_syn-1);
deg_max = alpha_max*180/pi

y_axis = linspace(-40e-3, 40e-3, 1000);

%% y-axis
width_theory = zeros(numel(z_f_), numel(sigma_));
interpn_yaxis = linspace(-40e-3, 40e-3, 10000);
z_idxx = 1;
for z_idx = z_f_
    z_f = z_idx;
    
    figure(20+round(z_f*1e2));
    for s_idx = 1:numel(sigma_)
        sigma_tmp = sigma_(s_idx);
        c_0 = sigma_tmp * sqrt(2*pi*z_f/(z_f-1j*k_*(z_f+r)*r*sigma_tmp^2)); % c_2
        q_ = 1j*k_*(r*(y_axis-y_f))/z_f - ((k_*sigma_tmp*r)^2)*y_axis*(z_f+r)/(z_f*(z_f-1j*k_*sigma_tmp^2*(z_f+r)*r)); % q_2
        p_ = (k_*sigma_tmp*r)^2*(z_f+r)^2/(2*z_f*(z_f - 1j*k_*(sigma_tmp^2)*(z_f+r)*r)); % q_1
        r_ = 1j*k_*(y_axis.^2-y_f^2)/(2*z_f) + (k_*sigma_tmp*r*y_axis).^2/(2*z_f*(z_f-1j*k_*(sigma_tmp^2)*(z_f+r)*r)); % q_3
        E_h = abs(c_0 * exp(0.25*(q_.^2)/p_-r_) * sqrt(pi*(1-exp(-0.25*p_*alpha_max^2))/p_)/(sigma_tmp*sqrt(2*pi)));
        
        plot(y_axis*1e3, E_h/max(E_h), 'LineWidth',2, 'color', [0.15 0.15 0.15]*s_idx); hold on; % cartesian view
        
        %% probability width calculation
        norm_Eh = E_h/max(E_h);
        norm_Eh = interpn(y_axis, norm_Eh, interpn_yaxis, 'linear');
        peak_idx = find(norm_Eh == 1 , 1, 'first');
        left_ = norm_Eh(1:peak_idx-1);
        right_ = norm_Eh(peak_idx+1:end);
        
        y_left = interpn_yaxis(1:peak_idx-1);
        y_right = interpn_yaxis(peak_idx+1:end);
        
        left_half = find(left_ >= 0.5, 1, 'first');
        right_half = find(right_ >= 0.5, 1, 'last');
        
        width_theory(z_idxx, s_idx) = (y_right(right_half) - y_left(left_half))*1e3;
        
    end
    hold off; %grid on;
%     legend('\sigma = 0.1\circ', '\sigma = 0.2\circ', '\sigma = 0.5\circ',...
%         '\sigma = 1\circ', '\sigma = 2\circ', '\sigma = 5\circ',...
%         'location', 'northeastoutside');
    %     xlabel('[ mm ]');
    xlim([-5 5]);
    %     ylabel('E[h(\epsilon)]'); title(['(y_f, z_f) : (' num2str(y_f*1e3) ', ' num2str(z_f*1e3) ')']);
    set(gca,'FontSize',14);
    set(gcf, 'Position', [-2868 666 560 420]);
    z_idxx = z_idxx + 1;
end
%% plot width
% mean_width_theory = mean(width_theory);
mean_width_theory = width_theory(1,:);
figure(11312); 
plot(errors_deg, mean_width_theory);

% correlation with experiment
depth_roi = [3 6 9];
mean_width_rsaf = mean(width_sa(depth_roi, :));

figure(123142);
plot(errors_deg, mean_width_theory, 'LineWidth', 2.5, 'color', 'k', 'Marker', 'square', 'MarkerSize',10); hold on;
plot(errors_deg, mean_width_rsaf, 'LineWidth', 2.5, 'color', 'b', 'Marker', 'v', 'MarkerSize',10); hold off;
legend('rSAF-CW', 'rSAF-PW');
%% radial-axis
rad_axis = linspace(-5, 5, 1000); % deg
for z_idx = z_f_
    z_f = z_idx;
    
    y_axis_ = tand(rad_axis) * z_f;
    
    figure(round(z_f*1e2)*10);
    for s_idx = 1:numel(sigma_)
        sigma_tmp = sigma_(s_idx);
        c_0 = sigma_tmp * sqrt(2*pi*z_f/(z_f-1j*k_*(z_f+r)*r*sigma_tmp^2)); % c_2
        q_ = 1j*k_*(r*(y_axis_-y_f))/z_f - ((k_*sigma_tmp*r)^2)*y_axis_*(z_f+r)/(z_f*(z_f-1j*k_*sigma_tmp^2*(z_f+r)*r)); % q_2
        p_ = (k_*sigma_tmp*r)^2*(z_f+r)^2/(2*z_f*(z_f - 1j*k_*(sigma_tmp^2)*(z_f+r)*r)); % q_1
        r_ = 1j*k_*(y_axis_.^2-y_f^2)/(2*z_f) + (k_*sigma_tmp*r*y_axis_).^2/(2*z_f*(z_f-1j*k_*(sigma_tmp^2)*(z_f+r)*r)); % q_3
        E_h = abs(c_0 * exp(0.25*(q_.^2)/p_-r_) * sqrt(pi*(1-exp(-0.25*p_*alpha_max^2))/p_)/(sigma_tmp*sqrt(2*pi)));
        
%         plot(y_axis_*1e3, E_h/max(E_h), 'LineWidth',2, 'color', [0.15 0.15 0.15]*s_idx); hold on; % cartesian view
        plot(rad_axis, E_h/max(E_h), 'LineWidth',2, 'color', [0.15 0.15 0.15]*s_idx); hold on; % radial view
    end
    hold off; grid on;
%     legend('\sigma = 0.1\circ', '\sigma = 0.2\circ', '\sigma = 0.5\circ',...
%         '\sigma = 1\circ', '\sigma = 2\circ', '\sigma = 5\circ',...
%         'location', 'northeastoutside');
    %     xlabel('[ mm ]');
    xlim([-5 5]);
    %     ylabel('E[h(\epsilon)]'); title(['(y_f, z_f) : (' num2str(y_f*1e3) ', ' num2str(z_f*1e3) ')']);
    set(gca,'FontSize',14);
    set(gcf, 'Position', [-2868 666 560 420]);
end

%% validation of Far-field approximation
r = 5e-3;
deg = 45; % degree
alpha_ = sind(deg);
beta_ = cosd(deg);

R = linspace(r, r+100e-3, 1000);

y = R*alpha_;
z = R*beta_;

% focal point
y_f = 0;
% z_f = linspace(0.01e-3, z_max, 1000);
z_f = 10e-3;

ground_truth = sqrt((y-r*alpha_).^2+(z-r*beta_).^2) - sqrt((y_f-r*alpha_).^2+(z_f-r*beta_).^2);

approximation = (z-r*beta_).*(1+0.5*((y-r*alpha_)./(z-r*beta_)).^2) - (z_f-r*beta_).*(1+0.5*((y_f-r*alpha_)./(z_f-r*beta_)).^2);
% approximation = (y.^2 - y_f^2)./(2*(z_f-r*beta_)) - (r*(y-y_f)*alpha_)./(z_f-r*beta_);

figure(1);
plot(R*1e3, ground_truth); hold on;
% plot(R*1e3, z-z_f, 'color', 'k');
plot(R*1e3, approximation); hold off;
xlabel('mm');grid on; %xlim([r*1e2+0.05 z_max*1e2]);
legend('ground truth', 'approximation');

%% 2D
r = 5e-3;
degrees = linspace(0, 45, 46); % degree
R = linspace(r, r+100e-3, 100);

% focal point
y_f = 0;
% z_f = linspace(0.01e-3, z_max, 1000);
z_f = 30e-3;

error_map = zeros(numel(R), numel(degrees));

for d_idx = 1:numel(degrees)
    deg = degrees(d_idx);
    
    alpha_ = sind(deg);
    beta_ = cosd(deg);
    
    y = R*alpha_;
    z = R*beta_;
    
    ground_truth = sqrt((y-r*alpha_).^2+(z-r*beta_).^2) - sqrt((y_f-r*alpha_).^2+(z_f-r*beta_).^2);
    
    approximation = (z-r*beta_).*(1+0.5*((y-r*alpha_)./(z-r*beta_)).^2) - (z_f-r*beta_).*(1+0.5*((y_f-r*alpha_)./(z_f-r*beta_)).^2);
    % approximation = (y.^2 - y_f^2)./(2*(z_f-r*beta_)) - (r*(y-y_f)*alpha_)./(z_f-r*beta_);
    
    error = abs(ground_truth - approximation)/lambda_;
    error_map(:,d_idx) = error;
end

figure(1);
imagesc(degrees, R*1e3, error_map); 
xlabel('\theta [\circ]'); ylabel('R [mm]');
colorbar; colormap jet;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');