fc = 10e6;
c = 1540;

lambda_ = c/fc;

r = [5 10 15 20]*1e-3;

fov_oneSide = 66; % deg

% N_syn = 65;
% d_theta = 0.4724;

%% SNR
% N_syn = 61;
% d_theta = 0.3779; % 1/1.25

% N_syn = 97;
% d_theta = 0.3149; % 1/1.5

% N_syn = 113;
% d_theta = 0.2699; % 1/1.75

% N_syn = 131;
% d_theta = 0.2362; % 1/2

%% grating lobe
N_syn = 197;
d_theta = 0.1564;

% N_syn = 125;
% d_theta = 0.2451;

% N_syn = 117;
% d_theta = 0.2644;
%%
% d_alpha = d_theta*pi/180;
d_alpha = sind(d_theta);

z_f_tmp = 10e-3;
% z_f_tmp = [10 15 20 30 40 50 60 70]*1e-3;
% z_f_tmp = [10 20 30 40 50 60 70]*1e-3;
% z_f_tmp = [10 15 20]*1e-3;

% f_elv = 7e-3;
f_elv = 5e-3;

z_ff = z_f_tmp - f_elv;

% y_range = 150e-3;
% y_ = linspace(-y_range, y_range, 500);

for z_idx = 1:numel(z_ff)
    
    figure(101111);
    
    z_f = z_ff(z_idx);
    
    for r_idx = 1:numel(r)
        y_range = tand(fov_oneSide)*(z_f_tmp(z_idx)+r(r_idx));
        y_ = linspace(-y_range, y_range, 3000);
        
        r_tmp = r(r_idx) + f_elv;
        bp_ = sin(pi*r_tmp*N_syn*d_alpha*y_/(lambda_*z_f))./sin(pi*r_tmp*d_alpha*y_/(lambda_*z_f));
        
        subplot(numel(z_ff),numel(r),(z_idx-1)*numel(r)+r_idx);
        plot(atand(y_/(z_f_tmp(z_idx)+r(r_idx))),abs(bp_)); hold on;
%         plot(y_,abs(bp_)); hold on;
        title(['z = ' num2str(z_f_tmp(z_idx)*1e3) 'mm | r = ' num2str(r(r_idx)*1e3) 'mm']);
        xlabel('deg');
        hold off;
        xlim([0 fov_oneSide]);
    end
end
set(gcf,'Position',[561 40 826 1297]);

%%
z_f = 10e-3;

f_elv = 5e-3;
r_ = 15e-3;

% y_ = tand(deg_)*(r_ + z_f);
% y_ = 22.5e-3; % 20mm, -20dB
% y_ = 14.7e-3; % 15mm, -20dB
% y_ = 10.2e-3; % 10mm, -20dB

% gl_theta = 33.8;
% y_ = tand(gl_theta)*(r_+z_f); % 10mm, -20dB
% y_ = 38.5e-3; % 20mm, -20dB
% y_ = 30.8e-3; % 15mm, -20dB
y_ = 21.7e-3; % 10mm, -20dB

fc = 6.5e6;
c = 1540;
lambda_ = c/fc;

d_alpha = lambda_*(z_f-f_elv)/(y_*(r_+f_elv));
d_theta = asind(d_alpha);

deg_ = atand(y_/(z_f + r_));

disp(['grating lobe at ' num2str(deg_) ' deg']);
disp(['required d_theta: ' num2str(d_theta) ' deg']);


%%
fc = 6.5e6;
c = 1540;

lambda_ = c/fc;

r = 15e-3;

N_syn = [65 81 97 113 131];
d_theta = [0.4724 0.3779 0.3149 0.2699 0.2362];
% N_syn = [49 73 79 85 127];
% d_theta = [0.4724 0.3133 0.2908 0.2714 0.1815];

d_alpha = sind(d_theta);

% z_f_tmp = 50e-3;
% z_f_tmp = [10 15 20 30 40 50 60 70]*1e-3;
z_f_tmp = [10 20 30 40 50 60 70]*1e-3;

f_elv = 7e-3;

z_ff = z_f_tmp - f_elv;

fov_oneSide = 66; % deg

for z_idx = 1:numel(z_ff)
    
    y_range = tand(fov_oneSide)*(z_f_tmp(z_idx)+r);
    y_ = linspace(-y_range, y_range, 3000);
    
    figure(191901);
    
    z_f = z_ff(z_idx);
    
    r_tmp = r + f_elv;
    
    subplot(numel(z_ff),1,z_idx);
    for th_idx = 1:numel(d_theta)
        d_alpha_tmp = d_alpha(th_idx);
        N_syn_tmp = N_syn(th_idx);
        bp_ = sin(pi*r_tmp*N_syn_tmp*d_alpha_tmp*y_/(lambda_*z_f))./sin(pi*r_tmp*d_alpha_tmp*y_/(lambda_*z_f));
        
        plot(atand(y_/(z_f_tmp(z_idx)+r)),db(bp_/max(bp_))); hold on;
%         plot(y_,db(bp_/max(bp_))); hold on;
    end
    title(['z = ' num2str(z_f_tmp(z_idx)*1e3) 'mm']);
    xlabel('deg');
    hold off;
    xlim([0 fov_oneSide]); ylim([-40 0]);
%     xlim([0 y_range]);
end
set(gcf,'Position',[513 14 915 1323]);

%% fig5(c)
fc = 6.5e6;
c = 1540;

lambda_ = c/fc;

r = [5 10 15 20]*1e-3;

fov_oneSide = 60; % deg

N_syn = [53 49 45 40];
d_theta = 0.4724;

d_alpha = sind(d_theta);

z_f_tmp = [10 20 30 40 50 60 70]*1e-3;

f_elv = 5e-3;

z_ff = z_f_tmp - f_elv;

for z_idx = 1:numel(z_ff)
    
    figure(10111+z_idx); 
    
    z_f = z_ff(z_idx);
    
    for r_idx = 1:numel(r)
        y_range = tand(fov_oneSide)*(z_f_tmp(z_idx)+r(r_idx));
        y_ = linspace(-y_range, y_range, 5000);
        
        r_tmp = r(r_idx) + f_elv;
        bp_ = sin(pi*r_tmp*N_syn(r_idx)*d_alpha*y_/(lambda_*z_f))./sin(pi*r_tmp*d_alpha*y_/(lambda_*z_f));
        
        plot(atand(y_/(z_f_tmp(z_idx)+r(r_idx))),db(bp_/max(bp_))); hold on;
        title(['z = ' num2str(z_f_tmp(z_idx)*1e3) 'mm']);
        xlabel('deg');
        xlim([0 60]);
        ylim([-35 0]);
    end
    hold off;
    legend('5mm','10mm','15mm','20mm');
    set(gcf,'Position',[464 477 584 205]);
end










