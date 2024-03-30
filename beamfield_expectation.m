%%
load('/Volumes/Research/rSAF/Data_beamfield/FOV65_Radius5mm_Tx2048_Foc25mm/Probability_gauss/error_0.1/Conv.mat')
% load('/Volumes/Research/rSAF/Data_beamfield/FOV60_Radius10mm_Tx2047_Foc20mm/Probability_gauss/error_0.1/Conv.mat')
con_0_1 = prob_;

load('/Volumes/Research/rSAF/Data_beamfield/FOV65_Radius5mm_Tx2048_Foc25mm/Probability_gauss/error_0.1/SA.mat')
% load('/Volumes/Research/rSAF/Data_beamfield/FOV60_Radius10mm_Tx2047_Foc20mm/Probability_gauss/error_0.1/SA.mat')
rsaf_0_1 = prob_;

load('/Volumes/Research/rSAF/Data_beamfield/FOV65_Radius5mm_Tx2048_Foc25mm/Probability_gauss/error_0.2/Conv.mat')
% load('/Volumes/Research/rSAF/Data_beamfield/FOV60_Radius10mm_Tx2047_Foc20mm/Probability_gauss/error_0.2/Conv.mat')
con_0_2 = prob_;

load('/Volumes/Research/rSAF/Data_beamfield/FOV65_Radius5mm_Tx2048_Foc25mm/Probability_gauss/error_0.2/SA.mat')
% load('/Volumes/Research/rSAF/Data_beamfield/FOV60_Radius10mm_Tx2047_Foc20mm/Probability_gauss/error_0.2/SA.mat')
rsaf_0_2 = prob_;

load('/Volumes/Research/rSAF/Data_beamfield/FOV65_Radius5mm_Tx2048_Foc25mm/Probability_gauss/error_0.5/Conv.mat')
% load('/Volumes/Research/rSAF/Data_beamfield/FOV60_Radius10mm_Tx2047_Foc20mm/Probability_gauss/error_0.5/Conv.mat')
con_0_5 = prob_;

load('/Volumes/Research/rSAF/Data_beamfield/FOV65_Radius5mm_Tx2048_Foc25mm/Probability_gauss/error_0.5/SA.mat')
% load('/Volumes/Research/rSAF/Data_beamfield/FOV60_Radius10mm_Tx2047_Foc20mm/Probability_gauss/error_0.5/SA.mat')
rsaf_0_5 = prob_;

load('/Volumes/Research/rSAF/Data_beamfield/FOV65_Radius5mm_Tx2048_Foc25mm/Probability_gauss/error_1/Conv.mat')
% load('/Volumes/Research/rSAF/Data_beamfield/FOV60_Radius10mm_Tx2047_Foc20mm/Probability_gauss/error_1/Conv.mat')
con_1 = prob_;

load('/Volumes/Research/rSAF/Data_beamfield/FOV65_Radius5mm_Tx2048_Foc25mm/Probability_gauss/error_1/SA.mat')
% load('/Volumes/Research/rSAF/Data_beamfield/FOV60_Radius10mm_Tx2047_Foc20mm/Probability_gauss/error_1/SA.mat')
rsaf_1 = prob_;

load('/Volumes/Research/rSAF/Data_beamfield/FOV65_Radius5mm_Tx2048_Foc25mm/Probability_gauss/error_2/Conv.mat')
% load('/Volumes/Research/rSAF/Data_beamfield/FOV60_Radius10mm_Tx2047_Foc20mm/Probability_gauss/error_2/Conv.mat')
con_2 = prob_;

load('/Volumes/Research/rSAF/Data_beamfield/FOV65_Radius5mm_Tx2048_Foc25mm/Probability_gauss/error_2/SA.mat')
% load('/Volumes/Research/rSAF/Data_beamfield/FOV60_Radius10mm_Tx2047_Foc20mm/Probability_gauss/error_2/SA.mat')
rsaf_2 = prob_;

load('/Volumes/Research/rSAF/Data_beamfield/FOV65_Radius5mm_Tx2048_Foc25mm/Probability_gauss/error_5/Conv.mat')
% load('/Volumes/Research/rSAF/Data_beamfield/FOV60_Radius10mm_Tx2047_Foc20mm/Probability_gauss/error_5/Conv.mat')
con_5 = prob_;

load('/Volumes/Research/rSAF/Data_beamfield/FOV65_Radius5mm_Tx2048_Foc25mm/Probability_gauss/error_5/SA.mat')
% load('/Volumes/Research/rSAF/Data_beamfield/FOV60_Radius10mm_Tx2047_Foc20mm/Probability_gauss/error_5/SA.mat')
rsaf_5 = prob_;


load('/Volumes/Research/rSAF/Data_beamfield/FOV65_Radius5mm_Tx2048_Foc25mm/Probability_gauss_noDelay/error_0.1/SA.mat')
con_0_1_fil = prob_;
load('/Volumes/Research/rSAF/Data_beamfield/FOV65_Radius5mm_Tx2048_Foc25mm/Probability_gauss_noDelay/error_0.2/SA.mat')
con_0_2_fil = prob_;
load('/Volumes/Research/rSAF/Data_beamfield/FOV65_Radius5mm_Tx2048_Foc25mm/Probability_gauss_noDelay/error_0.5/SA.mat')
con_0_5_fil = prob_;
load('/Volumes/Research/rSAF/Data_beamfield/FOV65_Radius5mm_Tx2048_Foc25mm/Probability_gauss_noDelay/error_1/SA.mat')
con_1_fil = prob_;
load('/Volumes/Research/rSAF/Data_beamfield/FOV65_Radius5mm_Tx2048_Foc25mm/Probability_gauss_noDelay/error_2/SA.mat')
con_2_fil = prob_;
load('/Volumes/Research/rSAF/Data_beamfield/FOV65_Radius5mm_Tx2048_Foc25mm/Probability_gauss_noDelay/error_5/SA.mat')
con_5_fil = prob_;

clear prob_
%%
axis_ = linspace(-30, 30, 300);

d_idx = 4; % 30mm
% d_idx = 8; % 60mm
% norm_val = max(rsaf_0_1(d_idx,:));
norm_val = 100;

figure(d_idx*100);
plot(axis_(1:end),rsaf_0_1(d_idx,1:end)/norm_val); hold on;
plot(axis_(1:end-1),con_0_2(d_idx,2:end)/norm_val);
plot(axis_(1:end-1),con_0_5(d_idx,2:end)/norm_val);
plot(axis_(1:end-2),con_1(d_idx,3:end)/norm_val);
plot(axis_,con_2(d_idx,:)/norm_val);
plot(axis_,con_5(d_idx,:)/norm_val); hold off;
title('CON');
xlim([-5 5]); ylim([0 0.8]);
legend('0.1','0.2','0.5','1','2','5');

figure(d_idx*10);
plot(axis_(1:end-1),con_0_1(d_idx,2:end)/norm_val); hold on;
plot(axis_(1:end),rsaf_0_2(d_idx,1:end)/norm_val);
plot(axis_(1:end-1),rsaf_0_5(d_idx,2:end)/norm_val);
plot(axis_(1:end-1),rsaf_1(d_idx,2:end)/norm_val);
plot(axis_,rsaf_2(d_idx,:)/norm_val);
plot(axis_,rsaf_5(d_idx,:)/norm_val); hold off;
title('rSAF');
xlim([-5 5]); ylim([0 0.8]);
legend('0.1','0.2','0.5','1','2','5');

figure(d_idx);
plot(axis_(1:end-1),con_0_1(d_idx,2:end)/norm_val); hold on;
plot(axis_,con_0_2_fil(d_idx,:)/norm_val);
plot(axis_(1:end-1),con_0_5_fil(d_idx,2:end)/norm_val);
plot(axis_(1:end-1),con_1_fil(d_idx,2:end)/norm_val);
plot(axis_,con_2_fil(d_idx,:)/norm_val);
plot(axis_,con_5_fil(d_idx,:)/norm_val); hold off;
title('CON-fil');
xlim([-5 5]); ylim([0 0.8]);
legend('0.1','0.2','0.5','1','2','5');


%% width calculation
d_idx = 8;
con_prob = cat(1, rsaf_0_1(d_idx,:),con_0_2(d_idx,:),con_0_5(d_idx,:),con_1(d_idx,:),con_2(d_idx,:),con_5(d_idx,:));
rsaf_prob = cat(1, con_0_1(d_idx,:),rsaf_0_2(d_idx,:),rsaf_0_5(d_idx,:),rsaf_1(d_idx,:),rsaf_2(d_idx,:),rsaf_5(d_idx,:));
% d_idx = 7;
% con_prob = cat(1, rsaf_0_1(d_idx,:),con_0_2(d_idx,:),con_0_5(d_idx,:),con_1(d_idx,:),con_2(d_idx,:),con_5(d_idx,:));
% rsaf_prob = cat(1, con_0_1(d_idx,:),rsaf_0_2(d_idx,:),rsaf_0_5(d_idx,:),rsaf_1(d_idx,:),rsaf_2(d_idx,:),rsaf_5(d_idx,:));

angle_axiss = axis_;
interpn_angle = linspace(-30, 30, 10000);

width_con = zeros(1,6);
width_rsaf = zeros(1,6);
for k = 1:6
% for k = 4
    % CON
    prob_con_d = con_prob(k, :) / 100;
    
    con_interpn = interpn(angle_axiss, prob_con_d, interpn_angle, 'linear');
    Peak_val = max(con_interpn);
    Peak_idx = find(con_interpn == max(con_interpn));
    left_ = con_interpn(1:Peak_idx-1);
    left_angle = interpn_angle(1:Peak_idx-1);
    right_ = con_interpn(Peak_idx+1:end);
    right_angle = interpn_angle(Peak_idx+1:end);
    
    idx_left = find(left_>0.5*Peak_val, 1, 'first');
    left_half = left_(idx_left);
    left_half_angle = left_angle(idx_left);
    
    idx_right = find(right_>0.5*Peak_val, 1, 'last');
    right_half = right_(idx_right);
    right_half_angle = right_angle(idx_right);
    
    width = right_half_angle - left_half_angle;
    width_con(k) = width;
    
    % rSAF
    prob_sa_d = rsaf_prob(k, :) / 100;
    
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
    width_rsaf(k) = width;
end

%%
axis_err = [0.1 0.2 0.5 1 2 5];

width_con = [0.3780	0.5341	1.0621	2.3161	5.2055	13.4623];
width_rsaf = [0.2700	0.4020	0.6361	1.4472	2.8953	8.2589];

figure(11412);
plot(axis_err, width_con,'LineWidth',3,'color','k','Marker','s');hold on;
plot(axis_err, width_rsaf,'LineWidth',3,'color','b','Marker','v');
plot(axis_err, 2*axis_err,'LineWidth',2,'LineStyle','--','color',[0.75 0.75 0.75]);hold off;

%%
rsaf_bf = [0.2700	0.4020	0.6361	1.4472	2.8953	8.2589]';
rsaf_theory = [0.1680	0.3120	0.7441	1.4801	2.9523	7.3847]';

% linear fit
rsaf_theoryy = [ones(length(rsaf_theory),1) rsaf_theory];
coeff_ = rsaf_theoryy \ rsaf_bf;

fitted_line = rsaf_theoryy * coeff_;

figure(14130);
scatter(rsaf_theory(1), rsaf_bf(1),'Marker','v', 'MarkerEdgeColor',[0.15 0.15 0.15],'LineWidth',2,'SizeData',120); hold on;
scatter(rsaf_theory(2), rsaf_bf(2),'Marker','v', 'MarkerEdgeColor',[0.15 0.15 0.15]*2,'LineWidth',2,'SizeData',120);
scatter(rsaf_theory(3), rsaf_bf(3),'Marker','v', 'MarkerEdgeColor',[0.15 0.15 0.15]*3,'LineWidth',2,'SizeData',120);
scatter(rsaf_theory(4), rsaf_bf(4),'Marker','v', 'MarkerEdgeColor',[0.15 0.15 0.15]*4,'LineWidth',2,'SizeData',120);
scatter(rsaf_theory(5), rsaf_bf(5),'Marker','v', 'MarkerEdgeColor',[0.15 0.15 0.15]*5,'LineWidth',2,'SizeData',120);
scatter(rsaf_theory(6), rsaf_bf(6),'Marker','v', 'MarkerEdgeColor',[0.15 0.15 0.15]*6,'LineWidth',2,'SizeData',120);
plot(rsaf_theory, fitted_line, 'LineWidth',3, 'color', 'b'); hold off;
xlim([0 10]); ylim([0 10]);