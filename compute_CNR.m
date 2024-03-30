clc; clear; close all;

addpath('functions');
%% cyst position
radius_cyst = 4e-3;
deg_ = [-18 0 18];
depth_cysts = [20 35 50]*1e-3;
%%
dir_ = uigetdir('./data','');
%%
load([dir_ '/Parameters.mat']);
acoustic_ = stParam.stRFInfo;
bf_ = stParam.stBFInfo;
trans_ = stParam.stTRInfo;

imgpos_y = stParam.mImgY;
imgpos_z = stParam.mImgZ;

clear stParam;

% mid processing parameter
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
%%
dir_data = [dir_ '/errors_bf/error_0/Sample001/Element_64'];
file_ = dir(dir_data);
file_name = file_(end).name;

load([dir_data '/' file_name]);
beamformed_data = stSaveInfo.mBFedData;
env_data = mid_proc(beamformed_data, mid_, acoustic_, bf_);
[axis_y, axis_z, dsc_data] = ScanConverter_convex(env_data, dr, da, bf_.nRadius, height, width, dz, dy);

aROI = find(dsc_data ~= 50);
aOutlier = find(dsc_data == 50);
mOutput = zeros(size(dsc_data));
mOutput(aROI) = dsc_data(aROI);
mOutput_db = db(mOutput/max(mOutput(:)));
mOutput_db(aOutlier) = -30;

%%
% offset_ = 1e-3;
offset_ = 0;
% bg_deg = [15 11 7]; % REF
bg_deg = [11 9 8];
fraction_ = 1; % REF

% deg_ = [-18 0 18];
deg_ = [-17 0 17]; % rSAF

CNR_ = zeros(1, numel(depth_cysts));
Contrast_ = zeros(1, numel(depth_cysts));
for c_idx = 1:numel(depth_cysts)
    %% set roi region
    cyst_y = (bf_.nRadius + depth_cysts(c_idx)) * sind(deg_(3));
    cyst_z = (bf_.nRadius + depth_cysts(c_idx)) * cosd(deg_(3)) - bf_.nRadius;
    circle_cyst = [cyst_y-radius_cyst+offset_ cyst_z-radius_cyst+offset_ 2*fraction_*radius_cyst 2*fraction_*radius_cyst]*1e3;
    
%     bg_y = (bf_.nRadius + depth_cysts(c_idx)) * sind(deg_(3)+bg_deg(c_idx));
%     bg_z = (bf_.nRadius + depth_cysts(c_idx)) * cosd(deg_(3)+bg_deg(c_idx)) - bf_.nRadius;
    bg_y = (bf_.nRadius + depth_cysts(c_idx)) * sind(deg_(1));
    bg_z = (bf_.nRadius + depth_cysts(c_idx)) * cosd(deg_(1)) - bf_.nRadius;
    circle_bg = [bg_y-radius_cyst+offset_ bg_z-radius_cyst+offset_ 2*fraction_*radius_cyst 2*fraction_*radius_cyst]*1e3;
    %% check roi region
    figure(1);
    subplot(numel(depth_cysts),2,(c_idx-1)*2+1);
    imagesc(axis_y*1e3,(axis_z-bf_.nRadius)*1e3, mOutput_db); caxis([-60 0]); hold on;
    rectangle('Position',circle_cyst,'Curvature', [1 1],'EdgeColor','y');
    rectangle('Position',circle_bg,'Curvature', [1 1],'EdgeColor','y'); hold off;
    axis tight; axis equal;
    colormap gray;
    
    %% get roi mask
    mask_map = zeros(numel(axis_z), numel(axis_y));
    
    [grid_z, grid_y] = ndgrid(axis_z, axis_y);
    
    dist_map_cyst = sqrt((grid_z-bf_.nRadius-cyst_z).^2 + (grid_y-cyst_y).^2);
    dist_map_bg = sqrt((grid_z-bf_.nRadius-bg_z).^2 + (grid_y-bg_y).^2);
    
    mask_cyst = dist_map_cyst<fraction_*radius_cyst;
    mask_bg = dist_map_bg<fraction_*radius_cyst;
    
    mask_map(mask_cyst) = 1;
    mask_map(mask_bg) = 1;
    
    figure(1);
    subplot(numel(depth_cysts),2,c_idx*2);
    imagesc(axis_y*1e3, (axis_z-bf_.nRadius)*1e3,mask_map); hold on;
    rectangle('Position',circle_cyst,'Curvature', [1 1],'EdgeColor','r');
    rectangle('Position',circle_bg,'Curvature', [1 1],'EdgeColor','r'); hold off;
    axis equal; axis tight;
    
    %% extract roi
    cyst_env = dsc_data(mask_cyst);
    bg_env = dsc_data(mask_bg);
    
    CNR_(c_idx) = abs(mean(cyst_env) - mean(bg_env))/sqrt(std(cyst_env)^2 + std(bg_env)^2);
    Contrast_(c_idx) = 10*log10(sum(bg_env.^2)/sum(cyst_env.^2));
end
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf,'Position',[70 48 1147 1251]);

disp('Contrast');
disp(Contrast_);










