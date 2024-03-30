clc; clear; close all;
addpath('functions');

dir_ = uigetdir('./data','');
% dir_data = [dir_ '/errors_bf_AWGN_20'];
dir_data = [dir_ '/errors_bf'];
%%
load([dir_ '/Parameters.mat']);
%%
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

%%
sigma_case = [0 0.1 0.2 0.5 1 2 5];
error_list = dir(dir_data);
flag = 0;
if(strcmp(error_list(3).name, '.DS_Store')), flag = 1; end
error_list = error_list(3+flag:end);
assert(numel(sigma_case)==numel(error_list), 'Number of error case does not match');
%%
filter_type = 'boxcar';
%%
fil_ = {};
CON_ = []; CON_fil = []; rSAF_ = [];
for s_idx = 1:numel(sigma_case)
    disp(['>>> ' error_list(s_idx).name]);
%     dir_tmp = [dir_data '/' error_list(s_idx).name '/Sample001/Element_64/'];
    dir_tmp = [dir_data '/' error_list(s_idx).name '/Sample001/Element_77/'];
    bf_list_tmp = dir(dir_tmp);
    load([dir_tmp '/' bf_list_tmp(end-1).name]); % load CON image
    
    beamformed_data = stSaveInfo.mBFedData;
    [beamformed_data_fil, fil] = spatial_filtering(beamformed_data, sigma_case(s_idx), bf_.nDTheta, filter_type);
    fil_ = cat(3, fil_, fil);
    
    load([dir_tmp '/' bf_list_tmp(end).name]); % load rSAF image
%     stSaveInfo.mBFedData = InterpNan(stSaveInfo.mBFedData);
    % stack images
    CON_ = cat(3, CON_, beamformed_data);
    CON_fil = cat(3, CON_fil, beamformed_data_fil);
    rSAF_ = cat(3, rSAF_, stSaveInfo.mBFedData);
    
    env_data = mid_proc(beamformed_data_fil, mid_, acoustic_, bf_);
    
    [axis_y, axis_z, dsc_data] = ScanConverter_convex(env_data, dr, da, bf_.nRadius, height, width, dz, dy);
    
    aROI = find(dsc_data ~= 50);
    aOutlier = find(dsc_data == 50);
    mOutput = zeros(size(dsc_data));
    mOutput(aROI) = dsc_data(aROI);
    mOutput_db = db(mOutput/max(mOutput(:)));
    mOutput_db(aOutlier) = -30;
    
    figure(s_idx);
    imagesc(axis_y*1e3,(axis_z-bf_.nRadius)*1e3, mOutput_db); caxis([-65 0]);
    axis tight; axis equal;
    colormap gray; %colorbar;
    set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
    set(gcf,'Position',[818 45 1031 676]);
end

%% 1-d psf
dor = 3e-3;
depth = 22e-3;

dist_map = sqrt((imgpos_z-bf_.nRadius).^2 + imgpos_y.^2);

mask_map = logical(zeros(size(dist_map)));

mask_map(dist_map>depth-dor & dist_map<depth+dor) = true;

figure(3131);imagesc(mask_map);
title('mask map');

for s_idx = 1:numel(sigma_case)
    roi_con = zeros(size(mask_map));
    roi_con_fil = zeros(size(mask_map));
    roi_rsaf = zeros(size(mask_map));
    
    con_tmp = abs(CON_(:,:,s_idx));
    con_fil_tmp = abs(CON_fil(:,:,s_idx));
    rsaf_tmp = abs(rSAF_(:,:,s_idx));
    
    roi_con(mask_map) = con_tmp(mask_map);
    roi_con_fil(mask_map) = con_fil_tmp(mask_map);
    roi_rsaf(mask_map) = rsaf_tmp(mask_map);
    
    con_1d = mean(roi_con,1);
    con_fil_1d = mean(roi_con_fil,1);
    rsaf_1d = mean(roi_rsaf,1);
    
    figure(s_idx);
    plot(db(con_1d/max(con_1d)),'LineWidth',2); hold on;
    plot(db(con_fil_1d/max(con_fil_1d)),'LineWidth',2);
    plot(db(rsaf_1d/max(rsaf_1d)),'LineWidth',2); hold off;
    legend('CON', 'CON+fil', 'rSAF');
    ylim([-60 0])
end



