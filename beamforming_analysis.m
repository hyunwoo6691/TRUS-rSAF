clc; clear; close all;

%%
dir_ = uigetdir('./data','');

% data
dir_data = [dir_ '/errors_bf'];
folder_list = dir(dir_data);
flag = 0;
if(strcmp(folder_list(3).name,'.DS_Store')), flag = 1; end
folder_list = folder_list(3+flag:end);

% load parameters
load([dir_ '/Parameters.mat']);
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

%% show image for picking pixel
disp('>>> load ground truth data');
load([dir_data '/error_0/Sample001/Element_64/[Save]SA_VS5_Syn131.mat']);
beamformed_data = stSaveInfo.mBFedData; clear stSaveInfo
env_data = mid_proc(beamformed_data, mid_, acoustic_, bf_);

[axis_y, axis_z, dsc_data] = ScanConverter_convex(env_data, dr, da, bf_.nRadius, height, width, dz, dy);

aROI = find(dsc_data ~= 50);
aOutlier = find(dsc_data == 50);
mOutput = zeros(size(dsc_data));
mOutput(aROI) = dsc_data(aROI);
mOutput_db = db(mOutput/max(mOutput(:)));
mOutput_db(aOutlier) = -30;

%%
figure(1);
imagesc(axis_y*1e3,(axis_z-bf_.nRadius)*1e3, mOutput_db); caxis([-60 0]);
axis tight; axis equal;
% xlabel('Elevational [mm]'); ylabel('Axial [mm]');
colormap gray; %colorbar;
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf,'Position',[818 45 1031 676]);

roi = drawpoint;
%% calculate the index in beamforming coordinate
target_point = roi.Position * 1e-3; % y, z

theta_ = atand(target_point(1)/(target_point(2)+bf_.nRadius));
r_ = sqrt(target_point(1)^2 + (target_point(2)+bf_.nRadius)^2);

imgpos_r = sqrt(imgpos_y.^2 + imgpos_z.^2);
imgpos_theta = atand(imgpos_y./imgpos_z);

imgpos_r = imgpos_r(:,1);
imgpos_theta = imgpos_theta(1,:);

r_idx = find(abs(imgpos_r - r_) == min(abs(imgpos_r - r_)));
theta_idx = find(abs(imgpos_theta - theta_) == min(abs(imgpos_theta - theta_)));

target_idx = [r_idx, theta_idx]; % z, y

%%
num_syn = 131;
beam_idx = max(theta_idx-(num_syn-1)*0.5,1):min(theta_idx+(num_syn-1)*0.5,bf_.nScline);
%%
num_beam = bf_.nScline;
beam_data = zeros(numel(folder_list), num_beam);

beamformed_data_ = [];
for f_idx = 1:numel(folder_list)
    folder_tmp = folder_list(f_idx).name;
    disp(['>>> ' folder_tmp]);
    
    dir_tmp = [dir_data '/' folder_tmp '/Sample001/Element_64'];
    file_ = dir(dir_tmp);
    load([dir_tmp '/' file_(end).name]);
    
    beamformed_data = stSaveInfo.mBFedData;
    beamformed_data_ = cat(3, beamformed_data_, beamformed_data);
    
    beamformed_vol = stSaveInfo.vBFedData; % sample * beam * scanline
    syn_reg_vol = stSaveInfo.vSynReg; % sample * scanline * beam
    
    clear stSaveInfo
    
    
    for k = 1:num_beam
        beam_data(f_idx,k) = beamformed_vol(target_idx(1),k,target_idx(2)) * syn_reg_vol(target_idx(1), target_idx(2), k);
    end
end
disp('    done');
%%
figure(1321);
for k = 1:numel(folder_list)
    plot(beam_idx, beam_data(k,beam_idx),'color',[0.13 0.13 0.13]*k, 'LineWidth',2); hold on;
end
hold off;
% line([theta_idx theta_idx], [max(beam_data(:)) min(beam_data(:))],'color','r'); hold off;

%%
mOutput_ = [];
for f_idx = 1:numel(folder_list)
    env_data = mid_proc(beamformed_data_(:,:,f_idx), mid_, acoustic_, bf_);
    
    [axis_y, axis_z, dsc_data_tmp] = ScanConverter_convex(env_data, dr, da, bf_.nRadius, height, width, dz, dy);
    
    aROI = find(dsc_data_tmp ~= 50);
    aOutlier = find(dsc_data_tmp == 50);
    mOutput = zeros(size(dsc_data_tmp));
    mOutput(aROI) = dsc_data_tmp(aROI);
    
    if(f_idx == 1), norm_val = max(mOutput(:)); end
    
    mOutput_ = cat(3, mOutput_, mOutput);
end

for f_idx = 1:numel(folder_list)
    mOutput_tmp = mOutput_(:,:,f_idx);
    
    mOutput_db = db(mOutput_tmp/norm_val);
    mOutput_db(aOutlier) = -30;
    
    figure(f_idx);
    imagesc(axis_y*1e3,(axis_z-bf_.nRadius)*1e3, mOutput_db); caxis([-60 0]);
    axis tight; axis equal;
    colormap gray; %colorbar;
    set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
    set(gcf,'Position',[818 45 1031 676]);
end
