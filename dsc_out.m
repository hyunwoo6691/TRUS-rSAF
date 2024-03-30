clear; close all; 
clc;

%%
% max_syn = [49 23 11 7 19]; % f_ele varies
% max_syn = [23 37 51 61 65]; % h varies
% max_syn = [53 49 45]; % R varies
syn_ = 7;

% num_syn = [3 5 7 9 11 13 15 17 19 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65];

% syn_idx = find(num_syn == syn_);
syn_idx = 2;
%%
dir_ = uigetdir('./data', '');

dir_tmp = [dir_ '/errors_bf_/error_0/Sample001/Element_64'];
file_list = dir(dir_tmp);
flag = 0;
if(strcmp(file_list(3).name,'.DS_Store')), flag = 1; end
file_list = file_list(3+flag:end);

load([dir_tmp '/' file_list(syn_idx).name]);
%%
load([dir_ '/Parameters.mat']);
%%
acoustic_ = stParam.stRFInfo;
bf_ = stParam.stBFInfo;
trans_ = stParam.stTRInfo;

imgpos_y = stParam.mImgY;
imgpos_z = stParam.mImgZ;

clear stParam;

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
% mid processing paramete
mid_.nTGC_Atten = 0.5;                % [dB]

mid_.nDCRType = 'high';
mid_.nDCRTap = 128;                   % BPF tap #
mid_.nDCRFcut = 1e6;

beamformed_data = stSaveInfo.mBFedData;
env_data = mid_proc(beamformed_data, mid_, acoustic_, bf_);

[axis_y, axis_z, dsc_data] = ScanConverter_convex(env_data, dr, da, bf_.nRadius, height, width, dz, dy);

aROI = find(dsc_data ~= 50);
aOutlier = find(dsc_data == 50);
mOutput = zeros(size(dsc_data));
mOutput(aROI) = dsc_data(aROI);
mOutput_db = db(mOutput/max(mOutput(:)));
mOutput_db(aOutlier) = -30;

% figure(1321);
figure;
imagesc(axis_y*1e3,(axis_z-bf_.nRadius)*1e3, mOutput_db); clim([-60 0]);
axis tight; axis equal; 
% xlabel('Elevational [mm]'); ylabel('Axial [mm]');
colormap gray; %colorbar;
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf,'Position',[493 1115 388 539]);
% set(gcf,'Position',[100 100 330 392]);
% xlim([-20 20]);



