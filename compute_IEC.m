clc; clear; close all;

addpath('functions');
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
beamformed_data = InterpNan(beamformed_data);
env_data = mid_proc(beamformed_data, mid_, acoustic_, bf_);
[axis_y, axis_z, dsc_data] = ScanConverter_convex(env_data, dr, da, bf_.nRadius, height, width, dz, dy);

aROI = find(dsc_data ~= 50);
aOutlier = find(dsc_data == 50);
mOutput = zeros(size(dsc_data));
mOutput(aROI) = dsc_data(aROI);
mOutput = mOutput/max(mOutput(:));

% db scale
mOutput_db = db(mOutput);
mOutput_db(aOutlier) = -30;

%% cyst position
deg_ = [-18 0 18];
depth_cysts = [20 35 50]*1e-3;
offset_ = 1e-3;

depth_cysts = depth_cysts + offset_;

IEC = zeros(1, numel(depth_cysts));
InEn = zeros(1, numel(depth_cysts));

roi_size = [7 30]*1e-3;


for c_idx = 1:numel(depth_cysts)
%     %% set roi region
    cyst_y = (bf_.nRadius + depth_cysts(c_idx)) * sind(deg_(2));
    cyst_z = (bf_.nRadius + depth_cysts(c_idx)) * cosd(deg_(2)) - bf_.nRadius;
    roi_ = [cyst_y-roi_size(2) cyst_z-roi_size(1) 2*roi_size(2) 2*roi_size(1)]*1e3;
    
    %% check roi region
    figure(1); subplot(numel(depth_cysts),2,(c_idx-1)*2+1);
    imagesc(axis_y*1e3,(axis_z-bf_.nRadius)*1e3, mOutput_db); caxis([-60 0]); hold on;
    rectangle('Position',roi_,'EdgeColor','y'); hold off;
    axis tight; axis equal;
    colormap gray;
    
    %% get roi mask
    hor_rng = [-roi_size(2) roi_size(2)] + cyst_y;
    ver_rng = [-roi_size(1) roi_size(1)] + cyst_z;
    
    hor_idx = [find(abs(hor_rng(1)-axis_y)==min(abs(hor_rng(1)-axis_y))) ...
        find(abs(hor_rng(2)-axis_y)==min(abs(hor_rng(2)-axis_y)))];
    ver_idx = [find(abs(ver_rng(1)-axis_z+bf_.nRadius)==min(abs(ver_rng(1)-axis_z+bf_.nRadius))) ...
        find(abs(ver_rng(2)-axis_z+bf_.nRadius)==min(abs(ver_rng(2)-axis_z+bf_.nRadius)))];
    
    roi_data = mOutput(ver_idx(1):ver_idx(2), hor_idx(1):hor_idx(2));
    
    figure(1); subplot(numel(depth_cysts),2,c_idx*2);
    imagesc(axis_y(hor_idx(1):hor_idx(2))*1e3, (axis_z(ver_idx(1):ver_idx(2))-bf_.nRadius)*1e3, roi_data); hold on;
    rectangle('Position',roi_,'EdgeColor','r'); hold off;
    axis equal; axis tight;
    
    %% compute IEC
    InEn(c_idx) = entropy(roi_data); % Information Entropy
%     InEn(c_idx) = entropy(roi_data/max(roi_data(:))); % Information Entropy
%     InEn = entropy(mOutput); % Information Entropy
    
    C_b = 0; % average contrast
    [numRow, numCol] = size(roi_data);
%     [numRow, numCol] = size(mOutput);
    for row_idx = 1:numRow-1
        for col_idx = 1:numCol-1
            tmp = abs(roi_data(row_idx, col_idx)-roi_data(row_idx,col_idx+1));
%             tmp = abs(mOutput(row_idx, col_idx)-mOutput(row_idx,col_idx+1));
            C_b = C_b + tmp;
        end
    end
    C_b = C_b / ((numRow-1)*(numCol-1));
    
    IEC(c_idx) = InEn(c_idx) * C_b;
end
set(gcf,'Position',[561 692 1147 1251]);

disp('InEn');
disp(InEn);



