clc; clear all; close all;

addpath('functions');
bSave = 1;
%% Load data
dir_ = uigetdir('./data','Select folder');
load([dir_ '/stSaveInfo.mat']);

% load info file
stRFInfo = stSaveInfo.stRFInfo;
stTRInfo = stSaveInfo.stTRInfo;
stTxInfo = stSaveInfo.stTxInfo;
clear stSaveInfo;

% load data list
dir_data = [dir_ '/Raw data'];
data_list = dir(dir_data);
flag = 0;
if(strcmp(data_list(3).name, '.DS_Store')), flag = 1; end
data_list = data_list(3+flag:end);


dir_save = [dir_ '/errors_bf'];
mkdir(dir_save);
%
stBFInfo.nDth = 72e-3;
stBFInfo.nDthSpl = ceil(stBFInfo.nDth/stRFInfo.nUnitDis*2);
stBFInfo.nFnum = 1; % receive f-number
stBFInfo.nRadius = stTRInfo.nRadius;
stBFInfo.nCh = 96;
stBFInfo.nFOV = stTRInfo.nFOV;
stBFInfo.nScline = stTRInfo.nNoEle;
stBFInfo.sDirection = 'lateral';
stBFInfo.sWindow = 'none';
stBFInfo.sMode = 'Conv';

d_theta = stTRInfo.nFOV/(stTRInfo.nNoEle-1);

% nDelayOff = 2*(sqrt((0.5*stTRInfo.nHeight)^2+stTRInfo.nEleFocus^2) - stTRInfo.nEleFocus)/stRFInfo.nC;
nDelayOff = 0;

fov = stBFInfo.nFOV;
scline = stBFInfo.nScline;
radius = stBFInfo.nRadius;
depth = stBFInfo.nDth;
depth_spl = stBFInfo.nDthSpl;

aSclineTheta = linspace(-0.5*(stTRInfo.nNoEle-1), 0.5*(stTRInfo.nNoEle-1), stTRInfo.nNoEle)*d_theta; % Ground truth transmitted angle

% Beamforming grid
aDth = linspace(radius, radius+depth, depth_spl);

mDth = repmat(aDth', 1, scline);
mTheta = repmat(aSclineTheta, numel(aDth), 1);

mImgZ = mDth .* cosd(mTheta);
mImgX = mDth .* sind(mTheta);

% linear array for wrong image (ncLUS)
aDth = linspace(radius, radius+depth, depth_spl);
aX = linspace(-0.5*(stTRInfo.nNoEle-1), 0.5*(stTRInfo.nNoEle-1), stTRInfo.nNoEle)*stTRInfo.nPitch;

[mImgZ, mImgX] = ndgrid(aDth, aX);


if(bSave)
    stParam.stRFInfo = stRFInfo;
    stParam.stBFInfo = stBFInfo;
    stParam.stTRInfo = stTRInfo;
    stParam.stTxInfo = stTxInfo;
    stParam.mImgY = mImgX;
    stParam.mImgZ = mImgZ;
    save([dir_ '/Parameters'], 'stParam', '-v7.3');
end
%% load data
vRcvData = [];
for t = 1:scline
    if(mod(t,round(scline/4))==0), disp(['      ' num2str(round(100*t/scline)) '%...']); end
    % load raw data
    load([dir_data '/' data_list(t).name]); % get planewave transmitted at 0 deg
    vRcvData = cat(3, vRcvData, mRcvData);
end

%% Beamforming
switch stBFInfo.sMode
    case 'Conv'
        [mBFedData] = fLateral_Conv( vRcvData, stRFInfo, stTRInfo, stBFInfo, mImgX, mImgZ, nDelayOff);
    case 'SA'
        if(strcmp(stTRInfo.sType, 'linear'))
            [mBFedData, vBFedData, vSynReg, vTD] = fLateral_SA(vRcvData, stRFInfo, stTRInfo, stBFInfo, stSAInfo, stTxInfo, 0, mImgX, mImgZ, nDelayOff);
        else % convex
            
        end
end
%% MID & DSC
% mid processing paramete
mid_.nTGC_Atten = 0.5;                % [dB]

mid_.nDCRType = 'high';
mid_.nDCRTap = 128;                   % BPF tap #
mid_.nDCRFcut = 1e6;

% dsc parameter
scanline_theta = linspace(-0.5*stBFInfo.nFOV, 0.5*stBFInfo.nFOV, stBFInfo.nScline); % Ground truth transmitted angle
depth_ = linspace(stBFInfo.nRadius, stBFInfo.nRadius+stBFInfo.nDth, stBFInfo.nDthSpl);

da = abs(scanline_theta(1)-scanline_theta(2));
dr = abs(depth_(1)-depth_(2));
view_depth = stBFInfo.nDth + (stBFInfo.nRadius*(1-cosd(0.5*stBFInfo.nFOV)));
view_width = 2* (stBFInfo.nRadius+stBFInfo.nDth)*sind(0.5*stBFInfo.nFOV);

dz = 1e-4;
dy = 1e-4;
height = round(view_depth / dz);
width = round(view_width / dy);

env_data = mid_proc(mBFedData, mid_, stRFInfo, stBFInfo);
disp('>>> MID done');


[axis_y, axis_z, dsc_data] = ScanConverter_convex(env_data, dr, da, stBFInfo.nRadius, height, width, dz, dy);

aROI = find(dsc_data ~= 50);
aOutlier = find(dsc_data == 50);
mOutput = zeros(size(dsc_data));
mOutput(aROI) = dsc_data(aROI);
mOutput_db = db(mOutput/max(mOutput(:)));
mOutput_db(aOutlier) = -30;
figure(1);
imagesc(axis_y*1e3,(axis_z-stBFInfo.nRadius)*1e3, mOutput_db); caxis([-40 0]);
axis tight; axis equal; 
% xlabel('Elevational [mm]'); ylabel('Axial [mm]');
colormap gray; %colorbar;
xlim([-30 30]); ylim([15 35]);
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf,'Position',[283 118 707 340]);

%% compute resolution
% % set region of measurement
% % Coordinate setting
% nNoPT = 7;
% nInitialPosY = 0e-3;
% nInitialPosZ = 10e-3;
% nSpacing = 10e-3;
% 
% % ROI position
% aTargetPosY = linspace(nInitialPosY, nInitialPosY, nNoPT);
% aTargetPosZ = linspace(nInitialPosZ, nInitialPosZ+(nNoPT-1)*nSpacing, nNoPT);
% 
% nROISize_Hor = 60e-3;
% nROISize_Ver = 8e-3;
% 
% mROIPosY = zeros(nNoPT, 2);
% mROIPosZ = zeros(nNoPT,2);
% 
% mROIPosY(:,1) = aTargetPosY - 0.5*nROISize_Hor; % left
% mROIPosY(:,2) = aTargetPosY + 0.5*nROISize_Hor; % right
% 
% mROIPosZ(:,1) = aTargetPosZ - 0.5*nROISize_Ver; % upper
% mROIPosZ(:,2) = aTargetPosZ + 0.5*nROISize_Ver; % bottom
% 
% dB = [-6 -6];
% % dB = [-12 -12];
% FWHM_CON = zeros(1, nNoPT);
% 
% for p_idx = 1:nNoPT
%     nLft = mROIPosY(p_idx,1);
%     nRgt = mROIPosY(p_idx,2);
%     nUp = mROIPosZ(p_idx,1);
%     nDn = mROIPosZ(p_idx,2);
%     
%     lIdx = find(abs(axis_y-nLft)==min(abs(axis_y-nLft)));
%     rIdx = find(abs(axis_y-nRgt)==min(abs(axis_y-nRgt)));
%     uIdx = find(abs(axis_z-stBFInfo.nRadius-nUp)==min(abs(axis_z-stBFInfo.nRadius-nUp)));
%     dIdx = find(abs(axis_z-stBFInfo.nRadius-nDn)==min(abs(axis_z-stBFInfo.nRadius-nDn)));
%     
%     roi_ = dsc_data(uIdx:dIdx, lIdx:rIdx);
%     roi_db = zeros(size(roi_));
%     idx_roi = find(roi_ ~= 50);
%     idx_reject = find(roi_ == 50);
%     roi_db(idx_roi) = db(roi_(idx_roi)/max(roi_(idx_roi)));
%     roi_db(idx_reject) = -30;
%     
%     resol_tmp = measure_spatial_resolution(roi_db,axis_y(lIdx:rIdx)*1e3, (axis_z(uIdx:dIdx)-stBFInfo.nRadius)*1e3, dB,1321, nNoPT, p_idx);
%     
%     FWHM_CON(p_idx) = resol_tmp;
% end
% figure(1321);
% set(gcf, 'Position', [743 82 575 971]);