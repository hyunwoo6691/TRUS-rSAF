clc; close all; clear;

addpath('functions');
%%
[fileName, dir_] = uigetfile('','');

load([dir_ fileName]);

%%
vRcvData = stSaveInfo.Dumped_data;
numScline = stSaveInfo.num_scline;
trans_ = stSaveInfo.Trans;
movingRange = stSaveInfo.moving_range;
movingRangeErr = stSaveInfo.moving_range_err;
errRange = stSaveInfo.err_range;
errors = stSaveInfo.errors;

clear stSaveInfo
%%
stRFInfo.nFc = 6.5e6;
stRFInfo.nC = 1540;
stRFInfo.nFs = stRFInfo.nFc * 4;
stRFInfo.nPitch = trans_.spacingMm*1e-3;
stRFInfo.nLambda = stRFInfo.nC / stRFInfo.nFc;
stRFInfo.nUnitDis = stRFInfo.nC / stRFInfo.nFs;

stTRInfo.nHeight = 5e-3;
stTRInfo.nEleFocus = 20e-3;

stTxInfo.nNoTx = 128;

stBFInfo.sDirection = 'elevational'; % 'elevational' or 'lateral'

stBFInfo.sMode = 'SA'; % Conv: conventional,      SA: synthetic aperture

stBFInfo.sWindow = 'boxcar';
stBFInfo.nScline = stTxInfo.nNoTx;
stBFInfo.nDth = 100e-3;
stBFInfo.nDthSpl = ceil(stBFInfo.nDth/stRFInfo.nUnitDis*2) ;
stBFInfo.nFnum = 1; % receive f-number

stSAInfo.nNoSyn = 5;

if (strcmp(stBFInfo.sDirection, 'elevational'))
    if(stTRInfo.nEleFocus == 0)
        nNaturalFocalDth = stTRInfo.nHeight^2/(4*stRFInfo.nLambda);
        stSAInfo.nVSPos = nNaturalFocalDth;
    else
        stSAInfo.nVSPos = stTRInfo.nEleFocus;
    end
    
    stBFInfo.nRadius = 10e-3;
    stBFInfo.nFOV = 60;          % deg
    stBFInfo.nCh = 1;
else % lateral
    stSAInfo.nVSPos = stTxInfo.nTxFocus;
    stBFInfo.nCh = 64;
end

stMID.nTGC_Atten = 0.5;                                              % [dB]

stMID.nDCRType = 'high';
stMID.nDCRTap = 128;                                             % BPF tap #
stMID.nDCRFcut = 1e6;

nDelayOff = 0;

bSave = 0;
%% Beamforming grid
aSclineTheta = linspace(-0.5*stBFInfo.nFOV, 0.5*stBFInfo.nFOV, stBFInfo.nScline);

aDth = linspace(stBFInfo.nRadius, stBFInfo.nRadius+stBFInfo.nDth, stBFInfo.nDthSpl);

mDth = repmat(aDth', 1, stBFInfo.nScline);
mTheta = repmat(aSclineTheta, numel(aDth), 1);

mImgZ = mDth .* cosd(mTheta);
mImgY = mDth .* sind(mTheta);

%%
ele_idx = 96;
[mBFedData, vBFedData, vSynReg] = fElevational_SA_apod(vRcvData, stRFInfo, stBFInfo, stTRInfo, stSAInfo, aSclineTheta, mImgY, mImgZ, nDelayOff, ele_idx);
                           
%%
env_data = mid_proc(mBFedData, stMID, stRFInfo, stBFInfo);

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

[axis_y, axis_z, dsc_data] = ScanConverter_convex(env_data, dr, da, stBFInfo.nRadius, height, width, dz, dy);

aROI = find(dsc_data ~= 50);
aOutlier = find(dsc_data == 50);
mOutput = zeros(size(dsc_data));
mOutput(aROI) = dsc_data(aROI);
mOutput_db = db(mOutput/max(mOutput(:)));
mOutput_db(aOutlier) = -30;

% figure(1321);
figure;
imagesc(axis_y*1e3,(axis_z-stBFInfo.nRadius)*1e3, mOutput_db); caxis([-40 0]);
axis tight; axis equal; 
% xlabel('Elevational [mm]'); ylabel('Axial [mm]');
colormap gray; %colorbar;
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf,'Position',[493 1115 388 539]);
% set(gcf,'Position',[100 100 330 392]);
%ylim([0 85]); xlim([-20 20]);









