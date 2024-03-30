clc; clear; close all;

addpath('functions');

bDisplay = 0;
%%
dir_ = uigetdir('./data', 'Select folder');
folder_list = dir(dir_);

flag = 0;
if(strcmp(folder_list(3).name, '.DS_Store')), flag = 1; end
folder_list = folder_list(3+flag:end);

bSave = 1;

%% Setup beamforming parameter
% load info file
load([dir_ '/stSaveInfo.mat']);
stRFInfo = stSaveInfo.stRFInfo;
stTRInfo = stSaveInfo.stTRInfo;
stTxInfo = stSaveInfo.stTxInfo;
nRadius = stSaveInfo.nRadius;
nFOV_Theta_fieldII = stSaveInfo.FOV_theta;
aPosRng_fieldII = fliplr(stSaveInfo.aPosRng);
clear stSaveInfo;

% load data list
dir_data = [dir_ '/Raw data'];
data_list = dir(dir_data);
flag = 0;
if(strcmp(data_list(3).name, '.DS_Store')), flag = 1; end
data_list = data_list(3+flag:end);

dir_save = [dir_ '/errors_bf'];
mkdir(dir_save);

stBFInfo.nDth = 72e-3;
stBFInfo.nDthSpl = ceil(stBFInfo.nDth/stRFInfo.nUnitDis*2);
stBFInfo.nFnum = 1; % receive f-number
stBFInfo.nRadius = nRadius;
stBFInfo.nCh = 128;
stBFInfo.nScline = 128; % number of scanline in lateral direction
stBFInfo.nLatScline = 128; % number of scanline in lateral direction
stBFInfo.sDirection = 'volume';
stBFInfo.sWindow = 'boxcar';

nDelayOff = 2*(sqrt((0.5*stTRInfo.nHeight)^2+stTRInfo.nEleFocus^2) - stTRInfo.nEleFocus)/stRFInfo.nC;

% mid processing paramete
mid_.nTGC_Atten = 0.5;                % [dB]

mid_.nDCRType = 'high';
mid_.nDCRTap = 128;                   % BPF tap #
mid_.nDCRFcut = 1e6;

%% Geometry setup
% transducer geometry
aElePosX = linspace(-0.5*(stTRInfo.nNoEle-1),0.5*(stTRInfo.nNoEle-1),stTRInfo.nNoEle)*stTRInfo.nPitch; % lateral
aElePosZ = zeros(1, stTRInfo.nNoEle);

% Beamforming grid
aDth = linspace(0, stBFInfo.nDth, stBFInfo.nDthSpl); % axial
aXAxis = linspace(aElePosX(1), aElePosX(end), stBFInfo.nLatScline); % lateral

[mImgZ, mImgX] = ndgrid(aDth, aXAxis);
%% save the paramters
if(bSave)
    stParam.stRFInfo = stRFInfo;
    stParam.stBFInfo = stBFInfo;
    stParam.stTRInfo = stTRInfo;
    stParam.stTxInfo = stTxInfo;
    stParam.mImgX = mImgX;
    stParam.mImgZ = mImgZ;
    save([dir_ '/Parameters'], 'stParam', '-v7.3');
end

%% beamforming
vRcvData = [];
for t = 1:stTxInfo.num_tx_plane
    % load raw data
    load([dir_data '/PW_' num2str(t,'%.3d') '.mat']);
    vRcvData = cat(3, vRcvData, mRcvData);
end
%%
% mBFedData -> not normalized
[mBFedData, mSynReg] = fLateral_PWSTF(vRcvData, stRFInfo, stTRInfo, stBFInfo, stTxInfo, aElePosX, aElePosZ, aDth', aXAxis, nDelayOff);

mBFedData = mBFedData./mSynReg;
    
%%
backendTmp = mid_proc(mBFedData, mid_, stRFInfo, stBFInfo);

%%
figure(13);
imagesc(aXAxis*1e3, aDth*1e3, db(backendTmp/max(backendTmp(:))));
colormap gray; caxis([-55 0]);
axis equal; axis tight;
set(gcf,'Position',[242 412 472 630]);
