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

num_syn = 131;

mode = 'SA';
% mode = 'Conv';

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

stBFInfo.nDth = 85e-3;
stBFInfo.nDthSpl = ceil(stBFInfo.nDth/stRFInfo.nUnitDis*2);
stBFInfo.nFnum = 1; % receive f-number
stBFInfo.nRadius = nRadius;
stBFInfo.nCh = stTRInfo.nNoEle;
stBFInfo.nFOV = nFOV_Theta_fieldII;
% stBFInfo.nDTheta = 0.4724;
stBFInfo.nDTheta = 0.2362;
stBFInfo.nRadialScline = numel(aPosRng_fieldII); % number of scanline in radial direction
stBFInfo.nLatScline = stTRInfo.nNoEle; % number of scanline in lateral direction
stBFInfo.sDirection = 'volume';
stBFInfo.sWindow = 'boxcar';
stBFInfo.sMode = mode;

stSAInfo.nNoSyn = num_syn;
if(stTRInfo.nEleFocus == 0)
    nNaturalFocalDth = stTRInfo.nHeight^2/(4*stRFInfo.nLambda);
    stSAInfo.nVSPos = nNaturalFocalDth;
else
    stSAInfo.nVSPos = stTRInfo.nEleFocus;
end

nDelayOff = 2*(sqrt((0.5*stTRInfo.nHeight)^2+stTRInfo.nEleFocus^2) - stTRInfo.nEleFocus)/stRFInfo.nC;

%% Geometry setup

% transducer geometry
aElePosX = linspace(-0.5*(stTRInfo.nNoEle-1),0.5*(stTRInfo.nNoEle-1),stTRInfo.nNoEle)*stTRInfo.nPitch; % lateral
aElePosZ = zeros(1, stTRInfo.nNoEle);

% Beamforming grid
aDth = linspace(stBFInfo.nRadius, stBFInfo.nRadius+stBFInfo.nDth, stBFInfo.nDthSpl); % axial
aXAxis = linspace(aElePosX(1), aElePosX(end), stBFInfo.nLatScline); % lateral
aSclineTheta = aPosRng_fieldII; % radial

% voxel
mDth = repmat(aDth', 1, stBFInfo.nRadialScline);
mTheta = repmat(aSclineTheta, numel(aDth), 1);

mImgZ = mDth .* cosd(mTheta);
mImgY = mDth .* sind(mTheta);

vImgX = zeros(numel(aDth), stTRInfo.nNoEle, numel(aSclineTheta));
vImgY = zeros(numel(aDth), stTRInfo.nNoEle, numel(aSclineTheta));
vImgZ = zeros(numel(aDth), stTRInfo.nNoEle, numel(aSclineTheta));

for scIdx = 1:numel(aSclineTheta)
    vImgX(:,:,scIdx) = repmat(aElePosX, numel(aDth),1);
    vImgZ(:,:,scIdx) = repmat(mImgZ(:,scIdx), 1, stTRInfo.nNoEle);
    vImgY(:,:,scIdx) = repmat(mImgY(:,scIdx), 1, stTRInfo.nNoEle);
end

%% save the paramters
if(bSave)
    stParam.stRFInfo = stRFInfo;
    stParam.stBFInfo = stBFInfo;
    stParam.stTRInfo = stTRInfo;
    stParam.stTxInfo = stTxInfo;
    stParam.vImgX = vImgX;
    stParam.vImgY = vImgY;
    stParam.vImgZ = vImgZ;
    save([dir_ '/Parameters'], 'stParam', '-v7.3');
end

%% beamforming
disp(['>>> Mode : ' mode]);

vBFedData = [];
if(strcmp(mode, 'Conv'))
    
    % CON
    dir_saveVol = [dir_save '/VolumetricBeamforming']; mkdir(dir_saveVol);
    for scIdx = 1:numel(aSclineTheta)
        currentTime = clock;
        disp(['      Radial scan plane: ' num2str(round(100*scIdx/numel(aSclineTheta))) '%...'...
            num2str(currentTime(4)) ':' num2str(currentTime(5),'%.2d')]);
        
        vRcvData = getRFData(dir_data,data_list(scIdx).name,stTxInfo.num_tx_plane);
        % mBFedData -> not normalized
        [mBFedData, mSynReg] = fLateral_PWSTF(vRcvData, stRFInfo, stTRInfo, stBFInfo, stTxInfo, aElePosX, aElePosZ, (aDth-stBFInfo.nRadius)', aXAxis, nDelayOff);
       
        mBFedData = mBFedData./mSynReg;
        vBFedData = cat(3, vBFedData, mBFedData);
        if(bSave)
            save([dir_saveVol '/RadialPlane_' num2str(scIdx,'%.3d')],'mBFedData','-v7.3');
        end
    end
    
else
    % rSAF
    dir_saveLat = [dir_save '/LateralBeamforming']; mkdir(dir_saveLat);
    dir_saveEle = [dir_save '/ElevationalBeamforming']; mkdir(dir_saveEle);
    dir_saveVol = [dir_save '/VolumetricBeamforming']; mkdir(dir_saveVol);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% lateral direction beamforming %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('     Lateral beamforming');
    vBFedDataLat = []; vSynRegLat = [];
    for scIdx = 1:numel(aSclineTheta)
        
        currentTime = clock;
        disp(['      Radial scan plane: ' num2str(round(100*scIdx/numel(aSclineTheta))) '%...'...
            num2str(currentTime(4)) ':' num2str(currentTime(5),'%.2d')]);
        
        vRcvData = getRFData(dir_data,data_list(scIdx).name,stTxInfo.num_tx_plane);
        [mBFedData, mSynReg] = fLateral_PWSTF(vRcvData, stRFInfo, stTRInfo, stBFInfo, stTxInfo, aElePosX, aElePosZ, (aDth-stBFInfo.nRadius)', aXAxis, nDelayOff);
        
        mBFedData = mBFedData./mSynReg;
        vBFedDataLat = cat(3, vBFedDataLat, mBFedData);
        vSynRegLat = cat(3, vSynRegLat, mSynReg);
        
        if(bSave)
            save([dir_saveLat '/LateralPlane_' num2str(scIdx,'%.3d')],'mBFedData','-v7.3');
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% elevational direction beamforming %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('     Elevational beamforming');
    
    for ele_idx = 1:stTRInfo.nNoEle
        currentTime = clock;
        disp(['      Element idx: ' num2str(ele_idx) '... '...
            num2str(currentTime(4)) ':' num2str(currentTime(5),'%.2d')]);
        
        mBFedDataLat = squeeze(vBFedDataLat(:,ele_idx,:));
        
        [mBFedData, vBFedDataVol, vSynReg, vApod] = fElevational_SA_Vol(mBFedDataLat, stRFInfo, stBFInfo, stTRInfo, stSAInfo, aSclineTheta, mImgY, mImgZ, nDelayOff);
        vBFedData = cat(3, vBFedData, mBFedData);
        
        if(bSave)
            save([dir_saveEle '/Element_' num2str(ele_idx,'%.3d')],'mBFedData','-v7.3');
        end
    end
    
    if(bSave)
        % save by radial scanning plane
        for scIdx = 1:numel(aSclineTheta)
            mBFedData = squeeze(vBFedData(:,scIdx,:));
            save([dir_saveVol '/RadialPlane_' num2str(scIdx,'%.3d')],'mBFedData','-v7.3');
        end
    end
end

%%
disp('>>> all process completed');


