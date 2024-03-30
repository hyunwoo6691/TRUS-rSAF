clc; clear; close all;

addpath('functions');
%%
dir_ = uigetdir('./data','');

%%
load([dir_ '/rRotate.mat']);
load([dir_ '/rTR.mat']);

load([dir_ '/posLat_r' num2str(rRotate*1e3,'%.2d') 'mm.mat']);
posLat = pos*1e-3;
load([dir_ '/posAxl_r' num2str(rRotate*1e3,'%.2d') 'mm.mat']);
posAxl = pos*1e-3;
load([dir_ '/posRot_r' num2str(rRotate*1e3,'%.2d') 'mm.mat']);
posRot = pos;
clear pos

load([dir_ '/origin.mat']);
origin_(1:2) = origin_(1:2)*1e-3;

load([dir_ '/scanningAngle.mat']);
%% offset correction
posLat = origin_(1) - posLat;
posAxl = origin_(2) - posAxl;
posRot = origin_(3) - posRot;

%% ground truth trajectory
posGT = [0; -(rRotate-rTR)]+(rRotate-rTR)*[sind(scanningAngle(1,:)); cosd(scanningAngle(1,:))];
%%
figure(1);
subplot(3,1,1);
plot(posGT(1,:),'LineWidth',2,'color','k'); hold on;
plot(posLat,'LineWidth',2,'color','r'); hold off;

subplot(3,1,2);
plot(posGT(2,:),'LineWidth',2,'color','k'); hold on;
plot(posAxl,'LineWidth',2,'color','r'); hold off;

subplot(3,1,3);
plot(scanningAngle(1,:),'LineWidth',2,'color','k'); hold on;
plot(posRot(2:end),'LineWidth',2,'color','r'); hold off;

%% parameter setting for beamforming
bSave = 1;
num_syn = [1 3 5 7 9 11 13 15];
% num_syn = 1;
no_samples = 1;
mode = 'SA';

% load info file
stRFInfo.nFc = 6.5e6;
stRFInfo.nFs = 4*stRFInfo.nFc;
stRFInfo.nC = 1480;
stRFInfo.nLambda = stRFInfo.nC / stRFInfo.nFc;
stRFInfo.nUnitDis = stRFInfo.nC / stRFInfo.nFs;

stTRInfo.nNoEle = 128;
stTRInfo.nPitch = 430e-6;
stTRInfo.nWidth = 420e-6;
stTRInfo.nHeight = 5e-3;
stTRInfo.nKerf = stTRInfo.nPitch - stTRInfo.nWidth;
stTRInfo.nEleFocus = 20e-3;

sDirection = 'elevational'; % lateral or elavational

stTxInfo.type = 'plane_wave';
stTxInfo.cycle = 1;
stTxInfo.num_tx_vol = size(scanningAngle,2); % volume scan - 0.4724
stTxInfo.num_tx_plane = 31; % number of tx at each angle
stTxInfo.max_steering = 15; % [deg]
stTxInfo.apodization = 'boxcar';
stTxInfo.tx_angles = linspace(-stTxInfo.max_steering, stTxInfo.max_steering, stTxInfo.num_tx_plane);

stBFInfo.nDth = 85e-3;
stBFInfo.nDthSpl = ceil(stBFInfo.nDth/stRFInfo.nUnitDis*2);
stBFInfo.nFnum = 1; % receive f-number
stBFInfo.nRadius = rRotate;
stBFInfo.nCh = 1;
stBFInfo.nFOV = abs(scanningAngle(1,1)-scanningAngle(1,end));
stBFInfo.nDTheta = abs(scanningAngle(1,1)-scanningAngle(1,2));
stBFInfo.nScline = size(scanningAngle,2);
stBFInfo.sDirection = 'elevational';
stBFInfo.sWindow = 'boxcar';
stBFInfo.sMode = mode;

nDelayOff = 2*(sqrt((0.5*stTRInfo.nHeight)^2+stTRInfo.nEleFocus^2) - stTRInfo.nEleFocus)/stRFInfo.nC;

fov = stBFInfo.nFOV;
scline = stBFInfo.nScline;
radius = stBFInfo.nRadius;
depth = stBFInfo.nDth;
depth_spl = stBFInfo.nDthSpl;

aSclineTheta = scanningAngle(1,:);
nDelta_theta = abs(aSclineTheta(1)-aSclineTheta(2));

% Beamforming grid
aDth = linspace(radius, radius+depth, depth_spl);

mDth = repmat(aDth', 1, scline);
mTheta = repmat(aSclineTheta, numel(aDth), 1);

mImgZ = mDth .* cosd(mTheta);
mImgY = mDth .* sind(mTheta);

if(bSave)
    stParam.stRFInfo = stRFInfo;
    stParam.stBFInfo = stBFInfo;
    stParam.stTRInfo = stTRInfo;
    stParam.stTxInfo = stTxInfo;
    stParam.mImgY = mImgY;
    stParam.mImgZ = mImgZ;
    save([dir_ '/Parameters'], 'stParam', '-v7.3');
end

%% load data list
dir_data = [dir_ '/raw data'];
data_list = dir(dir_data);
flag = 0;
if(strcmp(data_list(3).name, '.DS_Store')), flag = 1; end
data_list = data_list(3+flag:end);

dir_save = [dir_ '/errors_bf'];

mkdir(dir_save);

%% load raw data
disp('>>> load raw data');
pwIdx = 16;

vRcvData_ = [];
for dIdx = 1:numel(data_list)
    dataTmp = data_list(dIdx).name;
    load([dir_data '/' dataTmp]);
    
    vRcvData_ = cat(3, vRcvData_, vRcvData(:,:,pwIdx));
end
clear vRcvData
%% Beamforming
for synIdx = 1:numel(num_syn)
    disp(['    syn' num2str(num_syn(synIdx))]);
    %%%%%%%%%%%%%%%
    stSAInfo.nNoSyn = num_syn(synIdx);
    %%%%%%%%%%%%%%%
    if(stTRInfo.nEleFocus == 0)
        nNaturalFocalDth = stTRInfo.nHeight^2/(4*stRFInfo.nLambda);
        stSAInfo.nVSPos = nNaturalFocalDth;
    else
        stSAInfo.nVSPos = stTRInfo.nEleFocus;
    end
    
    for ele_idx = 96
        [mBFedData, vBFedData, vSynReg] = fElevational_SA_apod(vRcvData_, stRFInfo, stBFInfo, stTRInfo, stSAInfo, aSclineTheta, mImgY, mImgZ, nDelayOff, ele_idx);
        
        dir_error = [dir_save '/error_0/Sample001/Element_' num2str(ele_idx,'%.3d')];
        if(~exist(dir_error, 'dir'))
            mkdir(dir_error);
        end
        
        stSaveInfo.mBFedData = mBFedData;
        stSaveInfo.vSynReg = vSynReg;
        save([dir_error '/[Save]' stBFInfo.sMode '_VS' num2str(stTRInfo.nEleFocus*1e3) '_Syn' num2str(stSAInfo.nNoSyn,'%.3d') ],'stSaveInfo', '-v7.3');
        clear stSaveInfo
    end
    
end







