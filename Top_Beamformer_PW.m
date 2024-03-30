clc; clear all; close all;

addpath('02.Functions');
%% Load data
sFolderPath = uigetdir('','Select folder');
load([sFolderPath '/vRcvData.mat']);
load([sFolderPath '/stSaveInfo.mat']);
%% Parameter setting
stRFInfo = stSaveInfo.stRFInfo;
stTRInfo = stSaveInfo.stTRInfo;
stTxInfo = stSaveInfo.stTxInfo;

stBFInfo.sDirection             = 'lateral';              % 'elevational' or 'lateral'

% 'lateral' --> 'Conv' or 'PWSTF';
% 'elevational' --> 'Conv' or 'SA';
stBFInfo.sMode                   = 'PWSTF';

stBFInfo.sWindow               = 'boxcar';
stBFInfo.nScline                  = 500;
stBFInfo.nDth                     = 72e-3;
stBFInfo.nDthSpl                = ceil(stBFInfo.nDth/stRFInfo.nUnitDis*2) ;
stBFInfo.nFnum                  = 1; % receive f-number


if (strcmp(stBFInfo.sDirection, 'elevational'))
    if(stTRInfo.nEleFocus == 0)
        nNaturalFocalDth = stTRInfo.nHeight^2/(4*stRFInfo.nLambda);
        stSAInfo.nVSPos = nNaturalFocalDth;
        stSAInfo.nNoSyn = 16;
    end
    
    stBFInfo.nRadius = stSaveInfo.nRadius;
    stBFInfo.nFOV = stSaveInfo.nFOV_Theta; % [deg]
    stBFInfo.nCh = 1;
else % lateral
    stPWSTF.nNoSyn = size(vRcvData, 3);
    stBFInfo.nCh = stTRInfo.nNoEle;
end

stMID.nTGC_Atten = 0.5;                                              % [dB]

stMID.nDCRType = 'high';
stMID.nDCRTap = 128;                                             % BPF tap #
stMID.nDCRFcut = 1e6;
stMID.nDCRF1 = stRFInfo.nFc - 5e6;                                              % [Hz]            if BANDPASS
stMID.nDCRF2 = stRFInfo.nFc + 5e6;                                              % [Hz]           if BANDPASS

stMID.nDemodFreq = stRFInfo.nFc;
stMID.nDemodTap = 128;


nDelayOff = 0;

bSave = 0;
%% Beamforming grid
switch stBFInfo.sDirection
    case 'elevational'
        aSclineTheta = linspace(-0.5*stBFInfo.nFOV, 0.5*stBFInfo.nFOV, stBFInfo.nScline);
        
        aDth = linspace(stBFInfo.nRadius, stBFInfo.nRadius+stBFInfo.nDth, stBFInfo.nDthSpl);
        
        mDth = repmat(aDth', 1, stBFInfo.nScline);
        mTheta = repmat(aSclineTheta, numel(aDth), 1);
        
        mImgZ = mDth .* cosd(mTheta);
        mImgY = mDth .* sind(mTheta);
    case 'lateral'
        aEleIdx = linspace(-0.5*(stTRInfo.nNoEle-1), 0.5*(stTRInfo.nNoEle-1), stTRInfo.nNoEle);
        aElePosX = aEleIdx * stTRInfo.nPitch;
        aElePosZ = zeros(1, stTRInfo.nNoEle);
        
        aXAxis = linspace(aElePosX(1), aElePosX(end), stBFInfo.nScline);
        aZAxis = linspace(0, stBFInfo.nDthSpl * stRFInfo.nUnitDis / 2, stBFInfo.nDthSpl);
        [mImgZ, mImgX] = ndgrid(aZAxis, aXAxis);
end
%% Beamforming
switch stBFInfo.sDirection
    case 'lateral'
        switch stBFInfo.sMode
            case 'Conv'
                [mBFedData] = fLateral_Conv( vRcvData, stRFInfo, stTRInfo, stBFInfo, mImgX, mImgZ, nDelayOff);
            case 'PWSTF'
                [mBFedData, vBFedData, vSynReg] = fLateral_PWSTF(vRcvData, stRFInfo, stTRInfo, stBFInfo, stTxInfo, stPWSTF, aElePosX, aElePosZ, mImgX, mImgZ, nDelayOff);
        end
    case 'elevational'
        switch stBFInfo.sMode
            case 'Conv'
                [mBFedData] = fElevational_Conv(vRcvData, stRFInfo, stBFInfo, mImgY, mImgZ, nDelayOff);
            case 'SA'
                [mBFedData, vBFedData, vSynReg, mSlope] = fElevational_SA(vRcvData, stRFInfo, stBFInfo, stTRInfo, stSAInfo, aSclineTheta, mImgY, mImgZ, nDelayOff);
        end
end

%% MID
[mDCRData, Fil] = DCR(mBFedData, stMID, stRFInfo);
[mTGCOut, aTGCCurve] = fDTGC(mDCRData, stMID, stRFInfo, stBFInfo, size(mBFedData,1), stRFInfo.nUnitDis);
mBEedData = abs(hilbert(mTGCOut));

%% DSC
switch stBFInfo.sDirection
    case 'elevational'
        da = abs(aSclineTheta(1)-aSclineTheta(2));
        dr = abs(aDth(1)-aDth(2));
        nViewDth = stBFInfo.nDth + (stBFInfo.nRadius*(1-cosd(0.5*stBFInfo.nFOV)));
        nViewWidth = 2* (stBFInfo.nRadius+stBFInfo.nDth)*sind(0.5*stBFInfo.nFOV);
        
        dz = 1e-4;
        dy = 1e-4;
        nHeight = round(nViewDth / dz);
        nWidth = round(nViewWidth / dy);
        
        [aYAxis, aZAxis, mOutputTmp] = ScanConverter_convex(mBEedData, dr, da, stBFInfo.nRadius, nHeight, nWidth, dz, dy);
        
%         [aYAxis, aZAxis, mOutputTmp] = ScanConverter_convex(squeeze(vBFedData(:,40,:)), dr, da, stBFInfo.nRadius, nHeight, nWidth, dz, dy);
        
        aROI = find(mOutputTmp ~= 50);
        aOutlier = find(mOutputTmp == 50);
        mOutput = zeros(size(mOutputTmp));
        mOutput(aROI) = mOutputTmp(aROI);
        mImg = db(mOutput/max(mOutput(:)));
        mImg(aOutlier) = -30;
        
        figure(199);
        imagesc(aYAxis*1e3,(aZAxis-stBFInfo.nRadius)*1e3, mImg); caxis([-80 0]);
        axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]');colormap gray; colorbar;
        set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
        set(gcf,'Position',[818 45 1031 676]);
        
    case 'lateral'
        figure(199);
        imagesc(aXAxis*1e3,aZAxis*1e3, db(mBEedData/max(mBEedData(:)))); caxis([-50 0]);
        axis tight; axis equal; xlabel('Lateral [mm]'); ylabel('Axial [mm]');colormap gray; colorbar;
        set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
        set(gcf,'Position',[818 45 1031 676]);
        set(gcf,'Position',[448 1245 778 865]);
end

%% Synthetic region at specific scanline
if(strcmp(stBFInfo.sMode, 'SA')|| strcmp(stBFInfo.sMode, 'PWSTF'))
    %%%%%%%%%%%%
    nScline = round(0.5*stBFInfo.nScline); % center
    %%%%%%%%%%%%
    
    
    if(strcmp(stBFInfo.sDirection,'elevational'))
        aSynIdx = max(nScline-0.5*stSAInfo.nNoSyn,1):min(nScline+0.5*stSAInfo.nNoSyn-1, stBFInfo.nScline);
        
        mSynReg = sum(vSynReg(:,:,aSynIdx),3);
        
        [aYAxis, aZAxis, mOutputTmp] = ScanConverter_convex(mSynReg, dr, da, stBFInfo.nRadius, nHeight, nWidth, dz, dy);
        
        aROI = find(mOutputTmp ~= 50);
        aOutlier = find(mOutputTmp == 50);
        mOutput = zeros(size(mOutputTmp));
        mOutput(aROI) = mOutputTmp(aROI);
        mImg_Syn = mOutput;
        mImg_Syn(aOutlier) = -20;
        
        figure(200);
        imagesc(aYAxis*1e3,(aZAxis-stBFInfo.nRadius)*1e3, mImg_Syn); caxis([-20 stBFInfo.nScline]); hold on;
        line([stBFInfo.nRadius*sind(aSclineTheta(nScline)) (stBFInfo.nRadius+stBFInfo.nDth)*sind(aSclineTheta(nScline))]*1e3,...
            [stBFInfo.nRadius*cosd(aSclineTheta(nScline))-stBFInfo.nRadius (stBFInfo.nRadius+stBFInfo.nDth)*cosd(aSclineTheta(nScline))-stBFInfo.nRadius]*1e3...
            ,'Color', 'r', 'LineWidth', 2, 'LineStyle', '--'); hold off;
        axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]'); colorbar;
        title(['Overlapped beams when beamforming scanline #' num2str(nScline)]);
        set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
        set(gcf,'Position',[818 45 1031 676]);
    else        
        mImg_Syn = sum(vSynReg,3);
        figure(200);
        imagesc(aXAxis*1e3,aZAxis*1e3, mImg_Syn); caxis([0 stPWSTF.nNoSyn]); hold on;
        line([aXAxis(nScline) aXAxis(nScline)]*1e3, [aZAxis(1) aZAxis(end)]*1e3, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--'); hold off;
        axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]'); colorbar;
        title(['Overlapped beams when beamforming scanline #' num2str(nScline)]);
        set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
        set(gcf,'Position',[818 45 1031 676]);
    end
    
end

%% Save
if(bSave)
    disp('>>> Saving...');
    stSaveInfo.mBFedData = mBFedData;
    stSaveInfo.mBEedData = mBEedData;
    switch stBFInfo.sDirection
        case 'lateral'
            stSaveInfo.mImgX = mImgX;
            stSaveInfo.mImgZ = mImgZ;
        case 'elevational'
            stSaveInfo.mImgY = mImgY;
            stSaveInfo.mImgZ = mImgZ;
            stSaveInfo.aYAxis_DSC = aYAxis;
            stSaveInfo.aZAxis_DSC = aZAxis;
            stSaveInfo.mImg_DSC = mImg;
    end
    stSaveInfo.stRFInfo = stRFInfo;
    stSaveInfo.stBFInfo = stBFInfo;
    stSaveInfo.stTRInfo = stTRInfo;
    stSaveInfo.stTxInfo = stTxInfo;
    if(strcmp(stBFInfo.sMode, 'SA'))
        stSaveInfo.stSaInfo = stSAInfo;
        stSaveInfo.mImg_Syn = mImg_Syn;
    end
    
    sSaveDir = [sFolderPath '/Saved data'];
    mkdir(sSaveDir);
    switch stBFInfo.sMode
        case 'Conv'
            save([sSaveDir '/[Save_' stBFInfo.sDirection ']' stBFInfo.sMode '_Focus' num2str(stTxInfo.nTxFocus*1e3)],'stSaveInfo','-v7.3');
        case 'SA'
            save([sSaveDir '/[Save_' stBFInfo.sDirection ']' stBFInfo.sMode '_VS' num2str(stSAInfo.nVSPos*1e3) '_Syn' num2str(stSAInfo.nNoSyn,'%.3d') ],'stSaveInfo', '-v7.3');
    end
end

disp('>>> All process completed.');
