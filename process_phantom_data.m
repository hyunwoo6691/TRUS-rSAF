clc; clear; close all;
dir_ = uigetdir('/01.Data','');

scanline_data = dir(dir_);

addpath('./02.Functions');

b_save = 0;
sigma_ = 5;
%%
disp('>>> Data loading');
vRcvData = [];
% for k = 1:192
for k = 1:128
    tmp = ['scanline_' num2str(k,'%.3d')];
    load([dir_ '/' tmp '/raw_data.mat']);
    vRcvData = cat(3, vRcvData, RF);
end

clear RF k scanline_data tmp

%% select element
element = 64;

%
no_synthesize = 4;

mode = {'SA', 'Conv'};
nDelayOff = 0;

%
load([dir_ '/../stSaveInfo.mat']);
stRFInfo = stSaveInfo.stRFInfo;
stTRInfo = stSaveInfo.stTRInfo;
stTxInfo = stSaveInfo.stTxInfo;

clear stSaveInfo

% mid processing parameter
stMID.nTGC_Atten = 0.5;                % [dB]

stMID.nDCRType = 'high';
stMID.nDCRTap = 128;                   % BPF tap #
stMID.nDCRFcut = 1e6;
stMID.nDCRF1 = stRFInfo.nFc - 5e6;     % [Hz]            if BANDPASS
stMID.nDCRF2 = stRFInfo.nFc + 5e6;     % [Hz]           if BANDPASS

stMID.nDemodFreq = stRFInfo.nFc;
stMID.nDemodTap = 128;

stBFInfo.nDth = 72e-3;
% stBFInfo.nDth = 60e-3;
stBFInfo.nDthSpl = ceil(stBFInfo.nDth/stRFInfo.nUnitDis*2);
stBFInfo.nFnum = 1; % receive f-number

%%%%%%%%%%%%%%%
stSAInfo.nNoSyn = no_synthesize;
%%%%%%%%%%%%%%%

if(stTRInfo.nEleFocus == 0)
    nNaturalFocalDth = stTRInfo.nHeight^2/(4*stRFInfo.nLambda);
    stSAInfo.nVSPos = nNaturalFocalDth;
else
    stSAInfo.nVSPos = stTRInfo.nEleFocus;
    %       stSAInfo.nVSPos = 25e-3;
end

stBFInfo.nRadius = 5e-3;
stBFInfo.nCh = 1;

stBFInfo.nFOV = 60;
stBFInfo.nScline = 128;
% stBFInfo.nFOV = 90;
% stBFInfo.nScline = 192;

aSclineTheta = linspace(-0.5*stBFInfo.nFOV, 0.5*stBFInfo.nFOV, stBFInfo.nScline); % Ground truth transmitted angle
d_theta = abs(aSclineTheta(1)-aSclineTheta(2));

% aSclineTheta = aSclineTheta + 0.5*nDelta_theta;
% aSclineTheta = aSclineTheta(1:stBFInfo.nScline);

for m_idx = 1:numel(mode)
    disp(['      Mode : ' mode{m_idx}]);
    
    stBFInfo.sDirection = 'elevational';
    %             stBFInfo.sMode = 'SA';
    stBFInfo.sMode = mode{m_idx};
    stBFInfo.sWindow = 'boxcar';
    
    % Beamforming grid
    switch stBFInfo.sDirection
        case 'elevational'
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
    
    for ele_idx = element
        switch stBFInfo.sMode
            case 'Conv'
                switch stBFInfo.sDirection
                    case 'elevational'
                        [mBFedData] = fElevational_Conv(vRcvData, stRFInfo, stBFInfo, mImgY, mImgZ, nDelayOff, ele_idx);
                    case 'lateral'
                end
            case 'SA'
                switch stBFInfo.sDirection
                    case 'elevational'
                        [mBFedData, vBFedData, vSynReg, mSlope] = fElevational_SA(vRcvData, stRFInfo, stBFInfo, stTRInfo, stSAInfo, aSclineTheta, mImgY, mImgZ, nDelayOff, ele_idx);
                    case 'lateral'
                        if(strcmp(stTRInfo.sType, 'linear'))
                            [mBFedData, vBFedData, vSynReg] = fLateral_PWSTF(vRcvData, stRFInfo, stTRInfo, stBFInfo, stTxInfo, stPWSTF, mImgX, mImgZ, nDelayOff);
                        else % convex
                        end
                end
        end
        mBFedData(isnan(mBFedData)) = 0;
        
        [mBFedData, fil] = spatial_filtering(mBFedData,no_synthesize, sigma_, d_theta, 'boxcar'); % 'gauss', 'boxcar', 'none'
        
        % MID
        [mDCRData, Fil] = DCR(mBFedData, stMID, stRFInfo);
        [mTGCOut, aTGCCurve] = fDTGC(mDCRData, stMID, stRFInfo, stBFInfo, size(mBFedData,1), stRFInfo.nUnitDis);
        mBEedData = abs(hilbert(mTGCOut));
    end
    
    % plot image
    da = abs(aSclineTheta(1)-aSclineTheta(2));
    dr = abs(aDth(1)-aDth(2));
    nViewDth = stBFInfo.nDth + (stBFInfo.nRadius*(1-cosd(0.5*stBFInfo.nFOV)));
    nViewWidth = 2* (stBFInfo.nRadius+stBFInfo.nDth)*sind(0.5*stBFInfo.nFOV);
    
    dz = 1e-4; dy = 1e-4;
    nHeight = round(nViewDth / dz);
    nWidth = round(nViewWidth / dy);
    
    [aYAxis, aZAxis, mOutputTmp] = ScanConverter_convex(mBEedData, dr, da, stBFInfo.nRadius, nHeight, nWidth, dz, dy);

    aROI = find(mOutputTmp ~= 50);
    aOutlier = find(mOutputTmp == 50);
    mOutput = zeros(size(mOutputTmp));
    mOutput(aROI) = mOutputTmp(aROI);
    mImg = mOutput;
    mImg = db(mOutput/max(mOutput(:)));
    mImg(aOutlier) = -30;
    
    figure(element+m_idx);
    imagesc(aYAxis*1e3,(aZAxis-stBFInfo.nRadius)*1e3, mImg); caxis([-55 0]);
    axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]');colormap gray; colorbar;
    set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
    set(gcf,'Position',[818 45 1031 676]);
    
    % synthesized beam region
    scanline_ = 0.5*stBFInfo.nScline;
%     scanline_ = 96;
    beam_indices = scanline_-0.5*no_synthesize+1 : scanline_+0.5*no_synthesize;
    [aYAxis, aZAxis, m_syn_reg_tmp] = ScanConverter_convex(sum(vSynReg(:,:,beam_indices),3), dr, da, stBFInfo.nRadius, nHeight, nWidth, dz, dy);
%     [aYAxis, aZAxis, m_syn_reg_tmp] = ScanConverter_convex(sum(vSynReg(:,:,scanline_),3), dr, da, stBFInfo.nRadius, nHeight, nWidth, dz, dy);
    aROI = find(m_syn_reg_tmp ~= 50);
    aOutlier = find(m_syn_reg_tmp == 50);
    m_syn_reg = zeros(size(m_syn_reg_tmp));
    m_syn_reg(aROI) = m_syn_reg_tmp(aROI);
    
    figure(1010);
    imagesc(aYAxis*1e3,(aZAxis-stBFInfo.nRadius)*1e3, m_syn_reg);
    axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]');colormap default; colorbar;
    set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
    set(gcf,'Position',[818 45 1031 676]);
    
    % save image
    if(b_save)
        output.img_mag = mOutputTmp;
        output.z_axis = aZAxis;
        output.y_axis = aYAxis;
        
        save([dir_ '/img_' mode{m_idx} '_' num2str(element)],'output','-v7.3');
    end
    
end