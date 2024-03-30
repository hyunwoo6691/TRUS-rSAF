clear all; close all; clc;

sCurrentPath = pwd;
addpath(sCurrentPath);
addpath('02.Functions');
%% Load
sPath = uigetdir('','');
load([sPath '/Parameters']);

stRFInfo = stSaveInfo.stRFInfo;
stTRInfo = stSaveInfo.stTRInfo;
stTxInfo = stSaveInfo.stTxInfo;
nRadius = stSaveInfo.nRadius;
nFOV_Theta = stSaveInfo.nFOV_Theta;
nDth = stSaveInfo.nDth;
aPosRng = stSaveInfo.aPosRng;
stDSCInfo = stSaveInfo.stDSCInfo;

nSplDth = 250;
nSplTheta = 300;

aTheta = linspace(-0.5*nFOV_Theta, 0.5*nFOV_Theta, nSplTheta);
aDth = linspace(nRadius+1e-3, nRadius+nDth, nSplDth);

dr = abs(aDth(1)-aDth(2));
da = abs(aTheta(1) - aTheta(2));

aElePosY = nRadius * sind(aPosRng);
aElePosZ = nRadius * cosd(aPosRng);
%% SA parameter setting
nScline = 64;
nNoSyn = 32;

aTxIdx = max(nScline-0.5*nNoSyn,1):min(nScline+0.5*nNoSyn-1,stTxInfo.nNoTx);

nOuterEleY = aElePosY(aTxIdx(end));
nOuterEleZ = aElePosZ(aTxIdx(end));

mIntensity_syn = zeros(nSplTheta, nSplDth);

vPress = zeros(nSplTheta, nSplDth, numel(aTxIdx));
for pIdx = 1:numel(aTxIdx)
    disp(['>>> ' num2str(pIdx) '/' num2str(numel(aTxIdx))]); tic;
    sFileName = ['[Beamfield]Tx_Index_' num2str(aTxIdx(pIdx),'%.3d')];
    load([sPath '/' sFileName]);
    
    nEleCenY = aElePosY(aTxIdx(pIdx));
    nEleCenZ = aElePosZ(aTxIdx(pIdx));
        
    for zIdx = 1:nSplDth
        mTmp = squeeze(stPressure.vPressure(:,zIdx,:));
        nFocalPointY = 0;
        nFocalPointZ = aDth(zIdx);
        
        nDistDiff = abs(sqrt((nOuterEleY-nFocalPointY)^2+(nOuterEleZ-nFocalPointZ)^2) - sqrt((nEleCenY-nFocalPointY)^2+(nEleCenZ-nFocalPointZ)^2));
        nDelay_spl = round(nDistDiff/stRFInfo.nC*stRFInfo.nFs);
        
        mPress = mTmp(:, (nDelay_spl+1:end));
        aPress = sum(mPress.^2,2);
        vPress(:,zIdx,pIdx) = aPress;
    end
    clear stPressure; toc;
end
mPressure = sum(vPress,3);
mIntensity_syn = mPressure;
%% Image output
[aYAxis, aZAxis, mOutputTmp] = ScanConverter_convex(mIntensity_syn', dr, da, nRadius, stDSCInfo.nHeight, stDSCInfo.nWidth, stDSCInfo.dz, stDSCInfo.dy);

aROI = find(mOutputTmp ~= 50);
aOutlier = find(mOutputTmp == 50);
mOutput = zeros(size(mOutputTmp));
mOutput(aROI) = mOutputTmp(aROI);
mOutput(aOutlier) = 1e-2*min(mOutput(aROI));
mOutput_dB = db(mOutput/max(mOutput(:)));

figure(100); imagesc(aYAxis*1e3,(aZAxis-nRadius)*1e3, mOutput_dB); hold on;
line([0 0], [0, nDth]*1e3,'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
if(nNoSyn == 1)
    text(40, 110, 'Conv', 'FontWeight','bold', 'FontSize', 22, 'Color', 'w','HorizontalAlignment', 'center');
else
    text(40, 110, ['Syn ' num2str(nNoSyn)], 'FontWeight','bold', 'FontSize', 22, 'Color', 'w','HorizontalAlignment', 'center');
end
hold off;
axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]');
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf,'Position',[818 45 1031 676]);