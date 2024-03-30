clc; clear all; close all;

sCurrentPath = pwd;
addpath(sCurrentPath);
addpath('02.Functions');

[sFName, sFPath] = uigetfile('','');
load([sFPath sFName]);
%% 
vIntensity = stSaveInfo.vIntensity;
stRFInfo = stSaveInfo.stRFInfo;
stTRInfo = stSaveInfo.stTRInfo;
stTxInfo = stSaveInfo.stTxInfo;
stDSCInfo = stSaveInfo.stDSCInfo;
nRadius = stSaveInfo.nRadius;
nFOV_Theta = stSaveInfo.nFOV_Theta;

%% DSC for interested frame
%%%%%%%%%%%%
nDth = 120e-3;
%%%%%%%%%%%%

nOutlierdB = -150;
adB = [-6, -12, -20, -40, -60, -80];

nScline = round(0.5*stTxInfo.nNoTx);
nNoSyn = 16;
aBeam = max(nScline-0.5*nNoSyn,1):min(nScline+0.5*nNoSyn-1,stTxInfo.nNoTx);

nNoImg = 8; % number of images per figure;
nNoFigure = ceil(nNoSyn / nNoImg);
nNoRow = 2;

aTheta = linspace(-0.5*nFOV_Theta, 0.5*nFOV_Theta, size(vIntensity,1));
aDth = linspace(nRadius, nRadius+nDth, size(vIntensity,2));
dr = abs(aDth(1)-aDth(2));
da = abs(aTheta(1) - aTheta(2));

for bIdx = 1:nNoSyn
    mSynReg = vIntensity(:,:,aBeam(bIdx));
    
    [aYAxis, aZAxis, mOutputTmp] = ScanConverter_convex(mSynReg', dr, da, nRadius, stDSCInfo.nHeight, stDSCInfo.nWidth, stDSCInfo.dz, stDSCInfo.dy);
    
    aROI = find(mOutputTmp ~= 50);
    aOutlier = find(mOutputTmp == 50);
    mOutput = zeros(size(mOutputTmp));
    mOutput(aROI) = mOutputTmp(aROI);
    mOutput_dB = db(mOutput/max(mOutput(:)));
    mOutput_dB(aOutlier) = nOutlierdB;
    
    %     figure(bIdx);
    FigureIdx = ceil(bIdx/nNoImg);
    nSubplotPos = mod(bIdx,8);
    sText = ['Beam index ' num2str(aBeam(bIdx))];
    
    if(nSubplotPos == 0), nSubplotPos = nNoImg; end
    figure(FigureIdx); subplot(nNoRow, round(nNoImg/nNoRow), nSubplotPos); imagesc(aYAxis*1e3,(aZAxis-nRadius)*1e3, mOutput_dB); hold on;
    line([0 0], [0, nDth]*1e3,'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    text(4, 110, sText, 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'w');
    axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]');
    %     title(['Beam pattern for ' num2str(aBeam(bIdx)) 'th beam']);
    caxis([-150 0]);
    set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
    set(gcf,'Position',[27 17 1889 1088]);
    
    %     figure(10*bIdx);
    figure(1000*FigureIdx); subplot(nNoRow, round(nNoImg/nNoRow), nSubplotPos);
    [ContourLine, h] = contour(aYAxis*1e3,(aZAxis - nRadius)*1e3, mOutput_dB, adB, 'ShowText', 'on', 'LabelSpacing', 1000, 'LineWidth', 2); hold on;
    line([0 0], [0, nDth]*1e3,'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    text(2, 110, sText, 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
    set(gca, 'Ydir', 'reverse'); grid on;
    axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]');
    set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
    set(gcf,'Position',[27 17 1889 1088]);
end
%%
%%%%%%%%%%%%%% overlay
mIntensity_Sum = sum(vIntensity(:,:,aBeam),3);

[aYAxis, aZAxis, mOutputTmp] = ScanConverter_convex(mIntensity_Sum', dr, da, nRadius, stDSCInfo.nHeight, stDSCInfo.nWidth, stDSCInfo.dz, stDSCInfo.dy);

aROI = find(mOutputTmp ~= 50);
aOutlier = find(mOutputTmp == 50);
mOutput = zeros(size(mOutputTmp));
mOutput(aROI) = mOutputTmp(aROI);
mOutput_dB = db(mOutput/max(mOutput(:)));
mOutput_dB(aOutlier) = nOutlierdB;

figure; imagesc(aYAxis*1e3,(aZAxis-nRadius)*1e3, mOutput_dB); hold on;
line([0 0], [0, nDth]*1e3,'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]'); title(['Overlaid beam pattern for scanline ' num2str(nScline)]); caxis([-150 0]);
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf,'Position',[818 45 1031 676]);


adB_Overlay = [-6, -20, -40, -60, -80];
figure;
[ContourLine, h] = contour(aYAxis*1e3,(aZAxis - nRadius)*1e3, mOutput_dB, adB_Overlay, 'ShowText', 'on', 'LabelSpacing', 1000, 'LineWidth', 2);
set(gca, 'Ydir', 'reverse'); grid on;
axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]'); title(['Beam pattern contour']);
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf,'Position',[818 45 1031 676]);
