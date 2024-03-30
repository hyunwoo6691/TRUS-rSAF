% clc; clear all;
% close all;
%%
sCurrentPath = pwd;
addpath(sCurrentPath);
addpath('field_ii');
addpath('functions');

bSave = 0;

field_init(0);
%% Set parameter
stRFInfo.nFc                    = 6.5e6;
stRFInfo.nC                     = 1540;
stRFInfo.nFs                    = 16*stRFInfo.nFc;
% stRFInfo.nAttenFactor       = 0.5/(1e-2*1e6);   % Attenuation factor[dB/(m*Hz)] % 0.5 dB/(cm*MHz)
stRFInfo.nAttenFactor       = 0;   % Attenuation factor[dB/(m*Hz)] % 0.5 dB/(cm*MHz)
stRFInfo.nUnitDis             = stRFInfo.nC / stRFInfo.nFs;
stRFInfo.nLambda            = stRFInfo.nC / stRFInfo.nFc;

stTRInfo.nNoEle                 = 1;
stTRInfo.nPitch                  = 430-6;
stTRInfo.nWidth                 = 420e-6;
stTRInfo.nHeight                = 6e-3;
stTRInfo.nKerf                    = stTRInfo.nPitch - stTRInfo.nWidth;
stTRInfo.nElevFnum            = 1;
% stTRInfo.nEleFocus              = stTRInfo.nHeight * stTRInfo.nElevFnum;

stTRInfo.f_num = 1;
% stTRInfo.nEleFocus = stTRInfo.nHeight * stTRInfo.f_num;
stTRInfo.nEleFocus              = 5e-3;
% stTRInfo.nEleFocus            = 0;

aEleLocX = stTRInfo.nPitch * linspace(-0.5*(stTRInfo.nNoEle - 1), 0.5*(stTRInfo.nNoEle - 1), stTRInfo.nNoEle);
aEleLocZ = zeros(1, stTRInfo.nNoEle);

stTxInfo.nCycle                   = 1;
stTxInfo.sApod                   = 'boxcar';
stTxInfo.nNoTx                  = 1;

aTxPosY = linspace(-10e-3, 10e-3, stTxInfo.nNoTx);

dT                                  = 1 / stRFInfo.nFs;
nTc                                 = 1 / stRFInfo.nFc;

%% Natural focal depth
if(~stTRInfo.nEleFocus)
    nNaturalFocalDth = stTRInfo.nHeight^2/(4*stRFInfo.nLambda);
end
%%
aFocalPoint                      = [0 0 stTRInfo.nEleFocus];

%% Set field
set_field('c',stRFInfo.nC);             % Set speed of sound
set_field('Freq_att',stRFInfo.nAttenFactor);     % frequency dependency attenuation in [dB/(m*Hz)] around the center frequency
set_field('att', stRFInfo.nFc*stRFInfo.nAttenFactor);     % Freqency independent attenuation[dB/m]
set_field('att_f0',stRFInfo.nFc);                % Attenuation center frequency[Hz]
set_field('use_att',1);                 % attenuation on/off
set_field('fs',stRFInfo.nFs);        % set the sampling frequency
set_field('use_rectangles',1);          % use ractangles for apertures

%% TX: Generate Sources
nFarFieldDepth_x = 0.01e-3; % depth over nFarFieldDepth_x is assumed to be far field on x-z plane
% For accurate directivity pattern, the number of source is more than 5
nFarFieldDepth_y = 0.01e-3;     % depth over nFarFieldDepth_x is assumed to be far field on y-z plane

nMathEleSize_x = sqrt(nFarFieldDepth_x * 4 * stRFInfo.nLambda); % H < sqrt(4*lambda*z)
nMathEleSize_y = sqrt(nFarFieldDepth_y * 4 * stRFInfo.nLambda);

nSubDivNum_x = ceil(stTRInfo.nWidth / nMathEleSize_x);
nSubDivNum_y = ceil(stTRInfo.nHeight/ nMathEleSize_y);

if(~stTRInfo.nEleFocus)
    pTxAperture = xdc_linear_array(stTRInfo.nNoEle, stTRInfo.nWidth, stTRInfo.nHeight, stTRInfo.nKerf, nSubDivNum_x, nSubDivNum_y, aFocalPoint);
    pRxAperture = xdc_linear_array(stTRInfo.nNoEle, stTRInfo.nWidth, stTRInfo.nHeight, stTRInfo.nKerf, nSubDivNum_x, nSubDivNum_y, aFocalPoint);
else
    pTxAperture = xdc_focused_array (stTRInfo.nNoEle, stTRInfo.nWidth, stTRInfo.nHeight, stTRInfo.nKerf, stTRInfo.nEleFocus, nSubDivNum_x, nSubDivNum_y, aFocalPoint);
    pRxAperture = xdc_focused_array (stTRInfo.nNoEle, stTRInfo.nWidth, stTRInfo.nHeight, stTRInfo.nKerf, stTRInfo.nEleFocus, nSubDivNum_x, nSubDivNum_y, aFocalPoint);
end

xdc_baffle(pTxAperture,0);
xdc_baffle(pRxAperture,-1);
%% Show aperture
figure;
subplot(2,1,1); show_xdc(pTxAperture); title('Tx Aperture');
subplot(2,1,2); show_xdc(pRxAperture); title('Rx Aperture');

%% Impulse Response
at2 = 0:dT:2.3*nTc;
aTransImpulsResp = sin(2*pi*stRFInfo.nFc*at2);
aTransImpulsResp = aTransImpulsResp.*(hanning(numel(aTransImpulsResp))'); % Transducer's impulse response
xdc_impulse(pTxAperture,aTransImpulsResp); % Tx impulse response

%% Excitation Pulse
at1 = 0:dT:stTxInfo.nCycle*nTc;
aPulseSeq  = sin(2*pi*stRFInfo.nFc*at1);  % excitation pulse
xdc_excitation(pTxAperture,aPulseSeq);

figure;
subplot(2,1,1); plot(at1,aPulseSeq), title('Exitation Pulse')
subplot(2,1,2); plot(at2,aTransImpulsResp), title('Transducer`s Impulse response')

%% Array center position
nRadius = 10e-3;
nFOV_Theta = 60;
% nFOV_Theta = 132.2720;

%% Field dimension
aX = 0e-3;

nDth = 72e-3;
nSplDth = 250;
nSplTheta = 300;

offset = 0;
% offset = 2e-3;

% rotate_offset = 0; % deg
rotate_offset = -0.4724 * 1; % deg

aTheta = linspace(-0.5*nFOV_Theta, 0.5*nFOV_Theta, nSplTheta) + rotate_offset;
aDth = linspace(nRadius+offset, nRadius+nDth, nSplDth);

mDth = repmat(aDth', 1, nSplTheta);
mTheta = repmat(aTheta, nSplDth, 1);

mPosZ = mDth .* cosd(mTheta) - nRadius*cosd(mTheta);
mPosY = mDth .* sind(mTheta);

nDimX = numel(aX);
nDimY = nSplTheta;
nDimZ = nSplDth;

% for dsc
dr = abs(aDth(1)-aDth(2));
da = abs(aTheta(1) - aTheta(2));

%% Pressure calculation
for zIdx = 1:nDimZ
    dispstat(['>>> z=' num2str(zIdx) '/' num2str(nDimZ)]);
    for xIdx = 1:nDimX
        mCalcPoints = [repmat(aX(xIdx),nDimY,1), mPosY(zIdx,:)', mPosZ(zIdx,:)']; % (calcPointNum)x(3dim)
        [mP_tx,nStartTime]=calc_hp(pTxAperture, mCalcPoints);
%                             [mP_tx,nStartTime]=calc_hhp(pTxAperture, pRxAperture, mCalcPoints);
        
        %%% Zero-padding so that the first row = 0 sec
        nPadLen = round(nStartTime/dT);
        mP_tx_pad = padarray(mP_tx,[nPadLen 0],'pre');
        mTdim(xIdx,zIdx) = size(mP_tx_pad,1);
        cRawP{xIdx,zIdx} = mP_tx_pad;
    end
end

% Calculate vtPressure(x,y,z,t) from cRawP
nDimT = max(mTdim(:));
aT = 0:dT:dT*(nDimT-1);
vtPressure = zeros(nDimX,nDimY,nDimZ,nDimT);
for zIdx = 1:nDimZ
    for xIdx = 1:nDimX
        vtPressure(xIdx,:,zIdx,1:mTdim(xIdx,zIdx)) = cRawP{xIdx,zIdx}';
    end
end

%Calculate vIntensity(x,y,z) from vtPressure
vIntensityTmp = zeros(nDimX,nDimY,nDimZ);
for zIdx = 1:nDimZ
    for xIdx = 1:nDimX
        vIntensityTmp(xIdx,:,zIdx) = sum(squeeze(vtPressure(xIdx,:,zIdx,:)).^2,2)/mTdim(xIdx,zIdx); % time average
    end
end


%% DSC
nViewDth = nDth + (nRadius*(1-cosd(0.5*nFOV_Theta)));
nViewWidth = 2* (nRadius+nDth)*sind(0.5*nFOV_Theta);

dz = 1e-4;
dy = 1e-4;
nHeight = round(nViewDth / dz);
nWidth = round(nViewWidth / dy);

nOutlierdB = -150;


mIntensity_tmp = squeeze(vIntensityTmp(1,:,:));
% normalize at each radial depth
mIntensity_norm = zeros(size(mIntensity_tmp));
for z_idx = 1:nSplDth
    mIntensity_norm(:,z_idx) = mIntensity_tmp(:,z_idx)/max(mIntensity_tmp(:,z_idx));
end

mIntensity = mIntensity_tmp;

[axis_y, axis_z, mOutputTmp] = ScanConverter_convex(mIntensity', dr, da, nRadius, nHeight, nWidth, dz, dy);

aROI = find(mOutputTmp ~= 50);
aOutlier = find(mOutputTmp == 50);
mOutput = zeros(size(mOutputTmp));
mOutput(aROI) = mOutputTmp(aROI);
mOutput_dB = db(mOutput/max(mOutput(:)));
mOutput_dB(aOutlier) = nOutlierdB;

% for normalized bp
[axis_y, axis_z, mOutputTmp_norm] = ScanConverter_convex(mIntensity_norm', dr, da, nRadius, nHeight, nWidth, dz, dy);

aROI = find(mOutputTmp_norm ~= 50);
aOutlier = find(mOutputTmp_norm == 50);
mOutput = zeros(size(mOutputTmp_norm));
mOutput(aROI) = mOutputTmp_norm(aROI);
mOutput_dB_norm = db(mOutput/max(mOutput(:)));
mOutput_dB_norm(aOutlier) = nOutlierdB;

% figure(3);
figure;
imagesc(axis_y*1e3,(axis_z-nRadius)*1e3, mOutput_dB); hold on;
if(~stTRInfo.nEleFocus)
    scatter(0,nNaturalFocalDth*1e3, 'Marker','x','LineWidth',2, 'MarkerEdgeColor','r','SizeData',50);
    text(0,nNaturalFocalDth*1e3, 'Natural focal depth','VerticalAlignment','top','FontSize',15,'Color','w','FontWeight','bold');
end
axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]');
title(['Elevational aperture size: ' num2str(stTRInfo.nHeight*1e3) 'mm']); hold off;
caxis([-150 0]); colorbar; colormap default;
% xlim([-10 10]);
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf, 'Position', [378 105 726 535]);

% adB_Overlay = [-6, -20, -40, -60, -80];
% subplot(1,2,2);
% [ContourLine, h] = contour(axis_y*1e3,(axis_z - nRadius)*1e3, mOutput_dB, adB_Overlay, 'ShowText', 'on', 'LabelSpacing', 1000, 'LineWidth', 2); hold on;
% if(~stTRInfo.nEleFocus)
%     scatter(0,nNaturalFocalDth*1e3, 'Marker','x','LineWidth',2, 'MarkerEdgeColor','r','SizeData',50);
%     text(0,nNaturalFocalDth*1e3, 'Natural focal depth','VerticalAlignment','top','FontSize',15,'Color','r','FontWeight','bold');
% end
% set(gca, 'Ydir', 'reverse'); grid on; hold off;
% axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]'); title(['Beam pattern contour']);
% set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
% set(gcf,'Position',[818 45 1887 808]);
%%
if(bSave)
    bp.bp_ = mOutput;
    bp.y_axis = axis_y;
    bp.z_axis = axis_z;
    save(['./beamfield/aper_size' num2str(stTRInfo.nHeight*1e3) 'mm'],'bp','-v7.3');
end
%% Normalize at each depth and draw contour
adB_Overlay = [-6 -6];
% figure(3);
figure;
imagesc(axis_y*1e3,(axis_z-nRadius)*1e3, mOutput_dB_norm); hold on;
[ContourLine, h] = contour(axis_y*1e3,(axis_z - nRadius)*1e3, mOutput_dB_norm, adB_Overlay,...
                          'ShowText', 'off', 'LineWidth', 2, 'color','r');
axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]');
title(['Elevational aperture size: ' num2str(stTRInfo.nHeight*1e3) 'mm']); hold off;
caxis([-50 0]); colorbar; colormap default;
% xlim([-10 10]);
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf, 'Position', [378 105 726 535]);

%% Save

if(bSave)
    stDSCInfo.nViewDth = nViewDth;
    stDSCInfo.nViewWidth = nViewWidth;
    stDSCInfo.dz = dz;
    stDSCInfo.dy = dy;
    stDSCInfo.nHeight = nHeight;
    stDSCInfo.nWidth = nWidth;
    
    stSaveInfo.mIntensity = mIntensity;
    stSaveInfo.stRFInfo = stRFInfo;
    stSaveInfo.stTRInfo = stTRInfo;
    stSaveInfo.stTxInfo = stTxInfo;
    stSaveInfo.nRadius = nRadius;
    stSaveInfo.nFOV_Theta = nFOV_Theta;
    stSaveInfo.nDth = nDth;
    stSaveInfo.stDSCInfo = stDSCInfo;
    
    save([sSaveDir '/[Save]FOV' num2str(nFOV_Theta) '_Single_Radius' num2str(nRadius*1e3) 'mm_Focus' num2str(stTRInfo.nEleFocus*1e3) 'mm'],'stSaveInfo','-v7.3');
end

%%
field_end;