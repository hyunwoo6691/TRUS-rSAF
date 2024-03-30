clc; clear all; close all;

sCurrentPath = pwd;
addpath(sCurrentPath);
addpath('Field_II_ver_3_24_mac');
addpath('02.Functions');

field_init(0);
%% Set parameter
stRFInfo.nFc                    = 6.5e6;
stRFInfo.nC                     = 1540;
stRFInfo.nFs                    = 15*stRFInfo.nFc;
% stRFInfo.nAttenFactor       = 0.5/(1e-2*1e6);   % Attenuation factor[dB/(m*Hz)] % 0.5 dB/(cm*MHz)
stRFInfo.nAttenFactor       = 0;   % Attenuation factor[dB/(m*Hz)] % 0.5 dB/(cm*MHz)
stRFInfo.nUnitDis             = stRFInfo.nC / stRFInfo.nFs;
stRFInfo.nLambda            = stRFInfo.nC / stRFInfo.nFc;

stTRInfo.nNoEle                 = 1;
stTRInfo.nPitch                  = 430-6;
stTRInfo.nWidth                 = 420e-6;
stTRInfo.nHeight                = 5e-3;
stTRInfo.nKerf                    = stTRInfo.nPitch - stTRInfo.nWidth;
stTRInfo.nElevFnum            = 1;
% stTRInfo.nEleFocus              = stTRInfo.nHeight * stTRInfo.nElevFnum;
stTRInfo.nEleFocus              = 10e-3;

aEleLocX = stTRInfo.nPitch * linspace(-0.5*(stTRInfo.nNoEle - 1), 0.5*(stTRInfo.nNoEle - 1), stTRInfo.nNoEle);
aEleLocZ = zeros(1, stTRInfo.nNoEle);

stTxInfo.nCycle                   = 1;
stTxInfo.sApod                   = 'boxcar';

dT                                  = 1 / stRFInfo.nFs;
nTc                                 = 1 / stRFInfo.nFc;

aFocalPoint                      = [0 0 20]*1e-3;

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
nFarFieldDepth_y = 0.5e-3;     % depth over nFarFieldDepth_x is assumed to be far field on y-z plane

nMathEleSize_x = sqrt(nFarFieldDepth_x * 4 * stRFInfo.nLambda); % H < sqrt(4*lambda*z)
nMathEleSize_y = sqrt(nFarFieldDepth_y * 4 * stRFInfo.nLambda);

nSubDivNum_x = ceil(stTRInfo.nWidth / nMathEleSize_x);
nSubDivNum_y = ceil(stTRInfo.nHeight/ nMathEleSize_y);

pTxAperture = xdc_focused_array (stTRInfo.nNoEle, stTRInfo.nWidth, stTRInfo.nHeight, stTRInfo.nKerf, stTRInfo.nEleFocus, nSubDivNum_x, nSubDivNum_y, aFocalPoint);
pRxAperture = xdc_focused_array (stTRInfo.nNoEle, stTRInfo.nWidth, stTRInfo.nHeight, stTRInfo.nKerf, stTRInfo.nEleFocus, nSubDivNum_x, nSubDivNum_y, aFocalPoint);

xdc_baffle(pTxAperture,0);
xdc_baffle(pRxAperture,-1);
%% Show aperture
% figure(1);
% subplot(2,1,1); show_xdc(pTxAperture); title('Tx Aperture');
% subplot(2,1,2); show_xdc(pRxAperture); title('Rx Aperture');

%% Impulse Response
at2 = 0:dT:2.3*nTc;
aTransImpulsResp = sin(2*pi*stRFInfo.nFc*at2);
aTransImpulsResp = aTransImpulsResp.*(hanning(numel(aTransImpulsResp))'); % Transducer's impulse response
xdc_impulse(pTxAperture,aTransImpulsResp); % Tx impulse response

%% Excitation Pulse
at1 = 0:dT:stTxInfo.nCycle*nTc;
aPulseSeq  = sin(2*pi*stRFInfo.nFc*at1);  % excitation pulse
xdc_excitation(pTxAperture,aPulseSeq);

figure(2)
subplot(2,1,1); plot(at1,aPulseSeq), title('Exitation Pulse'); grid on;
subplot(2,1,2); plot(at2,aTransImpulsResp), title('Transducer`s Impulse response'); grid on;

%% Array center position
nNoTx = 16;
nOuterPos = -1e-3;
aPosRng = linspace(nOuterPos, -1*nOuterPos, nNoTx);

%% Loop for pressure calculation at each position
nDimY = 200;
nDimZ = 800;

vIntensity = zeros(nDimY,nDimZ,nNoTx);
mY = zeros(nNoTx,nDimY);
for pIdx = 1:nNoTx
    disp(['>>> ' num2str(pIdx) '/' num2str(nNoTx) '...']);
    %% Field Dimension
    aX = 0e-3;
    aY = linspace(aPosRng(pIdx)-10e-3, aPosRng(pIdx)+10e-3, nDimY);
    aZ = linspace(0, 120e-3, nDimZ);
    
    nDimX = numel(aX);
    nDimY = numel(aY);
    nDimZ = numel(aZ);
    
    %% Pressure calculation
    for zIdx = 1:nDimZ
%         dispstat(['>>> z=' num2str(zIdx) '/' num2str(nDimZ)]);
        for xIdx = 1:nDimX
            mCalcPoints = [repmat(aX(xIdx),nDimY,1), aY', repmat(aZ(zIdx),nDimY,1)]; % (calcPointNum)x(3dim)
            [mP_tx,nStartTime]=calc_hp(pTxAperture, mCalcPoints);
%             [mP_tx,nStartTime]=calc_hhp(pTxAperture, pRxAperture, mCalcPoints);
            
            %%% Zero-padding so that the first row = 0 sec
            nPadLen = round(nStartTime/dT);
            mP_tx_pad = padarray(mP_tx,[nPadLen 0],'pre');
            mTdim(xIdx,zIdx) = size(mP_tx_pad,1);
            cRawP{xIdx,zIdx} = mP_tx_pad;
        end
    end
    %%
    % Calculate vtPressure(x,y,z,t) from cRawP
    nDimT = max(mTdim(:));
    aT = 0:dT:dT*(nDimT-1);
    vtPressure = zeros(nDimX,nDimY,nDimZ,nDimT);
    for zIdx = 1:nDimZ
        for xIdx = 1:nDimX
            vtPressure(xIdx,:,zIdx,1:mTdim(xIdx,zIdx)) = cRawP{xIdx,zIdx}';
        end
    end
    
    % Calculate vIntensity(x,y,z) from vtPressure
    vIntensityTmp = zeros(nDimX,nDimY,nDimZ);
    for zIdx = 1:nDimZ
        for xIdx = 1:nDimX
            vIntensityTmp(xIdx,:,zIdx) = sum(squeeze(vtPressure(xIdx,:,zIdx,:)).^2,2)/mTdim(xIdx,zIdx); % time average
        end
    end
    
    mIntensity = squeeze(vIntensityTmp(1,:,:));
    mY(pIdx,:) = aY;
    vIntensity(:,:,pIdx) = mIntensity;
    
    cPressure{pIdx} = squeeze(vtPressure(1,:,:,:));
    clear cRawP;
    clear vtPressure;
    clear vIntensityTmp;
end
%% Single beam pattern
Idx = 1;
aYp = linspace(-10e-3, 10e-3, nDimY);
mIntensity_conv = vIntensity(:,:,(0.5*nNoTx));
% mIntensity = sum(vIntensity,3);
figure(1002);
imagesc(aZ*1e3, aYp*1e3, db(mIntensity_conv/max(mIntensity_conv(:)))); %caxis([-100 0]);
xlabel('Depth [mm]'); ylabel('Elevational [mm]'); colormap jet; 
axis tight; axis equal;
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf,'Position',[904 45 945 434]);

%% Synthetic transmit beam pattern
nSclineCenY = 0;
mIntensity_syn = zeros(nDimY, nDimZ);

for zIdx = 1:nDimZ
    disp(['>>> zIdx : ' num2str(zIdx) '/' num2str(nDimZ)]);
    vTmp = zeros(nDimY, max(mTdim), numel(aPosRng));
    for pIdx = 1:numel(aPosRng)
        vPressTmp = cPressure{pIdx};
        
        aEleCenY = aPosRng(pIdx);
%         nDelay = abs(abs(stTRInfo.nEleFocus + sqrt((aEleCenY-nSclineCenY)^2+(aZ(zIdx) - stTRInfo.nEleFocus)^2) - aZ(zIdx)) - abs(stTRInfo.nEleFocus + sqrt((nOuterPos-nSclineCenY)^2+(aZ(zIdx) - stTRInfo.nEleFocus)^2) - aZ(zIdx)));
        nDelay = abs(abs(sqrt((aEleCenY-nSclineCenY)^2+aZ(zIdx)^2) - aZ(zIdx)) - abs(sqrt((nOuterPos-nSclineCenY)^2+aZ(zIdx)^2) - aZ(zIdx)));
        nDelay_spl = round(nDelay /stRFInfo.nC * stRFInfo.nFs);
        
        mPressTmp = squeeze(vPressTmp(:,zIdx,:));
        mPress = mPressTmp(:, (nDelay_spl+1:end));
        vTmp(:,1:size(mPress,2),pIdx) = mPress;
    end
    
    mPressure = sum(vTmp,3);
%     mIntensity_syn(:,zIdx) = sum(mPressure.^2,2)/mTdim(zIdx);
    mIntensity_syn(:,zIdx) = sum(mPressure.^2,2);
end
%%
CropIdx = find(abs(aZ-stTRInfo.nEleFocus) == min(abs(aZ-stTRInfo.nEleFocus)));
CropIdx = 1;
mIntensity_sa = mIntensity_syn(:,CropIdx:end);
figure(1003);
% imagesc(aZ*1e3, aYp*1e3, db(mIntensity_syn/max(mIntensity_syn(:)))); %caxis([-100 0]);
imagesc(aZ(CropIdx:end)*1e3, aYp*1e3, db(mIntensity_sa/max(mIntensity_sa(:)))); %caxis([-100 0]);
xlabel('Depth [mm]'); ylabel('Elevational [mm]'); colormap jet; 
axis tight; axis equal;
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf,'Position',[904 45 945 434]);

%%
d_idx = 120;
a_intensity_conv = db(mIntensity_conv(:, d_idx)/max(mIntensity_conv(:, d_idx)));
a_intensity_sa = db(mIntensity_sa(:, d_idx)/max(mIntensity_sa(:, d_idx)));

figure;
plot(aYp*1e3, a_intensity_conv, 'LineWidth', 2); hold on;
plot(aYp*1e3, a_intensity_sa, 'LineWidth', 2); hold off; grid on;
