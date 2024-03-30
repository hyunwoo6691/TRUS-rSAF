clc; clear all; close all;

sCurrentPath = pwd;
addpath(sCurrentPath);
addpath('field_ii');
addpath('functions');

field_init(0);
%% Set parameter
stRFInfo.nFc                    = 6.5e6;
stRFInfo.nC                     = 1540;
stRFInfo.nFs                    = 15*stRFInfo.nFc;
% stRFInfo.nAttenFactor       = 0.5/(1e-2*1e6);   % Attenuation factor[dB/(m*Hz)] % 0.5 dB/(cm*MHz)
stRFInfo.nAttenFactor       = 0;   % Attenuation factor[dB/(m*Hz)] % 0.5 dB/(cm*MHz)
stRFInfo.nUnitDis             = stRFInfo.nC / stRFInfo.nFs;
stRFInfo.nLambda            = stRFInfo.nC / stRFInfo.nFc;

stTRInfo.nNoEle                 = 128;
stTRInfo.nPitch                  = 430e-6;
stTRInfo.nWidth                 = 420e-6;
stTRInfo.nHeight                = 5e-3;
stTRInfo.nKerf                    = stTRInfo.nPitch - stTRInfo.nWidth;
% stTRInfo.nElevFnum            = 1;
% stTRInfo.nEleFocus              = stTRInfo.nHeight * stTRInfo.nElevFnum;
stTRInfo.nEleFocus              = 20e-3;

aElePosX = stTRInfo.nPitch * linspace(-0.5*(stTRInfo.nNoEle - 1), 0.5*(stTRInfo.nNoEle - 1), stTRInfo.nNoEle);
aEleLocZ = zeros(1, stTRInfo.nNoEle);

stTxInfo.nCycle                   = 1;
stTxInfo.sApod                   = 'boxcar';

dT                                  = 1 / stRFInfo.nFs;
nTc                                 = 1 / stRFInfo.nFc;

aFocalPoint                      = [0 0 30]*1e-3;

%% Field Dimension
% aX = linspace(-0.5*(stTRInfo.nNoEle-1),0.5*(stTRInfo.nNoEle-1),stTRInfo.nNoEle)*stTRInfo.nPitch;
aX = linspace(-0.5*(stTRInfo.nNoEle-1),0.5*(stTRInfo.nNoEle-1),500)*stTRInfo.nPitch;
% aY = linspace(-10e-3, 10e-3, 800);
aY = 0;
aZ = linspace(10e-3, 120e-3, 200);

nDimX = numel(aX);
nDimY = numel(aY);
nDimZ = numel(aZ);

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
nFarFieldDepth_y = 0.1e-3;     % depth over nFarFieldDepth_x is assumed to be far field on y-z plane

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
% figure(1);
% subplot(2,1,1); show_xdc(pTxAperture); title('Tx Aperture');
% subplot(2,1,2); show_xdc(pRxAperture); title('Rx Aperture');

%% Time delay
nDelay = 0e-6; % sec
nDelay_spl = round(nDelay*stRFInfo.nFs); 
%% Impulse Response
at2 = 0:dT:2.3*nTc;
aTransImpulsResp = sin(2*pi*stRFInfo.nFc*at2);
aTransImpulsResp = aTransImpulsResp.*(hanning(numel(aTransImpulsResp))'); % Transducer's impulse response
% aTransImpulsResp = padarray(aTransImpulsResp.*(hanning(numel(aTransImpulsResp))'),[0,nDelay_spl],0,'pre'); % Transducer's impulse response
xdc_impulse(pTxAperture,aTransImpulsResp); % Tx impulse response

%% Excitation Pulse
at1 = 0:dT:stTxInfo.nCycle*nTc;
% aPulseSeq  = sin(2*pi*stRFInfo.nFc*at1);  % excitation pulse
aPulseSeq  = padarray(sin(2*pi*stRFInfo.nFc*at1),[0,nDelay_spl],0,'pre');  % excitation pulse
xdc_excitation(pTxAperture,aPulseSeq);

at1_pad = (-nDelay:dT:stTxInfo.nCycle*nTc) + nDelay;
at2_pad = (-nDelay:dT:2.3*nTc) + nDelay;

figure(2)
% subplot(2,1,1); plot(at1,aPulseSeq), title('Exitation Pulse')
subplot(2,1,2); plot(at2,aTransImpulsResp), title('Transducer`s Impulse response'); grid on;
subplot(2,1,1); plot(at1_pad,aPulseSeq), title('Exitation Pulse'); grid on;
% subplot(2,1,2); plot(at2_pad,aTransImpulsResp), title('Transducer`s Impulse response'); grid on;

%% Take back source information
mSrcPos = zeros(4, nSubDivNum_x*nSubDivNum_y*stTRInfo.nNoEle); % (x,y,z,eleIdx)x(SrcNum)
nSrcNum = nSubDivNum_x*nSubDivNum_y;
aTxApertureInfo = xdc_get(pTxAperture,'rect');

for eIdx = 1:stTRInfo.nNoEle
    for sIdx = 1:nSrcNum
        nSIdx = (eIdx-1)*nSrcNum + sIdx; % index of each source in the matrix 'mSrcPos'
        mSrcPos(1,nSIdx) = aTxApertureInfo(8, (eIdx-1)*nSubDivNum_x*nSubDivNum_y + sIdx); % x poistion of source
        mSrcPos(2,nSIdx) = aTxApertureInfo(9, (eIdx-1)*nSubDivNum_x*nSubDivNum_y + sIdx); % y poistion of source
        mSrcPos(3,nSIdx) = aTxApertureInfo(10, (eIdx-1)*nSubDivNum_x*nSubDivNum_y + sIdx);% z poistion of source
        mSrcPos(4,nSIdx) = eIdx; % element index of source
    end
end

%% apply delay
apod_window = ones(1, stTRInfo.nNoEle); % for plane wave, all element transmit

tx_angle = -15;
dist_tmp = aElePosX * tand(tx_angle);
dist_tmp = dist_tmp + abs(min(dist_tmp));
delay_curve = dist_tmp / stRFInfo.nC;

xdc_apodization(pTxAperture, 0, apod_window);
xdc_focus_times(pTxAperture, 0, delay_curve);
        
%% Pressure calculation
for zIdx = 1:nDimZ
    disp(['>>> z=' num2str(zIdx) '/' num2str(nDimZ)]);
    for yIdx = 1:nDimY
        mCalcPoints = [aX', repmat(aY(yIdx),nDimX,1) , repmat(aZ(zIdx),nDimX,1)]; % (calcPointNum)x(3dim)
        [mP_tx,nStartTime]=calc_hp(pTxAperture, mCalcPoints);
        %         [mP_tx,nStartTime]=calc_hhp(pTxAperture, pRxAperture, mCalcPoints);
        
        %%% Zero-padding so that the first row = 0 sec
        nPadLen = round(nStartTime/dT);
        mP_tx_pad = padarray(mP_tx,[nPadLen 0],'pre');
        mTdim(yIdx,zIdx) = size(mP_tx_pad,1);
        cRawP{yIdx,zIdx} = mP_tx_pad;
    end
end

%%
% Calculate vtPressure(x,y,z,t) from cRawP
nDimT = max(mTdim(:));
aT = 0:dT:dT*(nDimT-1);
vtPressure = zeros(nDimX,nDimY,nDimZ,nDimT);
for zIdx = 1:nDimZ
    for yIdx = 1:nDimY
        vtPressure(:,yIdx,zIdx,1:mTdim(yIdx,zIdx)) = cRawP{yIdx,zIdx}';
    end
end

% Calculate vIntensity(x,y,z) from vtPressure
vIntensityTmp = zeros(nDimX,nDimY,nDimZ);
for zIdx = 1:nDimZ
    for yIdx = 1:nDimY
        vIntensityTmp(:,yIdx,zIdx) = sum(squeeze(vtPressure(:,yIdx,zIdx,:)).^2,2)/mTdim(yIdx,zIdx); % time average
    end
end
%% Imaging
nNaturalFocalDth = stTRInfo.nHeight^2/(4*stRFInfo.nLambda);

figure(10);
mIntensity = squeeze(vIntensityTmp(:,1,:));
imagesc(aX*1e3,aZ*1e3,db(mIntensity'/max(mIntensity(:))));
xlabel('Depth [mm]'); ylabel('Elevational [mm]'); colormap default; colorbar; caxis([-80 0]);
axis tight; axis equal;
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf,'Position',[141 66 551 534]);
