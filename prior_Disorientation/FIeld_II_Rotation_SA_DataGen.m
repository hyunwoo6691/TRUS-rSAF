% Field II data generate simulation
% Lateral, elevational direction
% Elevation --> element rotating

clc; clear all; close all;

sCurrentPath = pwd;
addpath(sCurrentPath);
addpath('Field_II_ver_3_24_mac');
addpath('02.Functions');

sSaveDir = [sCurrentPath '/01.Data/'];

field_init(0);
%% Set parameter
stRFInfo.nFc                    = 6.5e6;
stRFInfo.nC                     = 1540;
stRFInfo.nFs                    = 4*stRFInfo.nFc;
stRFInfo.nAttenFactor       = 0.5/(1e-2*1e6);   % Attenuation factor[dB/(m*Hz)] % 0.5 dB/(cm*MHz)
stRFInfo.nUnitDis             = stRFInfo.nC / stRFInfo.nFs;
stRFInfo.nLambda            = stRFInfo.nC / stRFInfo.nFc;

stTRInfo.sType                  = 'linear';              % Transducer geometry : linear or convex
stTRInfo.nPitch                  = 430e-6;
stTRInfo.nWidth                 = 420e-6;
stTRInfo.nHeight                = 5e-3;
stTRInfo.nKerf                    = stTRInfo.nPitch - stTRInfo.nWidth;
stTRInfo.nElevFnum            = 1;
% stTRInfo.nEleFocus              = stTRInfo.nHeight * stTRInfo.nElevFnum;
stTRInfo.nEleFocus              = 30e-3;

sDirection                          = 'elevational';         % lateral or elavational

stTxInfo.nCycle                   = 1;
stTxInfo.nNoTx                  = 128;

if(strcmp(sDirection,'lateral'))
    stTRInfo.nNoEle                 = 128;
    
    stTxInfo.nTxFocus              = 30e-3;
    stTxInfo.nTxFnum              = 3;
    stTxInfo.sApod                   = 'boxcar';
else
    stTRInfo.nNoEle                 = 1;
    stTxInfo.nTxFocus               = 10e-3;
end

aElePosX = stTRInfo.nPitch * linspace(-0.5*(stTRInfo.nNoEle - 1), 0.5*(stTRInfo.nNoEle - 1), stTRInfo.nNoEle);
aElePosZ = zeros(1, stTRInfo.nNoEle);

nTc                                 = 1 / stRFInfo.nFc;
nUpsampleFactor             = 5;
nSample                         = 4096;

nSimulFreq                      = stRFInfo.nFs * nUpsampleFactor;
dT                                  = 1 /nSimulFreq;

aFocalPoint                      = [0 0 stTxInfo.nTxFocus];
%% Rotate angle (Elevational only)
if(strcmp(sDirection,'elevational'))
    %%%%%%%%%%%%
    nRadius = 5e-3;
    nFOV_Theta = 60;
	% error range between -0.5*nError_range ~ 0.5*nError_range
	nError_range = 1;
    mean_ = 0;
    sigma_ = 1; 
% 	aError =  2*nError_range* rand(1,stTxInfo.nNoTx) - nError_range;
	aError =  2*nError_range* normrnd(mean_, sigma_, [1,stTxInfo.nNoTx]) - nError_range;
    aPosRng_ground_truth = linspace(0.5* nFOV_Theta, -0.5*nFOV_Theta, stTxInfo.nNoTx);
	aPosRng = aPosRng_ground_truth + aError;
    %%%%%%%%%%%%
end

%% Set field
set_field('c',stRFInfo.nC);             % Set speed of sound
set_field('Freq_att',stRFInfo.nAttenFactor);     % frequency dependency attenuation in [dB/(m*Hz)] around the center frequency
set_field('att', stRFInfo.nFc*stRFInfo.nAttenFactor);     % Freqency independent attenuation[dB/m]
set_field('att_f0',stRFInfo.nFc);                % Attenuation center frequency[Hz]
set_field('use_att',1);                 % attenuation on/off
set_field('fs',nSimulFreq);        % set the sampling frequency
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
subplot(2,1,1); plot(at1,aPulseSeq), title('Exitation Pulse')
subplot(2,1,2); plot(at2,aTransImpulsResp), title('Transducer`s Impulse response')

nLag = round(((numel(aPulseSeq)+numel(aTransImpulsResp)-1)+numel(aTransImpulsResp)-1)/2);
%% Point target position
mScatXYZPos = ...
    [0, 0e-3, 10e-3;
    0, 0e-3, 20e-3;
    0, 0e-3, 30e-3;
    0, 0e-3, 40e-3;
    0, 0e-3, 50e-3;
    0, 0e-3, 60e-3;
    0, 0e-3, 70e-3;
    0, 0e-3, 80e-3;
    0, 0e-3, 90e-3;
    0, 0e-3, 100e-3];
mScatMag = 10*ones(size(mScatXYZPos,1),1);

mScatXYZPos_Rot = mScatXYZPos;

% figure;scatter3(mScatXYZPos_Rot(:,1), mScatXYZPos_Rot(:,2), mScatXYZPos_Rot(:,3),mScatMag,'filled');
figure;scatter(mScatXYZPos_Rot(:,2), mScatXYZPos_Rot(:,3),mScatMag,'filled');
% xlabel('Lateral'); 
xlabel('Elevational'); 
ylabel('Axial');
axis equal; axis tight; grid on;

%% Transmit focus (Lateral direction)
% #f = z/D
if(strcmp(sDirection,'lateral'))
    nAperSize = stTxInfo.nTxFocus/stTxInfo.nTxFnum;
    nNoEleUsed = ceil(nAperSize/stTRInfo.nPitch);
    mApod = zeros(stTRInfo.nNoEle);
    mDelay = zeros(stTRInfo.nNoEle);
    for eIdx = 1:stTRInfo.nNoEle
        aUsedEleIdx = max(eIdx-0.5*nNoEleUsed,1):min(eIdx+0.5*nNoEleUsed-1, stTRInfo.nNoEle);
        aUsedElePosX = aElePosX(max(eIdx-0.5*nNoEleUsed,1):min(eIdx+0.5*nNoEleUsed-1, stTRInfo.nNoEle));
        aDist_EleVS = sqrt((aUsedElePosX-aElePosX(eIdx)).^2 + stTxInfo.nTxFocus^2);
        aDelay = max(aDist_EleVS/stRFInfo.nC) - (aDist_EleVS / stRFInfo.nC);
        
        mApod(eIdx, aUsedEleIdx) = 1;
        mDelay(eIdx, aUsedEleIdx) = aDelay;
    end
end
%% Loop for pressure calculation at each position
vRcvData = zeros(nSample, stTRInfo.nNoEle, stTxInfo.nNoTx);

for pIdx = 1:stTxInfo.nNoTx
    disp(['>>> TX #: ' num2str(pIdx) '/' num2str(stTxInfo.nNoTx)]);
    switch sDirection
        case 'elevational'
            % Field rotation
            nTheta = aPosRng(pIdx);
            aTheta_Zero = atand(mScatXYZPos(:,2)./(nRadius+mScatXYZPos(:,3)));
            aTheta_Prime = aTheta_Zero + nTheta;
            
            aPosZ = ((mScatXYZPos(:,3)+nRadius)./cosd(aTheta_Zero)) .* cosd(aTheta_Prime) - nRadius;
            aPosY = ((mScatXYZPos(:,3)+nRadius)./cosd(aTheta_Zero)) .* sind(aTheta_Prime);
            
            mScatXYZPos_Rot(:,3) = aPosZ;
            mScatXYZPos_Rot(:,2) = aPosY;
        case 'lateral'
            mScatXYZPos_Rot = mScatXYZPos;
            xdc_apodization(pTxAperture, 0, mApod(pIdx,:));
            
            xdc_focus_times(pTxAperture, 0, mDelay(pIdx,:)); % Transmit delay
            
            xdc_center_focus(pRxAperture,[0 0 0]);% Set the origin for the dynamic focusing line.
            xdc_focus_times(pRxAperture,0,zeros(1,stTRInfo.nNoEle));
    end
    
    %%% Calculate the received signals from a collection of scatterers for all elements
    [mRF, nStartTime] = calc_scat_multi(pTxAperture, pRxAperture, mScatXYZPos_Rot, mScatMag);
    
    %%% Zero-padding (so that the first sample is captured at t = 0)
    nPadLen = round(nStartTime/dT);
    mRF_pad = padarray(mRF,[nPadLen 0],'pre');
    
    %%% Truncation or Zero-Padding (so that # of samples = RFInfo.nSample)
    nUpSample = nSample*nUpsampleFactor;
    nFirstSamIdx = 1 + nLag;
    nLastSamIdx = nFirstSamIdx + nUpSample -1;
    
    if(size(mRF_pad,1) > nLastSamIdx)
        mRF_trc = mRF_pad(nFirstSamIdx:nLastSamIdx,:);
    elseif(size(mRF_pad,1) < nLastSamIdx)
        mRF_pad2 = padarray(mRF_pad,[nLastSamIdx-size(mRF_pad,1) 0],'post');
        mRF_trc = mRF_pad2(nFirstSamIdx:nLastSamIdx,:);
    end
    
    %%% Decimate to lower sampling frequency (nFs)
    mRF_deci = mRF_trc(1:nUpsampleFactor:end,:);
    
    %%% Save RF data
    %         disp('Save RF data');
    mRcvData = mRF_deci;
    
    vRcvData(:,:,pIdx) = mRcvData;
end
%%
stSaveInfo.stRFInfo = stRFInfo;
stSaveInfo.stTRInfo = stTRInfo;
stSaveInfo.stTxInfo = stTxInfo;
if(strcmp(sDirection,'elevational'))
    stSaveInfo.nRadius = nRadius;
    stSaveInfo.nFOV_Theta = nFOV_Theta;
    stSaveInfo.aPosRng = aPosRng;
end

if(strcmp(sDirection, 'elevational'))
    sFolderName = ['[' sDirection ']Focus' num2str(stTRInfo.nEleFocus*1e3) 'mm_FOV' num2str(nFOV_Theta) '_R' num2str(nRadius*1e3) 'mm_Tx' num2str(stTxInfo.nNoTx) '_ErrorRange_' num2str(nError_range)];
else % lateral
    sFolderName = ['[' sDirection ',' stTRInfo.sType ']Focus' num2str(stTxInfo.nTxFocus*1e3) 'mm_Fnum' num2str(stTxInfo.nTxFnum) '_Tx' num2str(stTxInfo.nNoTx)];
end

mkdir([sSaveDir sFolderName]);

save([sSaveDir sFolderName '/vRcvData'],'vRcvData');
save([sSaveDir sFolderName '/stSaveInfo'], 'stSaveInfo');

%%
field_end;
