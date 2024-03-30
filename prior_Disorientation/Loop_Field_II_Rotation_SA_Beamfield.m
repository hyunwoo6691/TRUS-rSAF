clc; clear all; close all;
%%
sCurrentPath = pwd;
addpath(sCurrentPath);
addpath('Field_II_ver_3_24_mac');
addpath('02.Functions');

bSave = 1;

aEleFocus = [20e-3];
aNoTx = [2048];

Number = 1;
for a = 1:numel(aEleFocus)
    for b = 1:numel(aNoTx)
        disp(['>>> Elevation focus: ' num2str(aEleFocus(a)*1e3) 'mm, Number of Tx: ' num2str(aNoTx(b)) '......[' num2str(Number) '/' num2str(numel(aEleFocus)*numel(aNoTx)) ']']);
        field_init(0);
        %% Set parameter
        stRFInfo.nFc                    = 6.5e6;
        stRFInfo.nC                     = 1540;
        stRFInfo.nFs                    = 10*stRFInfo.nFc;
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
        stTRInfo.nEleFocus              = aEleFocus(a);
        
        aEleLocX = stTRInfo.nPitch * linspace(-0.5*(stTRInfo.nNoEle - 1), 0.5*(stTRInfo.nNoEle - 1), stTRInfo.nNoEle);
        aEleLocZ = zeros(1, stTRInfo.nNoEle);
        
        stTxInfo.nCycle                   = 1;
        stTxInfo.sApod                   = 'boxcar';
        stTxInfo.nNoTx                  = aNoTx(b);
        
        aTxPosY = linspace(-10e-3, 10e-3, stTxInfo.nNoTx);
        
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
        figure(1);
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
        
        figure(2)
        subplot(2,1,1); plot(at1,aPulseSeq), title('Exitation Pulse')
        subplot(2,1,2); plot(at2,aTransImpulsResp), title('Transducer`s Impulse response')
        
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
        
        %% Array center position
        nRadius = 5e-3;
        nFOV_Theta = 65;
        nNoTx = stTxInfo.nNoTx;
        aPosRng = linspace(-0.5*nFOV_Theta, 0.5*nFOV_Theta, nNoTx);
        
        %% Field dimension
        nDth = 105e-3;
        nSplDth = 250;
        nSplTheta = 300;
        
        aTheta = linspace(-0.5*nFOV_Theta, 0.5*nFOV_Theta, nSplTheta);
        aDth = linspace(nRadius+1e-3, nRadius+nDth, nSplDth);
        
        mDth = repmat(aDth', 1, nSplTheta);
        mTheta = repmat(aTheta, nSplDth, 1);
        
        mPosZ = mDth .* cosd(mTheta);
        mPosY = mDth .* sind(mTheta);
        
        dr = abs(aDth(1)-aDth(2));
        da = abs(aTheta(1) - aTheta(2));
        
        %DSC Info
        nViewDth = nDth + (nRadius*(1-cosd(0.5*nFOV_Theta)));
        nViewWidth = 2* (nRadius+nDth)*sind(0.5*nFOV_Theta);
        
        dz = 1e-4;
        dy = 1e-4;
        nHeight = round(nViewDth / dz);
        nWidth = round(nViewWidth / dy);
        %% Loop for pressure calculation at each position
        vIntensity = zeros(nSplTheta,nSplDth,nNoTx);
        mY = zeros(nNoTx,nSplTheta);
        for pIdx = 1:nNoTx
            disp(['>>> ' num2str(pIdx) '/' num2str(nNoTx)]);tic;
            %% Field Dimension
            aX = 0e-3;
            aTheta = linspace(-0.5*nFOV_Theta, 0.5*nFOV_Theta, nSplTheta) + aPosRng(pIdx);
            
            mTheta = repmat(aTheta, nSplDth, 1);
            
            mPosZ = mDth .* cosd(mTheta) - nRadius * cosd(mTheta);
            mPosY = mDth .* sind(mTheta);
            
            nDimX = numel(aX);
            nDimY = nSplTheta;
            nDimZ = nSplDth;
            
            %% Pressure calculation
            for zIdx = 1:nDimZ
                for xIdx = 1:nDimX
                    mCalcPoints = [repmat(aX(xIdx),nDimY,1), mPosY(zIdx,:)', mPosZ(zIdx,:)']; % (calcPointNum)x(3dim)
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
            vIntensity(:,:,pIdx) = mIntensity;
            toc;
        end
        
        %% Save
        sSaveDir = [sCurrentPath '/03.BeamField'];
        mkdir(sSaveDir);
        
        if(bSave)
            stDSCInfo.nViewDth = nViewDth;
            stDSCInfo.nViewWidth = nViewDth;
            stDSCInfo.dz = dz;
            stDSCInfo.dy = dy;
            stDSCInfo.nHeight = nHeight;
            stDSCInfo.nWidth = nWidth;
            
            stSaveInfo.vIntensity = vIntensity;
            stSaveInfo.stRFInfo = stRFInfo;
            stSaveInfo.stTRInfo = stTRInfo;
            stSaveInfo.stTxInfo = stTxInfo;
            stSaveInfo.nRadius = nRadius;
            stSaveInfo.nFOV_Theta = nFOV_Theta;
            stSaveInfo.nDth = nDth;
            stSaveInfo.stDSCInfo = stDSCInfo;
            
            save([sSaveDir '/[Save]FOV' num2str(nFOV_Theta) '_Tx' num2str(stTxInfo.nNoTx) '_Radius' num2str(nRadius*1e3) 'mm_Focus' num2str(stTRInfo.nEleFocus*1e3) 'mm'],'stSaveInfo','-v7.3');
        end
        Number = Number+1;
        field_end;
    end
end
%%