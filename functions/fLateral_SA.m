function [ mBFedData, vBFedData, vSynReg, vTD] = fLateral_SA( vData, stRFInfo, stTRInfo, stBFInfo, stSAInfo, stTxInfo, aTheta, mImgX, mImgZ, nDelayOff)
%% Element position
aElePosX = mImgX(1,:);
aElePosZ = mImgZ(1,:);
%% Virtual source Posation
if (strcmp(stTRInfo.sType,'linear'))
    % Linear
    aVSPosX = mImgX(1,:);
    aVSPosZ = stSAInfo.nVSPos + mImgZ(1,:);
    
    aBeamCenX = mImgX(1,:);
    aBeamCenZ = mImgZ(1,:);
else
    % Convex
    aVSPosX = (stTRInfo.nRadius + stSAInfo.nVSPos) * sind(aTheta);
    aVSPosZ = (stTRInfo.nRadius + stSAInfo.nVSPos) * cosd(aTheta);
        
    aBeamCenX = mImgX(1,:);
    aBeamCenZ = mImgZ(1,:);
end

%% Beamforming
mAccumData = zeros(size(mImgX));

mBFedDataTmp = zeros(size(mImgX));
mBFedData = zeros(size(mImgX));

vSynReg = zeros(size(mImgX,1), size(mImgX,2), size(vData,3));
vTD = zeros(size(mImgX,1), size(mImgX,2), size(vData,3));
vBFedData = zeros(size(mImgX,1), stTRInfo.nNoEle, stBFInfo.nScline);

tic;
% ############################### SA Beamformer for Multi beam #################################
for sIdx = 1: stBFInfo.nScline
    aNoSyn = zeros(size(mImgX,1), 1);
    nBeamNumber = 0; tic;
    aSclineSum = zeros(stBFInfo.nDthSpl,1);
    for bIdx = max(sIdx - 0.5*stSAInfo.nNoSyn, 1): min(sIdx + 0.5*stSAInfo.nNoSyn - 1, stBFInfo.nScline)
        aChSum = zeros(stBFInfo.nDthSpl,1);
        mDataTmp = vData(:,:,bIdx);
        
        % Tx Delay
        aLogicFront = sqrt((mImgX(:,sIdx) - aBeamCenX(bIdx)).^2 + (mImgZ(:,sIdx) - aBeamCenZ(bIdx)).^2) < aVSPosZ(bIdx);
        aLogicRear = sqrt((mImgX(:,sIdx) - aBeamCenX(bIdx)).^2 + (mImgZ(:,sIdx) - aBeamCenZ(bIdx)).^2) >= aVSPosZ(bIdx);
        
        aTDFront = aVSPosZ(bIdx) - sqrt((aVSPosX(bIdx) - mImgX(:, sIdx)).^2 + (aVSPosZ(bIdx) - mImgZ(:, sIdx)).^2);
        aTDRear = aVSPosZ(bIdx) + sqrt((aVSPosX(bIdx) - mImgX(:, sIdx)).^2 + (aVSPosZ(bIdx) - mImgZ(:, sIdx)).^2);
        
        aTD = aTDFront .* aLogicFront + aTDRear .* aLogicRear;
        
        vTD(:, bIdx, sIdx) = aTD;
        % Synthetic region
        if (strcmp(stTRInfo.sType, 'linear'))
            nAperSize = min((stTxInfo.nTxFocus / stTxInfo.nTxFnum), (stTRInfo.nNoEle*stTRInfo.nPitch));
            nTheta = 90 - atand(2*stTxInfo.nTxFocus / nAperSize); % [Deg]
            % Linear
            aEdgeLft = aVSPosX(bIdx) - abs(aVSPosZ(bIdx) - mImgZ(:,sIdx))*tand(nTheta);
            aEdgeRgt = aVSPosX(bIdx) + abs(aVSPosZ(bIdx) - mImgZ(:,sIdx))*tand(nTheta);
%             aSynReg = (mImgX(:,sIdx) >= aEdgeLft ) & (mImgX(:,sIdx) <= aEdgeRgt);
        else
            % Convex
            
        end
        vSynReg(:, sIdx, bIdx) = aSynReg;
        aNoSyn = aNoSyn + aSynReg;
        
        %%%% receive aperture
        if (strcmp(stTRInfo.sType, 'linear'))
            aAptSize = mImgZ(:,sIdx) / stBFInfo.nFnum + stTRInfo.nPitch;
        else
            aAptSize = (sqrt((aElePosX(sIdx) - mImgX(:,sIdx)).^2 + (aElePosZ(sIdx) - mImgZ(:,sIdx)).^2) - stTRInfo.nRadius )/ stBFInfo.nFnum + stTRInfo.nPitch;
        end
        
        for cIdx = max(bIdx - 0.5* stBFInfo.nCh, 1) : min(bIdx + 0.5 * stBFInfo.nCh - 1, stTRInfo.nNoEle)
            %%%%    Distance from imaging point to channel element
            aDistX = abs((mImgX(:,sIdx) - aElePosX(cIdx)));
            
            %%%%    Aperture growth
            switch stBFInfo.sWindow
                case 'hanning'
                    aApod = (aDistX <= 0.5 * aAptSize) .* (0.5 - 0.5*cos(2*pi*aDistX ./ aAptSize));
                case 'tuckey'
                    roll=0.5;
                    aApod = double(aDistX <= (aAptSize/2*(1-roll))) + ...
                        double(aDistX > (aAptSize/2*(1-roll))).*double(aDistX < (aAptSize/2)).* 0.5.*(1+cos(2*pi/roll*(aDistX./aAptSize-roll/2-1/2)));  % tukey window (nBFpoints x nChannel)
                case 'boxcar'
                    aApod = (aDistX <= 0.5 * aAptSize);
                otherwise
                    aApod = ones(1,numel(aAptSize));
            end
            
            %%%%    Rx Delay
            aRD = sqrt((mImgX(:,sIdx) - aElePosX(cIdx)).^2 + (mImgZ(:,sIdx) - aElePosZ(cIdx)).^2);
            
            %%%%    Tx + Rx
            aDelay = aTD + aRD;
            
            aAddr = aDelay / stRFInfo.nUnitDis + nDelayOff * stRFInfo.nFs;
            
            aLogic = (aAddr>0).*(aAddr<size(vData,1)-1);
            
            aInterpData = interpn(0:1:size(mDataTmp,1)-1,mDataTmp(:,cIdx),aAddr,'spline');
            
            %%%%    Accumulation (Channel sum)
            aChSum = aChSum + aInterpData .* aLogic .* aApod .* aSynReg;
        end
        nBeamNumber = nBeamNumber + 1;
        
        % Scanline synthesize
        aSclineSum = aSclineSum + aChSum;
        vBFedData(:, bIdx, sIdx) = aChSum;
    end
    mBFedData(:,sIdx) = aSclineSum ./ sqrt(aNoSyn);
    fprintf(['>>> [Lateral SA] Scan line [' num2str(sIdx) '/' num2str(stBFInfo.nScline) '] is done..!  '...
        num2str(numel(max(sIdx - 0.5*stSAInfo.nNoSyn, 1): min(sIdx + 0.5*stSAInfo.nNoSyn - 1, stBFInfo.nScline))) ' beams synthesized. ']);
    toc;
end

toc;

