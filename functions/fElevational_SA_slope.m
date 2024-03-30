function [mBFedData, vBFedData, vSynReg, mSlope] = fElevational_SA(vRcvData, stRFInfo, stBFInfo, stTRInfo, stSAInfo, aSclineTheta, mPosY, mPosZ, nDelayOff, idx)
[nNoSpl, nCh, nNoTx] = size(vRcvData);
%% Element position
% element_theta = linspace(-0.5*stBFInfo.nFOV, 0.5*stBFInfo.nFOV, stTRInfo.nNoEle);
element_theta = linspace(-0.5*stBFInfo.nFOV, 0.5*stBFInfo.nFOV, stBFInfo.nScline);
aElePosY = stBFInfo.nRadius * sind(element_theta);
aElePosZ = stBFInfo.nRadius * cosd(element_theta);

%% Virtual source position
aVSPosY = (stBFInfo.nRadius+stSAInfo.nVSPos)*sind(element_theta);
aVSPosZ = (stBFInfo.nRadius+stSAInfo.nVSPos)*cosd(element_theta);

%% Beamforming
vBFedData = zeros(stBFInfo.nDthSpl, nNoTx, stBFInfo.nScline);
mBFedData = zeros(stBFInfo.nDthSpl, stBFInfo.nScline);
vSynReg = zeros(stBFInfo.nDthSpl, stBFInfo.nScline, nNoTx);

mSlope = zeros(2,stBFInfo.nScline);

for s_idx = 1:stBFInfo.nScline
%     disp(['scanline: ' num2str(s_idx)]);
    %tic;
    aSumTmp = zeros(stBFInfo.nDthSpl, 1);
    aNoSyn = zeros(stBFInfo.nDthSpl,1);
    nNoBeam = 0;
    
    for b_idx = max(s_idx - 0.5*(stSAInfo.nNoSyn),1):min(s_idx + 0.5*(stSAInfo.nNoSyn), nNoTx) % 5 beams
%     for b_idx = max(s_idx - 0.5*stSAInfo.nNoSyn + 1,1):min(s_idx + 0.5*stSAInfo.nNoSyn, nNoTx)
        % Transmit delay
        aLogicFront = sqrt(mPosY(:,s_idx).^2+mPosZ(:,s_idx).^2) < sqrt(aVSPosY(b_idx)^2+aVSPosZ(b_idx)^2);
        aLogicRear = sqrt(mPosY(:,s_idx).^2+mPosZ(:,s_idx).^2) >= sqrt(aVSPosY(b_idx)^2+aVSPosZ(b_idx)^2);
        
        aTDFront = sqrt((aElePosY(b_idx) - aVSPosY(b_idx)).^2 + (aElePosZ(b_idx) - aVSPosZ(b_idx)).^2) - ...
            sqrt((aVSPosY(b_idx) - mPosY(:, s_idx)).^2 + (aVSPosZ(b_idx) - mPosZ(:, s_idx)).^2);
        aTDRear = sqrt((aElePosY(b_idx) - aVSPosY(b_idx)).^2 + (aElePosZ(b_idx) - aVSPosZ(b_idx)).^2) + ...
            sqrt((aVSPosY(b_idx) - mPosY(:, s_idx)).^2 + (aVSPosZ(b_idx) - mPosZ(:, s_idx)).^2);
        
        aTD = aTDFront .* aLogicFront + aTDRear .* aLogicRear;
        %         vTD(:, b_idx, s_idx) = aTD;
        %% Synthetic region
        nBeamCenY = aElePosY(b_idx);
        nBeamCenZ = aElePosZ(b_idx);
        
        nVSPosY = aVSPosY(b_idx);
        nVSPosZ = aVSPosZ(b_idx);
        
        nBeamPosR_Y = nBeamCenY + 0.5*stTRInfo.nHeight*cosd(element_theta(b_idx));
        nBeamPosR_Z = nBeamCenZ - 0.5*stTRInfo.nHeight*sind(element_theta(b_idx));
        
        nBeamPosL_Y = nBeamCenY - 0.5*stTRInfo.nHeight*cosd(element_theta(b_idx));
        nBeamPosL_Z = nBeamCenZ + 0.5*stTRInfo.nHeight*sind(element_theta(b_idx));
        
        nSlopeL_EdgeVS = (nBeamPosL_Z - nVSPosZ)/(nBeamPosL_Y - nVSPosY);
        nSlopeR_EdgeVS = (nBeamPosR_Z - nVSPosZ)/(nBeamPosR_Y - nVSPosY);
        
        aPosition = (mPosZ(:,s_idx)-nVSPosZ-nSlopeL_EdgeVS*(mPosY(:,s_idx)-nVSPosY)) .* (mPosZ(:,s_idx)-nVSPosZ-nSlopeR_EdgeVS*(mPosY(:,s_idx)-nVSPosY));
        
        if (nSlopeL_EdgeVS * nSlopeR_EdgeVS > 0)
            aMsk = aPosition < 0;
        else
            aMsk = aPosition > 0;
        end
        
        mSlope(1,b_idx) = nSlopeL_EdgeVS;
        mSlope(2,b_idx) = nSlopeR_EdgeVS;
        vSynReg(:, s_idx, b_idx) = aMsk;
        aNoSyn = aNoSyn + aMsk;
        %%
        %%%%    Rx Delay
        aRD = sqrt((mPosY(:,s_idx) - aElePosY(b_idx)).^2 + (mPosZ(:,s_idx) - aElePosZ(b_idx)).^2);
        
        %%%%    Tx + Rx
        aDelay = aTD + aRD;
        
        aAddr = aDelay / stRFInfo.nUnitDis + nDelayOff * stRFInfo.nFs;
        aAddr = max(min(aAddr,size(vRcvData,1)-1),1);
        aLogic = (aAddr>0).*(aAddr<size(vRcvData,1)-1);
        
        aInterpData = interpn(0:1:size(vRcvData,1)-1,vRcvData(:,idx,b_idx),aAddr,'spline');
        
        %%%%    Accumulation (Channel sum)
        aSumTmp = aSumTmp + aInterpData .* aLogic .* aMsk;
        
        vBFedData(:, b_idx, s_idx) = aInterpData;
        nNoBeam = nNoBeam +1;
    end
%     mBFedData(:,s_idx) = aSumTmp./sqrt(aNoSyn);% toc;
    mBFedData(:,s_idx) = aSumTmp./(aNoSyn);% toc;
end

