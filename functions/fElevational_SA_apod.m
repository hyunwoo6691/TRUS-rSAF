function [mBFedData, vBFedData, vSynReg, vApod] = fElevational_SA_apod(vRcvData, stRFInfo, stBFInfo, stTRInfo, stSAInfo, aSclineTheta, mPosY, mPosZ, nDelayOff, idx)
[nNoSpl, nCh, nNoTx] = size(vRcvData);
%% Element position
% element_theta = linspace(-0.5*stBFInfo.nFOV, 0.5*stBFInfo.nFOV, stBFInfo.nScline);
element_theta = aSclineTheta;
aElePosY = stBFInfo.nRadius * sind(element_theta);
aElePosZ = stBFInfo.nRadius * cosd(element_theta);

d_theta = abs(element_theta(1)-element_theta(2));
%% Virtual source position
aVSPosY = (stBFInfo.nRadius+stSAInfo.nVSPos)*sind(element_theta);
aVSPosZ = (stBFInfo.nRadius+stSAInfo.nVSPos)*cosd(element_theta);

%% Beamforming
% vBFedData = zeros(stBFInfo.nDthSpl, nNoTx, stBFInfo.nScline);
mBFedData = zeros(stBFInfo.nDthSpl, stBFInfo.nScline);
% vSynReg = zeros(stBFInfo.nDthSpl, stBFInfo.nScline, nNoTx);
% vApod = zeros(stBFInfo.nDthSpl, nNoTx, stBFInfo.nScline);

vBFedData = [];
vSynReg = [];
vApod = [];

for s_idx = 1:stBFInfo.nScline
    disp(['scanline: ' num2str(s_idx)]);
    %tic;
    aSumTmp = zeros(stBFInfo.nDthSpl, 1);
    aNoSyn = zeros(stBFInfo.nDthSpl,1);
    
    % aperture size
    aAptSize = (sqrt((aElePosY(s_idx) - mPosY(:,s_idx)).^2 + (aElePosZ(s_idx) - mPosZ(:,s_idx)).^2) - stBFInfo.nRadius )/ stBFInfo.nFnum + d_theta;
        
    beam_use = max(s_idx - 0.5*(stSAInfo.nNoSyn-1),1):min(s_idx + 0.5*(stSAInfo.nNoSyn-1), nNoTx);
    for b_idx = beam_use % include center beam
         %%%%    Tx delay
        aLogicFront = sqrt(mPosY(:,s_idx).^2+mPosZ(:,s_idx).^2) < sqrt(aVSPosY(b_idx)^2+aVSPosZ(b_idx)^2);
        aLogicRear = sqrt(mPosY(:,s_idx).^2+mPosZ(:,s_idx).^2) >= sqrt(aVSPosY(b_idx)^2+aVSPosZ(b_idx)^2);
        
        aTDFront = sqrt((aElePosY(b_idx) - aVSPosY(b_idx)).^2 + (aElePosZ(b_idx) - aVSPosZ(b_idx)).^2) - ...
            sqrt((aVSPosY(b_idx) - mPosY(:, s_idx)).^2 + (aVSPosZ(b_idx) - mPosZ(:, s_idx)).^2);
        aTDRear = sqrt((aElePosY(b_idx) - aVSPosY(b_idx)).^2 + (aElePosZ(b_idx) - aVSPosZ(b_idx)).^2) + ...
            sqrt((aVSPosY(b_idx) - mPosY(:, s_idx)).^2 + (aVSPosZ(b_idx) - mPosZ(:, s_idx)).^2);
        
        aTD = aTDFront .* aLogicFront + aTDRear .* aLogicRear;
        
        %% Synthetic region
        nBeamCenY = aElePosY(b_idx);
        nBeamCenZ = aElePosZ(b_idx);
        
        nVSPosY = aVSPosY(b_idx);
        nVSPosZ = aVSPosZ(b_idx);
        
        ele_edgeR_Y = nBeamCenY + 0.5*stTRInfo.nHeight*cosd(element_theta(b_idx));
        ele_edgeR_Z = nBeamCenZ - 0.5*stTRInfo.nHeight*sind(element_theta(b_idx));
        
        ele_edgeL_Y = nBeamCenY - 0.5*stTRInfo.nHeight*cosd(element_theta(b_idx));
        ele_edgeL_Z = nBeamCenZ + 0.5*stTRInfo.nHeight*sind(element_theta(b_idx));
        
        left_beam_vec = [(nVSPosY-ele_edgeL_Y), (nVSPosZ-ele_edgeL_Z)];
        right_beam_vec = [(nVSPosY-ele_edgeR_Y), (nVSPosZ-ele_edgeR_Z)];
        
        % left beam edge
        y_left = ((left_beam_vec(1)/left_beam_vec(2))*(mPosZ - ele_edgeL_Z) + ele_edgeL_Y);
        
        % right beam edge
        y_right = ((right_beam_vec(1)/right_beam_vec(2))*(mPosZ - ele_edgeR_Z) + ele_edgeR_Y);
        
        mask_ = ((mPosY >= y_left) & (mPosY <= y_right)) | ((mPosY <= y_left) & (mPosY >= y_right)); % upper beam | lower beam
        
        vSynReg(:, :, b_idx) = mask_;
        aNoSyn = aNoSyn + mask_(:,s_idx);
        aNoSyn(aNoSyn == 0) = 1;
        %%
        %%%%    Apodization
%         dist_to_scline = sqrt((aElePosY(b_idx) - mPosY(:,s_idx)).^2 + (aElePosZ(b_idx) - mPosZ(:,s_idx)).^2);
        dist_to_scline = repmat(sqrt((aElePosY(b_idx) - aElePosY(s_idx)).^2 + (aElePosZ(b_idx) - aElePosZ(s_idx)).^2),stBFInfo.nDthSpl,1);
        switch stBFInfo.sWindow
            case 'boxcar'
                aApod = (dist_to_scline <= 0.5 * aAptSize);
            case 'hanning'
%                 aApod = (dist_to_scline <= 0.5 * aAptSize) .* (0.5 - 0.5*cos(2*pi*dist_to_scline ./ aAptSize));
%                 aApod = (dist_to_scline <= 0.5 * aAptSize) .* (0.5*cos(2*pi*dist_to_scline ./ aAptSize));
                aApod = (dist_to_scline <= 0.5 * aAptSize) .* (0.5 - 0.5*cos(2*pi*(b_idx-beam_use(1)+1)./numel(beam_use)));
        end
        vApod(:,b_idx,s_idx) = aApod;
             
        %%%%    Rx Delay
        aRD = sqrt((mPosY(:,s_idx) - aElePosY(b_idx)).^2 + (mPosZ(:,s_idx) - aElePosZ(b_idx)).^2);
        
        %%%%    Tx + Rx
        aDelay = aTD + aRD;
        
        %%%%    Address
        aAddr = aDelay / stRFInfo.nUnitDis + nDelayOff * stRFInfo.nFs;
        aAddr = max(min(aAddr,size(vRcvData,1)-1),1);
        aLogic = (aAddr>0).*(aAddr<size(vRcvData,1)-1);
        
        aInterpData = interpn(0:1:size(vRcvData,1)-1,vRcvData(:,idx,b_idx),aAddr,'spline');
        
        %%%%    Accumulation (Channel sum)
        aSumTmp = aSumTmp + aInterpData .* aLogic .* mask_(:,s_idx) .* aApod;
        
        vBFedData(:, b_idx, s_idx) = aInterpData;
    end
    
    mBFedData(:,s_idx) = aSumTmp./(aNoSyn);% toc;
end
