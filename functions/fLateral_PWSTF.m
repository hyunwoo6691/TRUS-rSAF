function [mBFedData, mSynReg] = fLateral_PWSTF(vRcvData, stRFInfo, stTRInfo, stBFInfo, stTxInfo, aElePosX, aElePosZ, aDepth, aXAxis, nDelayOff)

[mPosZ, mPosX] = ndgrid(aDepth, aXAxis);

%% transmit angles
tx_angles = stTxInfo.tx_angles;
%% Beamforming
vBFedData = zeros(stBFInfo.nDthSpl, stBFInfo.nLatScline, stTxInfo.num_tx_plane);
vSynReg = zeros(stBFInfo.nDthSpl, stBFInfo.nLatScline, stTxInfo.num_tx_plane);

for t_idx = 1:numel(tx_angles) % at each angle
    mBFed_tmp = zeros(stBFInfo.nDthSpl, stBFInfo.nLatScline); tic;
    
    angle_ = tx_angles(t_idx);
    % select beam region
    beam_edge_lft = aElePosX(1) + aDepth*tand(angle_);
    beam_edge_rgt = aElePosX(end) + aDepth*tand(angle_);
    
    beam_region = (mPosX >= repmat(beam_edge_lft, 1, size(mPosX,2))) & (mPosX <= repmat(beam_edge_rgt, 1, size(mPosX,2)));
    
    % tx offset
    dist_tmp = aElePosX * tand(angle_);
    dist_tmp = dist_tmp + abs(min(dist_tmp));
    tx_offset = 0.5*max(dist_tmp)*cosd(angle_);
    
    % transmit distance
    tx_dist = mPosX(beam_region) * sind(angle_) + mPosZ(beam_region) * cosd(angle_) + tx_offset;
    
    % receive aperture size
    aAptSize = mPosZ(beam_region) / stBFInfo.nFnum + stTRInfo.nPitch;
    for c_idx = 1:stBFInfo.nCh
        % receive distance
        rx_dist = sqrt((mPosX(beam_region)-aElePosX(c_idx)).^2 + (mPosZ(beam_region)-aElePosZ(c_idx)).^2);
        
        aDelay = tx_dist + rx_dist;
        
        aDistX = abs((mPosX(beam_region) - aElePosX(c_idx)));
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
        
        aAddr = aDelay / stRFInfo.nUnitDis + nDelayOff * stRFInfo.nFs;
        aAddr = max(min(aAddr,size(vRcvData,1)-1),1);
        aLogic = (aAddr>0).*(aAddr<size(vRcvData,1)-1);
        
        aInterpData = interpn(0:1:size(vRcvData,1)-1,vRcvData(:,c_idx,t_idx), aAddr, 'spline');
        mBFed_tmp(beam_region) = mBFed_tmp(beam_region) + aInterpData .* aLogic .* aApod;
    end
    vBFedData(:,:,t_idx) = mBFed_tmp;
    vSynReg(:,:,t_idx) = beam_region;
%     disp(['>>> [Lateral, PWSTF] Tx angle: ' num2str(angle_) 'degree | [' num2str(t_idx) '/' num2str(numel(tx_angles)) ']']); toc
end

mSynReg = sum(vSynReg,3);
% mBFedData = sum(vBFedData,3)./mSynReg;
mBFedData = sum(vBFedData,3);
end

