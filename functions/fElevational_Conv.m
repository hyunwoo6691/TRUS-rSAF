function [mBFedData] = fElevational_Conv(vRcvData, stRFInfo, stBFInfo, mPosY, mPosZ, nDelayOff, idx)
[nNoSpl, nCh, nNoTx] = size(vRcvData);

aElePosY = mPosY(1,:);
aElePosZ = mPosZ(1,:);

mBFedData = zeros(stBFInfo.nDthSpl, stBFInfo.nScline);
for sIdx = 1:stBFInfo.nScline
    aTD = sqrt((mPosY(:,sIdx) - aElePosY(sIdx)).^2 + (mPosZ(:,sIdx) - aElePosZ(sIdx)).^2);
    aRD = sqrt((mPosY(:,sIdx) - aElePosY(sIdx)).^2 + (mPosZ(:,sIdx) - aElePosZ(sIdx)).^2);
    aDelay = aTD + aRD;
    aAddr = aDelay / stRFInfo.nUnitDis + nDelayOff*stRFInfo.nFs;
    aAddr = max(min(aAddr,nNoSpl-1),0);
    aLogic = (aAddr>0).*(aAddr<nNoSpl-1);
    aInterpData = interpn(0:1:nNoSpl-1,vRcvData(:,idx,sIdx),aAddr, 'spline');
    mBFedData(:,sIdx) = aInterpData .*aLogic;
    %fprintf('>>> [Elevational Conv] Scanline : %d / %d done. \n',sIdx,stBFInfo.nScline);
end

end
