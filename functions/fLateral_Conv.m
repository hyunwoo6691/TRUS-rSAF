function [mBFedData] = fLateral_Conv( vData, stRFInfo, stTRInfo, stBFInfo, mImgX, mImgZ, nDelayOff)
%% Element Info
aElePosX = mImgX(1,:);
aElePosZ = mImgZ(1,:);
%% Beamforming
disp('>>> Beamforming...');
mBFedData = zeros(size(mImgZ));

for sIdx = 1:stBFInfo.nScline
    aTD = sqrt((mImgX(:,sIdx) - aElePosX(sIdx)).^2 + (mImgZ(:,sIdx) - aElePosZ(sIdx)).^2);
    
    %%%    RX APERTURE SIZE %%%%%%%%%%%%%%%%
    if (strcmp(stTRInfo.sType,'linear'))
        aAptSize = mImgZ(:,sIdx) / stBFInfo.nFnum + stTRInfo.nPitch;
    else
        aAptSize = ((sqrt((aElePosX(sIdx) - mImgX(:,sIdx)).^2 + (aElePosZ(sIdx) - mImgZ(:,sIdx)).^2))/ stBFInfo.nFnum + stTRInfo.nPitch)';
    end

    for cIdx = max(sIdx - 0.5*stBFInfo.nCh +1,1):min(sIdx + 0.5*stBFInfo.nCh, stTRInfo.nNoEle)
        %%%%    Distance from imaging point to channel element
        if(strcmp(stTRInfo.sType, 'linear'))
            aDistX = abs((mImgX(:,sIdx) - aElePosX(cIdx)));
        else % convex
            aDistX = sqrt((mImgX(:,sIdx)- mImgX(:,cIdx)).^2 + (mImgZ(:,sIdx)- mImgZ(:,cIdx)).^2);
        end
        
        %%%%    Aperture growth
        switch stBFInfo.sWindow
            case 'hanning'
                aApod = (aDistX <= 0.5 * aAptSize) .* (0.5 - 0.5*cos(2*pi*aDistX ./ aAptSize));
            case 'tuckey'
                roll=0.5;
                aApod = double(aDistX <= (aAptSize/2*(1-roll))) + ...
                    double(aDistX > (aAptSize/2*(1-roll))).*double(aDistX < (aAptSize/2)).* 0.5.*(1+cos(2*pi/roll*(aDistX./aAptSize-roll/2-1/2)));  % tukey window (nBFpoints x nChannel)
            case 'boxcar'
                aApod = (aDistX <= 0.5 * aAptSize');
            otherwise
                aApod = ones(numel(aAptSize),1);
        end
        
        %%%%    Rx Delay
        aRD = sqrt((mImgX(:,sIdx) - aElePosX(cIdx)).^2 + (mImgZ(:,sIdx) - aElePosZ(cIdx)).^2);
        
        %%%%    Delay summation
        aDelay = aTD + aRD;
        
        %%%%    Address & thresholding        
        aAddr = aDelay / stRFInfo.nUnitDis + nDelayOff* stRFInfo.nFs;
        
        aLogic = (aAddr>0).*(aAddr<size(vData,1)-1);
        
        aInterpData = interpn(0:1:size(vData,1)-1,vData(:,cIdx,sIdx),aAddr,'spline');
        %%%%    Accumulation
        mBFedData(:,sIdx) = mBFedData(:,sIdx) + aInterpData .* aLogic .* aApod;
    end
    fprintf('>>> [Lateral Conv] Scanline : %.3d / %d  done..!\n', sIdx, stBFInfo.nScline);
end

end

