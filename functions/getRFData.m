function vRcvData = getRFData(dir_data,folderName,numAngle)
vRcvData = [];
for t = 1:numAngle
    % load raw data
    load([dir_data '/' folderName '/PW_' num2str(t,'%.3d') '.mat']);
    vRcvData = cat(3, vRcvData, mRcvData);
end

end

