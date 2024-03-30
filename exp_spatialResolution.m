% Coordinate setting
nNoPT = 7;
nInitialPosY = 0e-3;
nInitialPosZ = 10e-3;
nSpacing = 10e-3;

% ROI position
aTargetPosY = linspace(nInitialPosY, nInitialPosY, nNoPT);
aTargetPosZ = linspace(nInitialPosZ, nInitialPosZ+(nNoPT-1)*nSpacing, nNoPT);

nROISize_Hor = 60e-3;
nROISize_Ver = 10e-3;

mROIPosY = zeros(nNoPT, 2);
mROIPosZ = zeros(nNoPT,2);

mROIPosY(:,1) = aTargetPosY - 0.5*nROISize_Hor; % left
mROIPosY(:,2) = aTargetPosY + 0.5*nROISize_Hor; % right

mROIPosZ(:,1) = aTargetPosZ - 0.5*nROISize_Ver; % upper
mROIPosZ(:,2) = aTargetPosZ + 0.5*nROISize_Ver; % bottom

dB = [-6 -6];
% dB = [-12 -12];
%%
FWHM_rSAF = zeros(1, nNoPT);

for p_idx = 1:nNoPT
    nLft = mROIPosY(p_idx,1);
    nRgt = mROIPosY(p_idx,2);
    nUp = mROIPosZ(p_idx,1);
    nDn = mROIPosZ(p_idx,2);
    
    lIdx = find(abs(axis_y-nLft)==min(abs(axis_y-nLft)));
    rIdx = find(abs(axis_y-nRgt)==min(abs(axis_y-nRgt)));
    uIdx = find(abs(axis_z-bf_.nRadius-nUp)==min(abs(axis_z-bf_.nRadius-nUp)));
    dIdx = find(abs(axis_z-bf_.nRadius-nDn)==min(abs(axis_z-bf_.nRadius-nDn)));
    
    roi_ = dsc_data(uIdx:dIdx, lIdx:rIdx);
    roi_db = zeros(size(roi_));
    idx_roi = find(roi_ ~= 50);
    idx_reject = find(roi_ == 50);
    roi_db(idx_roi) = db(roi_(idx_roi)/max(roi_(idx_roi)));
    roi_db(idx_reject) = -30;
    
    resol_tmp = measure_spatial_resolution(roi_db,axis_y(lIdx:rIdx)*1e3, (axis_z(uIdx:dIdx)-bf_.nRadius)*1e3,dB, 234, nNoPT, p_idx);
    
    FWHM_rSAF(1, p_idx) = resol_tmp;
end
figure(234);
set(gcf, 'Position', [743 82 575 971]);

