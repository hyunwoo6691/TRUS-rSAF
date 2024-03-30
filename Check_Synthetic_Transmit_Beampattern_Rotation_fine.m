clear; close all; clc;

sCurrentPath = pwd;
addpath(sCurrentPath);
addpath('functions');
%% Load
dir_ = uigetdir('','');
load([dir_ '/Parameters.mat']);

stRFInfo = stSaveInfo.stRFInfo;
stTRInfo = stSaveInfo.stTRInfo;
stTxInfo = stSaveInfo.stTxInfo;
nRadius = stSaveInfo.nRadius;
nFOV_Theta = stSaveInfo.nFOV_Theta;
nDth = stSaveInfo.nDth;
aPosRng = stSaveInfo.aPosRng;
stDSCInfo = stSaveInfo.stDSCInfo;
clear stSaveInfo;

%%
% for beam pattern
nSplDth = 250;
nSplTheta = 300;

% SA parameter
nScline = 65;
nNoSyn = 5;
% nNoSyn = 3;
if(nNoSyn ~= 0)
    file_name = 'SA';
else
    file_name = 'Conv';
end

errors_ = [0.1 0.2 0.5 1 2 5];
% errors_ = 0;
%%
for e_idx = 1:numel(errors_)
    % error
    error_range = errors_(e_idx);
    disp(['Error : ' num2str(error_range) 'deg']);
    
    no_sample = 100;
    depth_ = (1:1:10)*1e-2;
    prob_ = zeros(numel(depth_), nSplTheta);
    
    for s_idx = 1:no_sample
        disp(['>>> Sample ' num2str(s_idx)]); tic;
        %% Select the angles from 2048 transmits
        nFOV = 60;
        no_scline = 129;
        
        aScline_theta = linspace(-0.5*nFOV, 0.5*nFOV, no_scline); % Ground truth transmitted angle
        errors = normrnd(0, error_range, [1, no_scline]); % Gaussian distribution, mean: 0, std: nError_range            
        
        aScline_theta_err = aScline_theta + errors;
        
        a_closest_angles = zeros(1, no_scline);
        tx_idx = zeros(1, no_scline);
        for t = 1:no_scline
            angle_tmp = aScline_theta_err(t);
            idx = find(abs(aPosRng - angle_tmp) == min(abs(aPosRng - angle_tmp)));
            a_closest_angles(t) = aPosRng(idx);
            tx_idx(t) = idx;
        end
        
        aTheta = linspace(-0.5*nFOV, 0.5*nFOV, nSplTheta);
        aDth = linspace(nRadius+1e-3, nRadius+nDth, nSplDth);
        
        dr = abs(aDth(1)-aDth(2));
        da = abs(aTheta(1) - aTheta(2));
        
        aElePosY = nRadius * sind(aScline_theta);
        aElePosZ = nRadius * cosd(aScline_theta);
        
        if(nNoSyn ~= 0)
            aTxIdx = max(nScline-0.5*(nNoSyn-1),1):min(nScline+0.5*(nNoSyn-1),no_scline);
        else % conventional
            aTxIdx = 0.5*(no_scline+1);
        end
        
        nOuterEleY = aElePosY(aTxIdx(end));
        nOuterEleZ = aElePosZ(aTxIdx(end));
        
        vPress = zeros(nSplTheta, nSplDth, numel(aTxIdx));
        m_dist_diff = zeros(nSplDth, numel(aTxIdx));
        m_delay = zeros(nSplDth, numel(aTxIdx));
        for pIdx = 1:numel(aTxIdx)
            sFileName = ['[Beamfield]Tx_Index_' num2str(tx_idx(aTxIdx(pIdx)),'%.4d')];
            disp(['    [' num2str(pIdx) '/' num2str(numel(aTxIdx)) ']: ' sFileName]);
            load([dir_ '/' sFileName]);
            
            nEleCenY = aElePosY(aTxIdx(pIdx));
            nEleCenZ = aElePosZ(aTxIdx(pIdx));
            
            for zIdx = 1:nSplDth
                mTmp = squeeze(stPressure.vPressure(:,zIdx,:)); % x * time
                nFocalPointY = 0;
                nFocalPointZ = aDth(zIdx);
                
                nDistDiff = abs(sqrt((nOuterEleY-nFocalPointY)^2+(nOuterEleZ-nFocalPointZ)^2) - sqrt((nEleCenY-nFocalPointY)^2+(nEleCenZ-nFocalPointZ)^2));
                nDelay_spl = round(nDistDiff/stRFInfo.nC*stRFInfo.nFs);
                
                m_dist_diff(zIdx, pIdx) = nDistDiff;
                m_delay(zIdx, pIdx) = nDelay_spl;
                
                mPress = mTmp(:, (nDelay_spl+1:end));
                aPress = sum(mPress.^2,2);
                vPress(:,zIdx,pIdx) = aPress;
            end
            clear stPressure;
        end
        mPressure = sum(vPress,3);
        mIntensity_syn = mPressure;
        
        for d_idx = 1:numel(depth_)
            depth = depth_(d_idx);
            idx = find(abs(aDth-depth) == min(abs(aDth-depth)));
            b_pattern = mIntensity_syn(:, idx);
            max_idx = find(b_pattern == max(b_pattern));
            prob_(d_idx, max_idx) = prob_(d_idx, max_idx) + 1;
        end
        toc;
    end
    
    %% Save
    dir_save = [dir_ '/Probability_gauss/error_' num2str(error_range)];
    mkdir(dir_save);
    save([dir_save '/' file_name], 'prob_');
    
end