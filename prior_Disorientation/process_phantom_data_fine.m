clc; clear; close all;
dir_ = uigetdir('./data','');

scanline_data = dir(dir_);

addpath('./functions');

no_samples = 1;
b_save = 0;
b_display = 1;
sigma_ = 0;
% sigma_ = [0.1 0.2 0.5 1 2 5];
%%
spatial_filter = 'gauss';
%% select element
close all;
element = 77;
% element = 45;

%
no_synthesize = 2;

mode = {'SA', 'Conv'};
nDelayOff = 0;

%
load([dir_ '/stSaveInfo.mat']);
stRFInfo = stSaveInfo.stRFInfo;
stTRInfo = stSaveInfo.stTRInfo;
stTxInfo = stSaveInfo.stTxInfo;

clear stSaveInfo

% mid processing parameter
stMID.nTGC_Atten = 0.5;                % [dB]

stMID.nDCRType = 'high';
stMID.nDCRTap = 128;                   % BPF tap #
stMID.nDCRFcut = 1e6;
stMID.nDCRF1 = stRFInfo.nFc - 5e6;     % [Hz]            if BANDPASS
stMID.nDCRF2 = stRFInfo.nFc + 5e6;     % [Hz]           if BANDPASS

stMID.nDemodFreq = stRFInfo.nFc;
stMID.nDemodTap = 128;

stBFInfo.nDth = 72e-3;
% stBFInfo.nDth = 60e-3;
stBFInfo.nDthSpl = ceil(stBFInfo.nDth/stRFInfo.nUnitDis*2);
stBFInfo.nFnum = 1; % receive f-number

%%%%%%%%%%%%%%%
stSAInfo.nNoSyn = no_synthesize;
%%%%%%%%%%%%%%%

if(stTRInfo.nEleFocus == 0)
    nNaturalFocalDth = stTRInfo.nHeight^2/(4*stRFInfo.nLambda);
    stSAInfo.nVSPos = nNaturalFocalDth;
else
    stSAInfo.nVSPos = stTRInfo.nEleFocus;
    %       stSAInfo.nVSPos = 25e-3;
end

stBFInfo.nRadius = 5e-3;
stBFInfo.nCh = 1;

% %%
fov_original = 60;
scanning_original = 2048;
angle_original = linspace(-0.5*fov_original, 0.5*fov_original, scanning_original);

angle_transmit = fliplr(angle_original(342:1707)); % actually transmitted angle

stBFInfo.nFOV = 40;
scanline_original = linspace(-0.5*fov_original, 0.5*fov_original, 128);
a_transmitted_angles = scanline_original(22:107); % ground truth angles (desired angle)
stBFInfo.nScline = numel(a_transmitted_angles);

for s_idx = 1:no_samples
    tic;
    disp(['>>> ' num2str(s_idx,'%.3d') '/' num2str(no_samples,'%.3d') '...']);
    for e_idx = 1:numel(sigma_)
        disp(['    Error : ' num2str(sigma_(e_idx)) ]);
        
        disp('      Data loading...');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% RF data extract for beamforming %%%%%%%
        aSclineTheta = linspace(-0.5*stBFInfo.nFOV, 0.5*stBFInfo.nFOV, stBFInfo.nScline); % Ground truth transmitted angle
        d_theta = abs(aSclineTheta(1)-aSclineTheta(2));
        
        %         nError_range = 0.1;
        nError_range = sigma_(e_idx);
        aError = normrnd(0, nError_range, [1, stBFInfo.nScline]); % Gaussian distribution, mean: 0, std: nError_range
        
        a_used_angles = a_transmitted_angles + aError; % actually transmitted angle
        
        a_closest_angles = zeros(1, stBFInfo.nScline);
        vRcvData = [];
        for t = 1:stBFInfo.nScline
            %             for t = 1:numel(a_used_angles)
            angle_tmp = a_used_angles(t);
            idx = find(abs(angle_transmit - angle_tmp) == min(abs(angle_transmit - angle_tmp)));
            a_closest_angles(t) = angle_transmit(idx);
            % load raw data
            tmp = ['scanline_' num2str(idx,'%.4d')];
            load([dir_ '/raw_data/' tmp '/raw_data.mat']);
            vRcvData = cat(3, vRcvData, RF);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for m_idx = 1:numel(mode)
            disp(['      Mode : ' mode{m_idx}]);
            
            stBFInfo.sDirection = 'elevational';
            %             stBFInfo.sMode = 'SA';
            stBFInfo.sMode = mode{m_idx};
            stBFInfo.sWindow = 'boxcar';
            
            % Beamforming grid
            switch stBFInfo.sDirection
                case 'elevational'
                    aDth = linspace(stBFInfo.nRadius, stBFInfo.nRadius+stBFInfo.nDth, stBFInfo.nDthSpl);
                    
                    mDth = repmat(aDth', 1, stBFInfo.nScline);
                    mTheta = repmat(aSclineTheta, numel(aDth), 1);
                    
                    mImgZ = mDth .* cosd(mTheta);
                    mImgY = mDth .* sind(mTheta);
                case 'lateral'
                    aEleIdx = linspace(-0.5*(stTRInfo.nNoEle-1), 0.5*(stTRInfo.nNoEle-1), stTRInfo.nNoEle);
                    aElePosX = aEleIdx * stTRInfo.nPitch;
                    aElePosZ = zeros(1, stTRInfo.nNoEle);
                    
                    aXAxis = linspace(aElePosX(1), aElePosX(end), stBFInfo.nScline);
                    aZAxis = linspace(0, stBFInfo.nDthSpl * stRFInfo.nUnitDis / 2, stBFInfo.nDthSpl);
                    [mImgZ, mImgX] = ndgrid(aZAxis, aXAxis);
            end
            
            for ele_idx = element
                switch stBFInfo.sMode
                    case 'Conv'
                        switch stBFInfo.sDirection
                            case 'elevational'
                                [mBFedData] = fElevational_Conv(vRcvData, stRFInfo, stBFInfo, mImgY, mImgZ, nDelayOff, ele_idx);
                            case 'lateral'
                        end
                    case 'SA'
                        switch stBFInfo.sDirection
                            case 'elevational'
                                [mBFedData, vBFedData, vSynReg, mSlope] = fElevational_SA(vRcvData, stRFInfo, stBFInfo, stTRInfo, stSAInfo, aSclineTheta, mImgY, mImgZ, nDelayOff, ele_idx);
                            case 'lateral'
                                if(strcmp(stTRInfo.sType, 'linear'))
                                    [mBFedData, vBFedData, vSynReg] = fLateral_PWSTF(vRcvData, stRFInfo, stTRInfo, stBFInfo, stTxInfo, stPWSTF, mImgX, mImgZ, nDelayOff);
                                else % convex
                                end
                        end
                end
                mBFedData(isnan(mBFedData)) = 0;
                
%                 [mBFedData, fil] = spatial_filtering(mBFedData, sigma_, d_theta, spatial_filter); % 'gauss', 'boxcar', 'none'
                
                % MID
                [mDCRData, Fil] = DCR(mBFedData, stMID, stRFInfo);
                [mTGCOut, aTGCCurve] = fDTGC(mDCRData, stMID, stRFInfo, stBFInfo, size(mBFedData,1), stRFInfo.nUnitDis);
                mBEedData = abs(hilbert(mTGCOut));
                
%                 [mBEedData, fil] = spatial_filtering(mBEedData, sigma_, d_theta, spatial_filter); % 'gauss', 'boxcar', 'none'
%                 [mBEedData, fil] = spatial_filtering(mBEedData, 1, d_theta, spatial_filter); % 'gauss', 'boxcar', 'none'
                 
            end
            
            % plot image
            da = abs(aSclineTheta(1)-aSclineTheta(2));
            dr = abs(aDth(1)-aDth(2));
            nViewDth = stBFInfo.nDth + (stBFInfo.nRadius*(1-cosd(0.5*stBFInfo.nFOV)));
            nViewWidth = 2* (stBFInfo.nRadius+stBFInfo.nDth)*sind(0.5*stBFInfo.nFOV);
            
            dz = 1e-4; dy = 1e-4;
            nHeight = round(nViewDth / dz);
            nWidth = round(nViewWidth / dy);
            
            [aYAxis, aZAxis, mOutputTmp] = ScanConverter_convex(mBEedData, dr, da, stBFInfo.nRadius, nHeight, nWidth, dz, dy);
            
            if(b_display)
                aROI = find(mOutputTmp ~= 50);
                aOutlier = find(mOutputTmp == 50);
                mOutput = zeros(size(mOutputTmp));
                mOutput(aROI) = mOutputTmp(aROI);
                
%                 if(strcmp(mode{m_idx},'SA')), norm_val = max(mOutput(:)); end
                norm_val = max(mOutput(:));
            
                mImg = mOutput;
                mImg = db(mOutput/norm_val);
                mImg(aOutlier) = -30;
                
                figure(element+m_idx);
                imagesc(aYAxis*1e3,(aZAxis-stBFInfo.nRadius)*1e3, mImg); caxis([-60 0]);
                axis tight; axis equal;
                xlabel('Elevational [mm]'); ylabel('Axial [mm]'); title(mode{m_idx});
                colormap gray; colorbar;
                set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
                set(gcf,'Position',[47 90 675 676]);
            
                % synthesized beam region
                scanline_ = 0.5*stBFInfo.nScline;
                %     scanline_ = 96;
                beam_indices = scanline_-0.5*no_synthesize+1 : scanline_+0.5*no_synthesize;
                [aYAxis, aZAxis, m_syn_reg_tmp] = ScanConverter_convex(sum(vSynReg(:,:,beam_indices),3), dr, da, stBFInfo.nRadius, nHeight, nWidth, dz, dy);
                %     [aYAxis, aZAxis, m_syn_reg_tmp] = ScanConverter_convex(sum(vSynReg(:,:,scanline_),3), dr, da, stBFInfo.nRadius, nHeight, nWidth, dz, dy);
                aROI = find(m_syn_reg_tmp ~= 50);
                aOutlier = find(m_syn_reg_tmp == 50);
                m_syn_reg = zeros(size(m_syn_reg_tmp));
                m_syn_reg(aROI) = m_syn_reg_tmp(aROI);
                
                figure(1010);
                imagesc(aYAxis*1e3,(aZAxis-stBFInfo.nRadius)*1e3, m_syn_reg);
                axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]');colormap default; colorbar;
                set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
                set(gcf,'Position',[818 45 1031 676]);
            end
            
            % save image
            if(b_save)
                dir_save = [dir_ '/Errors'];
                if(~exist(dir_save, 'dir'))
                    mkdir(dir_save);
                end
                
                dir_error = [dir_save '/error_' num2str(nError_range) '/Sample' num2str(s_idx,'%.3d') '/Element_' num2str(ele_idx)];
                if(~exist(dir_error, 'dir'))
                    mkdir(dir_error);
                end
                
                output.img_mag = mOutputTmp;
                output.z_axis = aZAxis;
                output.y_axis = aYAxis;
                
                save([dir_error '/img_' mode{m_idx} '_' num2str(element) '_sf_' spatial_filter],'output','-v7.3');
            end
            
        end
    end
    toc;
end