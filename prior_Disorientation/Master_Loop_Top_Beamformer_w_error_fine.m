clc; clear; close all;

addpath('02.Functions');

bDisplay = 0;
%%
dir_ = uigetdir('/01.Data', 'Select master folder');
folder_list = dir(dir_);

if(strcmp(folder_list(3).name, '.DS_Store'))
    tmp = folder_list(4:end);
    folder_list = tmp;
    clear tmp;
end

nDelayOff = 0;

bSave = 1;

% no_synthesize = [32 32 16 16 16 16 8 8 8 8 8 4 4 4 4 4 4 4];
% for IUS (10mm, 15mm, 20mm, 25mm)
no_synthesize = [16 8 4 4];

assert(numel(folder_list) == numel(no_synthesize));

no_samples = 100;
error_case = [2 5];
mode = {'Conv', 'SA'};
%% Setup beamforming parameter
for syn_idx = 1:(numel(no_synthesize)-1) % no 25mm, already done
    disp(['>>> Folder name: ' folder_list(syn_idx).name]);
    dir_tmp = [dir_ '/' folder_list(syn_idx).name];
    
    % load info file
    load([dir_tmp '/stSaveInfo.mat']);
    stRFInfo = stSaveInfo.stRFInfo;
    stTRInfo = stSaveInfo.stTRInfo;
    stTxInfo = stSaveInfo.stTxInfo;
    nRadius = stSaveInfo.nRadius;
    nFOV_Theta_fieldII = stSaveInfo.nFOV_Theta;
    aPosRng_fieldII = fliplr(stSaveInfo.aPosRng);
    clear stSaveInfo;
    
    % mid processing parameter
    stMID.nTGC_Atten = 0.5;                % [dB]
    
    stMID.nDCRType = 'high';
    stMID.nDCRTap = 128;                   % BPF tap #
    stMID.nDCRFcut = 1e6;
    stMID.nDCRF1 = stRFInfo.nFc - 5e6;     % [Hz]            if BANDPASS
    stMID.nDCRF2 = stRFInfo.nFc + 5e6;     % [Hz]           if BANDPASS
    
    stMID.nDemodFreq = stRFInfo.nFc;
    stMID.nDemodTap = 128;
    
    % load data file
    load([dir_tmp '/vRcvData.mat']);
    vRcvData_fieldII = vRcvData;
    clear vRcvData;
    
    dir_save = [dir_tmp '/Errors'];
    mkdir(dir_save);
    
    for s_idx = 1:no_samples
        disp(['>>> ' num2str(s_idx,'%.3d') '/' num2str(no_samples) '...']);
        for e_idx = 1:numel(error_case)
            disp(['    Error : ' num2str(error_case(e_idx)) ]);
            for m_idx = 1:numel(mode)
                disp(['      Mode : ' mode{m_idx}]);
                
                stBFInfo.sDirection = 'elevational';
                %             stBFInfo.sMode = 'SA';
                stBFInfo.sMode = mode{m_idx};
                stBFInfo.sWindow = 'boxcar';
                
                stBFInfo.nDth = 120e-3;
                stBFInfo.nDthSpl = ceil(stBFInfo.nDth/stRFInfo.nUnitDis*2);
                stBFInfo.nFnum = 1; % receive f-number
                
                %%%%%%%%%%%%%%%
                %             stSAInfo.nNoSyn = 32;
                stSAInfo.nNoSyn = no_synthesize(syn_idx);
                %%%%%%%%%%%%%%%
                
                if(stTRInfo.nEleFocus == 0)
                    nNaturalFocalDth = stTRInfo.nHeight^2/(4*stRFInfo.nLambda);
                    stSAInfo.nVSPos = nNaturalFocalDth;
                else
                    stSAInfo.nVSPos = stTRInfo.nEleFocus;
                end
                
                stBFInfo.nRadius = nRadius;
                stBFInfo.nCh = 1;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%% RF data extract for beamforming %%%%%%%
                stBFInfo.nFOV = 60;
                stBFInfo.nScline = 128;
                
                aSclineTheta = linspace(-0.5*stBFInfo.nFOV, 0.5*stBFInfo.nFOV, stBFInfo.nScline); % Ground truth transmitted angle
                nDelta_theta = abs(aSclineTheta(1)-aSclineTheta(2));
                
                %         nError_range = 0.1;
                nError_range = error_case(e_idx);
                aError = 2*nError_range * (rand(1, stBFInfo.nScline)-0.5);
                
                a_used_angles = aSclineTheta + aError; % actually transmitted angle
                
                a_closest_angles = zeros(1, stBFInfo.nScline);
                vRcvData = zeros(size(vRcvData_fieldII,1), 1, stBFInfo.nScline);
                for t = 1:stBFInfo.nScline
                    angle_tmp = a_used_angles(t);
                    idx = find(abs(aPosRng_fieldII - angle_tmp) == min(abs(aPosRng_fieldII - angle_tmp)));
                    a_closest_angles(t) = aPosRng_fieldII(idx);
                    vRcvData(:,:,t) = vRcvData_fieldII(:,:,idx);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Beamforming grid
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
                %% Beamforming
                switch stBFInfo.sMode
                    case 'Conv'
                        switch stBFInfo.sDirection
                            case 'elevational'
                                %                             disp('>>> Beamforming function: fElevational_Conv.mat');
                                [mBFedData] = fElevational_Conv(vRcvData, stRFInfo, stBFInfo, mImgY, mImgZ, nDelayOff,1);
                            case 'lateral'
                                %                             disp('>>> Beamforming function: fLateral_Conv.mat');
                                [mBFedData] = fLateral_Conv( vRcvData, stRFInfo, stTRInfo, stBFInfo, mImgX, mImgZ, nDelayOff);
                        end
                    case 'SA'
                        switch stBFInfo.sDirection
                            case 'elevational'
                                %                             disp('>>> Beamforming function: fElevational_SA.mat');
                                [mBFedData, vBFedData, vSynReg, mSlope] = fElevational_SA(vRcvData, stRFInfo, stBFInfo, stTRInfo, stSAInfo, aSclineTheta, mImgY, mImgZ, nDelayOff,1);
                            case 'lateral'
                                %                             disp('>>> Beamforming function: fLateral_SA.mat');
                                if(strcmp(stTRInfo.sType, 'linear'))
                                    [mBFedData, vBFedData, vSynReg, vTD] = fLateral_SA(vRcvData, stRFInfo, stTRInfo, stBFInfo, stSAInfo, stTxInfo, 0, mImgX, mImgZ, nDelayOff);
                                else % convex
                                    
                                end
                        end
                end
                %% MID
                [mDCRData, Fil] = DCR(mBFedData, stMID, stRFInfo);
                %     mDCRData = mBFedData;
                %                             BPFCheck( mBFedData, mDCRData, Fil, stRFInfo, stMID ); % BPF check
                [mTGCOut, aTGCCurve] = fDTGC(mDCRData, stMID, stRFInfo, stBFInfo, size(mBFedData,1), stRFInfo.nUnitDis);
                %     mTGCOut = mDCRData;
                
                %%%%%%%%%%%%%% QDM
                %     [mBEedData, aInphTmp, aInphData, aQDM_LPF, aCos, aSin ]= fQDM( mTGCOut, stRFInfo, stMID, stIQInfo );
                %%%%%%%%%%%%% Hilbert
                mHilbertOut = hilbert(mTGCOut);
                mBEedData = abs(sqrt(mTGCOut.^2+mHilbertOut.^2));
                %                 disp('>>> Hilbert demodulation');
                
                %                 disp('>>> MID done')
                %% DSC
                switch stBFInfo.sDirection
                    case 'elevational'
                        da = abs(aSclineTheta(1)-aSclineTheta(2));
                        dr = abs(aDth(1)-aDth(2));
                        nViewDth = stBFInfo.nDth + (stBFInfo.nRadius*(1-cosd(0.5*stBFInfo.nFOV)));
                        nViewWidth = 2* (stBFInfo.nRadius+stBFInfo.nDth)*sind(0.5*stBFInfo.nFOV);
                        
                        dz = 1e-4;
                        dy = 1e-4;
                        nHeight = round(nViewDth / dz);
                        nWidth = round(nViewWidth / dy);
                        
                        [aYAxis, aZAxis, mOutputTmp] = ScanConverter_convex(mBEedData, dr, da, stBFInfo.nRadius, nHeight, nWidth, dz, dy);
                        
                        aROI = find(mOutputTmp ~= 50);
                        aOutlier = find(mOutputTmp == 50);
                        mOutput = zeros(size(mOutputTmp));
                        mOutput(aROI) = mOutputTmp(aROI);
                        mImg = db(mOutput/max(mOutput(:)));
                        mImg(aOutlier) = -30;
                        
                        if(bDisplay)
                            figure(199);
                            imagesc(aYAxis*1e3,(aZAxis-stBFInfo.nRadius)*1e3, mImg); caxis([-80 0]);
                            axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]');colormap gray; colorbar;
                            set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
                            set(gcf,'Position',[818 45 1031 676]);
                        end
                        
                    case 'lateral'
                        if(bDisplay)
                            figure(199);
                            imagesc(aXAxis*1e3,aZAxis*1e3, db(mBEedData/max(mBEedData(:)))); caxis([-60 0]);
                            axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]');colormap gray; colorbar;
                            set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
                            set(gcf,'Position',[818 45 1031 676]);
                        end
                end
                
                %% Synthetic region at specific scanline
                if(strcmp(stBFInfo.sMode, 'SA'))
                    %%%%%%%%%%%%
                    nScline = round(0.5*stBFInfo.nScline); % center
                    %%%%%%%%%%%%
                    
                    aSynIdx = max(nScline-0.5*stSAInfo.nNoSyn,1):min(nScline+0.5*stSAInfo.nNoSyn-1, stBFInfo.nScline);
                    
                    mSynReg = sum(vSynReg(:,:,aSynIdx),3);
                    
                    if(strcmp(stBFInfo.sDirection,'elevational'))
                        [aYAxis, aZAxis, mOutputTmp] = ScanConverter_convex(mSynReg, dr, da, stBFInfo.nRadius, nHeight, nWidth, dz, dy);
                        
                        aROI = find(mOutputTmp ~= 50);
                        aOutlier = find(mOutputTmp == 50);
                        mOutput = zeros(size(mOutputTmp));
                        mOutput(aROI) = mOutputTmp(aROI);
                        mImg_Syn = mOutput;
                        mImg_Syn(aOutlier) = -20;
                        
                        if(bDisplay)
                            figure(200);
                            imagesc(aYAxis*1e3,(aZAxis-stBFInfo.nRadius)*1e3, mImg_Syn); caxis([-20 stBFInfo.nScline]); hold on;
                            line([stBFInfo.nRadius*sind(aSclineTheta(nScline)) (stBFInfo.nRadius+stBFInfo.nDth)*sind(aSclineTheta(nScline))]*1e3,...
                                [stBFInfo.nRadius*cosd(aSclineTheta(nScline))-stBFInfo.nRadius (stBFInfo.nRadius+stBFInfo.nDth)*cosd(aSclineTheta(nScline))-stBFInfo.nRadius]*1e3...
                                ,'Color', 'r', 'LineWidth', 2, 'LineStyle', '--'); hold off;
                            axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]'); colorbar;
                            title(['Overlapped beams when beamforming scanline #' num2str(nScline)]);
                            set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
                            set(gcf,'Position',[818 45 1031 676]);
                        end
                    else
                        mImg_Syn = mSynReg;
                        if(bDisplay)
                            figure(200);
                            imagesc(aXAxis*1e3,aZAxis*1e3, mImg_Syn); caxis([0 stBFInfo.nScline]); hold on;
                            line([aXAxis(nScline) aXAxis(nScline)]*1e3, [aZAxis(1) aZAxis(end)]*1e3, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--'); hold off;
                            axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]'); colorbar;
                            title(['Overlapped beams when beamforming scanline #' num2str(nScline)]);
                            set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
                            set(gcf,'Position',[818 45 1031 676]);
                        end
                    end
                    
                end
                
                %% Save
                if(bSave)
                    %                 disp('>>> Saving...');
                    %                     stSaveInfo.mBFedData = mBFedData;
                    %                     stSaveInfo.mBEedData = mBEedData;
                    %                     switch stBFInfo.sDirection
                    %                         case 'lateral'
                    %                             stSaveInfo.mImgX = mImgX;
                    %                             stSaveInfo.mImgZ = mImgZ;
                    %                         case 'elevational'
                    %                             stSaveInfo.mImgY = mImgY;
                    %                             stSaveInfo.mImgZ = mImgZ;
                    %                             stSaveInfo.aYAxis_DSC = aYAxis;
                    %                             stSaveInfo.aZAxis_DSC = aZAxis;
                    %                             stSaveInfo.mImg_DSC = mImg;
                    %                             stSaveInfo.a_closest_angles = a_closest_angles;
                    %                     end
                    %                     stSaveInfo.stRFInfo = stRFInfo;
                    %                     stSaveInfo.stBFInfo = stBFInfo;
                    %                     stSaveInfo.stTRInfo = stTRInfo;
                    %                     stSaveInfo.stTxInfo = stTxInfo;
                    
                    if(strcmp(stBFInfo.sMode, 'SA'))
                        stSaveInfo.stSaInfo = stSAInfo;
                        stSaveInfo.mImg_Syn = mImg_Syn;
                    end
                    
                    if(s_idx == 1 && e_idx == 1 && m_idx == 1)
                        stParam.stRFInfo = stRFInfo;
                        stParam.stBFInfo = stBFInfo;
                        stParam.stTRInfo = stTRInfo;
                        stParam.stTxInfo = stTxInfo;
                        switch stBFInfo.sDirection
                            case 'lateral'
                                stParam.mImgX = mImgX;
                                stParam.mImgZ = mImgZ;
                            case 'elevational'
                                stParam.mImgY = mImgY;
                                stParam.mImgZ = mImgZ;
                                stParam.aYAxis_DSC = aYAxis;
                                stParam.aZAxis_DSC = aZAxis;
                        end
                        save([dir_tmp '/Parameters'], 'stParam', '-v7.3');
                    end
                    stSaveInfo.mImg_DSC = mImg;
                    stSaveInfo.a_closest_angles = a_closest_angles;
                    
                    dir_error = [dir_save '/error_' num2str(nError_range) '/Sample' num2str(s_idx,'%.3d')];
                    mkdir(dir_error);
                    
                    
                    switch stBFInfo.sMode
                        case 'Conv'
                            save([dir_error '/[Save_' stBFInfo.sDirection ']' stBFInfo.sMode '_Focus' num2str(stTxInfo.nTxFocus*1e3)],'stSaveInfo','-v7.3');
                            %                             save([sSaveDir '/[Save_' stBFInfo.sDirection ']' stBFInfo.sMode '_Natural_Focus'],'stSaveInfo','-v7.3');
                        case 'SA'
                            save([dir_error '/[Save_' stBFInfo.sDirection ']' stBFInfo.sMode '_VS' num2str(stSAInfo.nVSPos*1e3) '_Syn' num2str(stSAInfo.nNoSyn,'%.3d') ],'stSaveInfo', '-v7.3');
                            %                             save([sSaveDir '/[Save_' stBFInfo.sDirection ']' stBFInfo.sMode '_Natural_Syn' num2str(stSAInfo.nNoSyn,'%.3d') ],'stSaveInfo', '-v7.3');
                    end
                end
            end
        end
    end
end

disp('>>> all process completed');







