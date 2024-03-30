clc; clear; close all;

addpath('functions');

bDisplay = 0;
nDelayOff = 0;

bSave = 1;
bAWGN = 1; dB_ = 20;

num_syn = [3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65];

error_case = 0;
mode = {'Conv', 'SA'};
% mode = {'SA'};

%%
dir_master = uigetdir('./data',''); %% master folder. e.g., rsaf_optim_fele
folder_list = dir(dir_master);
flag = 0;
if(strcmp(folder_list(3).name, '.DS_Store')), flag = 1; end
folder_list = folder_list(3+flag:end);

for folder_idx = 1:numel(folder_list)
    
    folder_tmp = folder_list(folder_idx).name;
    
    disp(['>>> ' folder_tmp]);
    dir_ = [dir_master '/' folder_tmp];
    
    for noise_ = 1:2
        % no_samples = 300;
        % error_case = [0.1 0.2 0.5 1 2 5];
        bAWGN = noise_;
        if(bAWGN)
            no_samples = 5;
        else
            no_samples = 1;
        end
        %% Setup beamforming parameter
        
        % load info file
        load([dir_ '/stSaveInfo.mat']);
        stRFInfo = stSaveInfo.stRFInfo;
        stTRInfo = stSaveInfo.stTRInfo;
        stTxInfo = stSaveInfo.stTxInfo;
        nRadius = stSaveInfo.nRadius;
        nFOV_Theta_fieldII = stSaveInfo.FOV_theta;
        aPosRng_fieldII = fliplr(stSaveInfo.aPosRng);
        clear stSaveInfo;
        
        % load data list
        dir_data = [dir_ '/Raw data'];
        data_list = dir(dir_data);
        flag = 0;
        if(strcmp(data_list(3).name, '.DS_Store')), flag = 1; end
        data_list = data_list(3+flag:end);
        
        if(bAWGN)
            %     dir_save = [dir_ '/errors_bf_AWGN'];
            dir_save = [dir_ '/errors_bf_AWGN_noMask'];
        else
            %     dir_save = [dir_ '/errors_bf'];
            dir_save = [dir_ '/errors_bf_noMask'];
        end
        mkdir(dir_save);
        
        stBFInfo.nDth = 72e-3;
        stBFInfo.nDthSpl = ceil(stBFInfo.nDth/stRFInfo.nUnitDis*2);
        stBFInfo.nFnum = 1; % receive f-number
        stBFInfo.nRadius = nRadius;
        stBFInfo.nCh = 1;
        stBFInfo.nFOV = 60;
        stBFInfo.nScline = 128;
        stBFInfo.sDirection = 'elevational';
        stBFInfo.sWindow = 'boxcar';
        
        fov = stBFInfo.nFOV;
        scline = stBFInfo.nScline;
        radius = stBFInfo.nRadius;
        depth = stBFInfo.nDth;
        depth_spl = stBFInfo.nDthSpl;
        
        aSclineTheta = linspace(-0.5*fov, 0.5*fov, scline); % Ground truth transmitted angle
        nDelta_theta = abs(aSclineTheta(1)-aSclineTheta(2));
        % Beamforming grid
        aDth = linspace(radius, radius+depth, depth_spl);
        
        mDth = repmat(aDth', 1, scline);
        mTheta = repmat(aSclineTheta, numel(aDth), 1);
        
        mImgZ = mDth .* cosd(mTheta);
        mImgY = mDth .* sind(mTheta);
        
        if(bSave)
            stParam.stRFInfo = stRFInfo;
            stParam.stBFInfo = stBFInfo;
            stParam.stTRInfo = stTRInfo;
            stParam.stTxInfo = stTxInfo;
            stParam.mImgY = mImgY;
            stParam.mImgZ = mImgZ;
            save([dir_ '/Parameters'], 'stParam', '-v7.3');
        end
        %%
        for s_idx = 1:no_samples
            disp(['>>> ' num2str(s_idx,'%.3d') '/' num2str(no_samples, '%.3d') '...']);
            for e_idx = 1:numel(error_case)
                disp(['    Error : ' num2str(error_case(e_idx)) ]);
                
                disp('      Data loading...');
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%% RF data extract for beamforming %%%%%%%
                nError_range = error_case(e_idx);
                aError = normrnd(0, nError_range, [1, scline]); % Gaussian distribution, mean: 0, std: nError_range
                
                a_transmitted_angles = linspace(-0.5*fov, 0.5*fov, scline); % angles should be transmitted
                a_used_angles = a_transmitted_angles + aError; % actually transmitted angle
                
                a_closest_angles = zeros(1, scline);
                vRcvData = [];
                for t = 1:scline
                    angle_tmp = a_used_angles(t);
                    idx = find(abs(aPosRng_fieldII - angle_tmp) == min(abs(aPosRng_fieldII - angle_tmp)));
                    a_closest_angles(t) = aPosRng_fieldII(idx);
                    % load raw data
                    load([dir_data '/' data_list(idx).name '/Scanline_001.mat']); % get planewave transmitted at 0 deg
                    %                 load([dir_data '/' data_list(idx).name '/Scanline_006.mat']); % get planewave transmitted at 0 deg
                    vRcvData = cat(3, vRcvData, mRcvData);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if(bAWGN)
                    vRcvData = awgn(vRcvData, dB_, 'measured');
                end
                
                for m_idx = 1:numel(mode)
                    disp(['      Mode : ' mode{m_idx}]); tic;
                    stBFInfo.sMode = mode{m_idx};
                    
                    if(strcmp(mode{m_idx}, 'Conv'))
                        % if it is CON
                        for ele_idx = 64
                            [mBFedData] = fElevational_Conv(vRcvData, stRFInfo, stBFInfo, mImgY, mImgZ, nDelayOff, ele_idx);
                            
                            dir_error = [dir_save '/error_' num2str(nError_range) '/Sample' num2str(s_idx,'%.3d') '/Element_' num2str(ele_idx)];
                            if(~exist(dir_error, 'dir'))
                                mkdir(dir_error);
                            end
                            
                            stSaveInfo.mBFedData = mBFedData;
                            save([dir_error '/[Save]' stBFInfo.sMode '_Focus' num2str(stTRInfo.nEleFocus*1e3)],'stSaveInfo','-v7.3');
                            
                        end
                    else
                        % if it is rSAF
                        for syn_idx = 1:numel(num_syn)
                            disp(['               syn' num2str(num_syn(syn_idx))]);
                            %%%%%%%%%%%%%%%
                            stSAInfo.nNoSyn = num_syn(syn_idx);
                            %%%%%%%%%%%%%%%
                            if(stTRInfo.nEleFocus == 0)
                                nNaturalFocalDth = stTRInfo.nHeight^2/(4*stRFInfo.nLambda);
                                stSAInfo.nVSPos = nNaturalFocalDth;
                            else
                                stSAInfo.nVSPos = stTRInfo.nEleFocus;
                            end
                            
                            for ele_idx = 64
                                [mBFedData, vBFedData, vSynReg] = fElevational_SA(vRcvData, stRFInfo, stBFInfo, stTRInfo, stSAInfo, aSclineTheta, mImgY, mImgZ, nDelayOff, ele_idx);
                                
                                dir_error = [dir_save '/error_' num2str(nError_range) '/Sample' num2str(s_idx,'%.3d') '/Element_' num2str(ele_idx)];
                                if(~exist(dir_error, 'dir'))
                                    mkdir(dir_error);
                                end
                                
                                stSaveInfo.mBFedData = mBFedData;
                                stSaveInfo.vSynReg = vSynReg;
                                save([dir_error '/[Save]' stBFInfo.sMode '_VS' num2str(stTRInfo.nEleFocus*1e3) '_Syn' num2str(stSAInfo.nNoSyn,'%.3d') ],'stSaveInfo', '-v7.3');
                            end
                        end
                    end
                end
            end
            toc;
        end
    end
end
disp('>>> all process completed');


