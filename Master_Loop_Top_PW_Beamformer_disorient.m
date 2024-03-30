clc; clear; close all;

addpath('functions');

bSave = 1;
bAWGN = 0; dB_ = 20;

error_case = [0 0.1 0.2 0.5 1 2 5];
%%
dir_master = uigetdir('./data', 'Select folder');
folder_list = dir(dir_master);
mask_ = zeros(1,numel(folder_list));
for f_idx = 1:numel(folder_list)
    mask_(f_idx) = folder_list(f_idx).isdir;
end
folder_list = folder_list(logical(mask_));
folder_list = folder_list(3:end);

%% load disorientation data
load([dir_master '/disorientations.mat']);
assert(numel(error_case) == size(disorientations_,1),'Error cases do not match');

%%
for f_idx = 1:numel(folder_list)
    folder_name = folder_list(f_idx).name;
    switch folder_list(f_idx).name
        case 'REF'
            mode = {'Conv', 'SA'};
%             d_theta = 0.4724;
            d_theta = 0.2362;
            num_syn = 7;
        case 'rSAF'
            mode = {'SA'};
            d_theta = 0.2362;
            num_syn = 131;
    end
    
    %% Setup beamforming parameter
    dir_ = [dir_master '/' folder_name];
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
        dir_save = [dir_ '/errors_bf_AWGN_' num2str(dB_)];
        %     dir_save = [dir_ '/errors_bf_AWGN_noTarget'];
    else
        dir_save = [dir_ '/errors_bf'];
    end
    mkdir(dir_save);
    
    stBFInfo.nDth = 72e-3;
    stBFInfo.nDthSpl = ceil(stBFInfo.nDth/stRFInfo.nUnitDis*2);
    stBFInfo.nFnum = 1; % receive f-number
    stBFInfo.nRadius = nRadius;
    stBFInfo.nCh = 1;
    stBFInfo.nFOV = nFOV_Theta_fieldII;
    stBFInfo.nDTheta = d_theta;
    stBFInfo.nScline = round(stBFInfo.nFOV/stBFInfo.nDTheta);
    stBFInfo.sDirection = 'elevational';
    stBFInfo.sWindow = 'boxcar';
    
    nDelayOff = 2*(sqrt((0.5*stTRInfo.nHeight)^2+stTRInfo.nEleFocus^2) - stTRInfo.nEleFocus)/stRFInfo.nC;
    
    fov = stBFInfo.nFOV;
    scline = stBFInfo.nScline;
    radius = stBFInfo.nRadius;
    depth = stBFInfo.nDth;
    depth_spl = stBFInfo.nDthSpl;
    
    aSclineTheta = linspace(-0.5*(scline-1), 0.5*(scline-1), scline) * stBFInfo.nDTheta; % Ground truth transmitted angle
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
    for s_idx = (1:size(disorientations_,3))
        disp(['>>> [' folder_name '] ' num2str(s_idx,'%.3d') '/' num2str(size(disorientations_,3), '%.3d') '...']);
        for e_idx = 1:numel(error_case)
            disp(['    Error : ' num2str(error_case(e_idx)) ]); tic;
            
            disp('      Data loading...');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%% RF data extract for beamforming %%%%%%%
            aError = disorientations_(e_idx,:);
            
            a_transmitted_angles = aSclineTheta; % angles should be transmitted
            a_used_angles = a_transmitted_angles + aError; % actually transmitted angle
            
            a_closest_angles = zeros(1, scline);
            vRcvData = [];
            for t = 1:scline
                if(mod(t,round(scline/4))==0), disp(['      ' num2str(round(100*t/scline)) '%...']); end
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
                vRcvData_tmp = awgn(vRcvData, dB_, 'measured');
                %             vRcvData = vRcvData_tmp - vRcvData; % only noise data
                vRcvData = vRcvData_tmp; % noise + signal
            end
            
            
            % Beamforming
            for m_idx = 1:numel(mode)
                mode_tmp = mode{m_idx};
                disp(['      Mode : ' mode_tmp]);
                stBFInfo.sMode = mode_tmp;
                
                if(strcmp(mode_tmp, 'Conv'))
                    % CON
                    for ele_idx = 64
                        [mBFedData] = fElevational_Conv(vRcvData, stRFInfo, stBFInfo, mImgY, mImgZ, nDelayOff, ele_idx);
                        
                        dir_error = [dir_save '/error_' num2str(error_case(e_idx)) '/Sample' num2str(s_idx,'%.3d') '/Element_' num2str(ele_idx)];
                        if(~exist(dir_error, 'dir'))
                            mkdir(dir_error);
                        end
                        
                        if(bSave)
                            stSaveInfo.mBFedData = mBFedData;
                            save([dir_error '/[Save]' stBFInfo.sMode '_Focus' num2str(stTRInfo.nEleFocus*1e3)],'stSaveInfo','-v7.3');
                            clear stSaveInfo
                        end
                    end
                else
                    % rSAF
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
                            [mBFedData, vBFedData, vSynReg] = fElevational_SA_apod(vRcvData, stRFInfo, stBFInfo, stTRInfo, stSAInfo, aSclineTheta, mImgY, mImgZ, nDelayOff, ele_idx);
                            
                            dir_error = [dir_save '/error_' num2str(error_case(e_idx)) '/Sample' num2str(s_idx,'%.3d') '/Element_' num2str(ele_idx)];
                            if(~exist(dir_error, 'dir'))
                                mkdir(dir_error);
                            end
                            
                            if(bSave)
                                stSaveInfo.mBFedData = mBFedData;
                                stSaveInfo.vSynReg = vSynReg;
                                %%%%%%%%%%%%
                                stSaveInfo.vBFedData = vBFedData;
                                %%%%%%%%%%%%
                                save([dir_error '/[Save]' stBFInfo.sMode '_VS' num2str(stTRInfo.nEleFocus*1e3) '_Syn' num2str(stSAInfo.nNoSyn,'%.3d') ],'stSaveInfo', '-v7.3');
                                clear stSaveInfo
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


