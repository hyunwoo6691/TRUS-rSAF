clc; clear; close all;

addpath('./functions');

no_samples = 1;
bSave = 1;
% error_case = 0;
error_case = [0 0.1 0.2 0.5 1 2 5];
%%
dir_ = uigetdir('./data','');

dir_data = [dir_ '/raw_data'];
scanline_data = dir(dir_data);
flag = 0;
if(strcmp(scanline_data(3).name,'.DS_Store')), flag = 1; end
scanline_data = scanline_data(3+flag:end);


dir_save = [dir_ '/errors_bf'];
mkdir(dir_save);
%%
spatial_filter = 'boxcar';
%% select element
close all;
element = 77;
% element = 45;

num_syn = 7;

mode = {'Conv','SA'};

d_theta = 0.4724;

%%
fov_original = 60;
scanning_original = 2048;
angle_original = linspace(-0.5*fov_original, 0.5*fov_original, scanning_original);

angle_transmit = fliplr(angle_original(342:1707)); % actually transmitted angle (fliplr : just image orientation))
%%
% acoustic parameter
stRFInfo.nFc = 6.5e6;
stRFInfo.nC = 1540;
stRFInfo.nFs = 4*stRFInfo.nFc;
stRFInfo.nUnitDis = stRFInfo.nC / stRFInfo.nFs;
stRFInfo.nLambda = stRFInfo.nC / stRFInfo.nFc;

% transducer parameter
stTRInfo.nPitch = 430e-6;
stTRInfo.nHeight = 5e-3;
stTRInfo.nEleFocus = 20e-3;
stTRInfo.nNoEle = 128;

% beamforming parameter
stBFInfo.nDth = 72e-3;
stBFInfo.nDthSpl = ceil(stBFInfo.nDth/stRFInfo.nUnitDis*2);
stBFInfo.nFnum = 1; % receive f-number
stBFInfo.nRadius = 10e-3;
stBFInfo.nCh = 1;
stBFInfo.nFOV = abs(angle_transmit(1)-angle_transmit(end));
stBFInfo.nDTheta = d_theta;
stBFInfo.nScline = round(stBFInfo.nFOV/stBFInfo.nDTheta);
stBFInfo.sDirection = 'elevational';
stBFInfo.sWindow = 'boxcar';

nDelayOff = 2*(sqrt((0.5*stTRInfo.nHeight)^2+stTRInfo.nEleFocus^2) - stTRInfo.nEleFocus)/stRFInfo.nC;


% aSclineTheta = linspace(-0.5*fov_original, 0.5*fov_original, 128);
% aSclineTheta = aSclineTheta(22:107); % ground truth angles (desired angle)


fov = stBFInfo.nFOV;
scline = stBFInfo.nScline;
% scline = numel(aSclineTheta);
radius = stBFInfo.nRadius;
depth = stBFInfo.nDth;
depth_spl = stBFInfo.nDthSpl;

% stBFInfo.nScline = scline;
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
    stParam.mImgY = mImgY;
    stParam.mImgZ = mImgZ;
    save([dir_ '/Parameters'], 'stParam', '-v7.3');
end
%%
for s_idx = 1:no_samples
    tic;
    disp(['>>> ' num2str(s_idx,'%.3d') '/' num2str(no_samples,'%.3d') '...']);
    for e_idx = 1:numel(error_case)
        disp(['    Error : ' num2str(error_case(e_idx)) ]);
        
        disp('      Data loading...');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% RF data extract for beamforming %%%%%%%
        nError_range = error_case(e_idx);
        aError = normrnd(0, nError_range, [1, scline]); % Gaussian distribution, mean: 0, std: nError_range
       
        a_transmitted_angles = aSclineTheta; % angles should be transmitted
        a_used_angles = a_transmitted_angles + aError; % actually transmitted angle
        
        a_closest_angles = zeros(1, scline);
        vRcvData = [];
        for t = 1:scline
            if(mod(t,round(scline/4))==0), disp(['      ' num2str(round(100*t/scline)) '%...']); end
            angle_tmp = a_used_angles(t);
            idx = find(abs(angle_transmit - angle_tmp) == min(abs(angle_transmit - angle_tmp)));
            a_closest_angles(t) = angle_transmit(idx(1));
            % load raw data
            load([dir_data '/' scanline_data(idx(1)).name '/raw_data.mat']); % get planewave transmitted at 0 deg
            %                 load([dir_data '/' data_list(idx).name '/Scanline_006.mat']); % get planewave transmitted at 0 deg
            vRcvData = cat(3, vRcvData, RF);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Beamforming
        for m_idx = 1:numel(mode)
            mode_tmp = mode{m_idx};
            disp(['      Mode : ' mode_tmp]);
            stBFInfo.sMode = mode_tmp;
            
            if(strcmp(mode_tmp, 'Conv'))
                % CON
                for ele_idx = element
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
                    
                    
                    for ele_idx = element
                        [mBFedData, vBFedData, vSynReg] = fElevational_SA_apod(vRcvData, stRFInfo, stBFInfo, stTRInfo, stSAInfo, aSclineTheta, mImgY, mImgZ, nDelayOff, ele_idx);
                        
                        dir_error = [dir_save '/error_' num2str(error_case(e_idx)) '/Sample' num2str(s_idx,'%.3d') '/Element_' num2str(ele_idx)];
                        if(~exist(dir_error, 'dir'))
                            mkdir(dir_error);
                        end
                        
                        if(bSave)
                            stSaveInfo.mBFedData = mBFedData;
                            stSaveInfo.vSynReg = vSynReg;
                            save([dir_error '/[Save]' stBFInfo.sMode '_VS' num2str(stTRInfo.nEleFocus*1e3) '_Syn' num2str(stSAInfo.nNoSyn,'%.3d') ],'stSaveInfo', '-v7.3');
                            clear stSaveInfo
                        end
                    end
                end
            end
        end
    end
    toc;
end