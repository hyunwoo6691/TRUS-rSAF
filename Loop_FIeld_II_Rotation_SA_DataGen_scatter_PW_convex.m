% Field II data generate simulation
% Lateral, elevational direction
% Elevation --> element rotating
% Data generation for 3D volumetric image

clc; clear; close all;

sCurrentPath = pwd;
addpath(sCurrentPath);
addpath('field_ii');
addpath('functions');

save_dir_master = [sCurrentPath '/data/'];

load_scatter = 0;
save_scatter = 0; scatter_name = 'prostate_mimic';

field_init(0);
%% Set parameter
stRFInfo.nFc = 6e6;
stRFInfo.nC = 1540;
stRFInfo.nFs = 4*stRFInfo.nFc;
stRFInfo.nAttenFactor = 0.5/(1e-2*1e6); % Attenuation factor[dB/(m*Hz)] % 0.5 dB/(cm*MHz)
stRFInfo.nUnitDis = stRFInfo.nC / stRFInfo.nFs;
stRFInfo.nLambda = stRFInfo.nC / stRFInfo.nFc;

stTRInfo.sType = 'convex'; % Transducer geometry : linear or convex
stTRInfo.nPitch = 210e-6;
stTRInfo.nWidth = 200e-6;
stTRInfo.nHeight = 7e-3;
stTRInfo.nKerf = stTRInfo.nPitch - stTRInfo.nWidth;
stTRInfo.nRadius = 10e-3;
stTRInfo.nEleFocus = 35e-3;
stTRInfo.nLatFocus = 25e-3;
stTRInfo.nNoEle = 128;
stTRInfo.nFOV = 149; % deg

d_theta = stTRInfo.nFOV/(stTRInfo.nNoEle-1);

sDirection = 'lateral'; % lateral or elavational

stTxInfo.type = 'conventional';
stTxInfo.cycle = 1;
stTxInfo.f_num = 3;
stTxInfo.num_tx_scline = stTRInfo.nNoEle;
stTxInfo.apodization = 'boxcar';

nTc = 1 / stRFInfo.nFc;
nUpsampleFactor = 5;
nSample = 2450; % 7.26 cm view depth

nSimulFreq = stRFInfo.nFs * nUpsampleFactor;
dT = 1 /nSimulFreq;

aFocalPoint = [0 0 0];

%% element position
ele_angle = linspace(-0.5*(stTRInfo.nNoEle-1), 0.5*(stTRInfo.nNoEle-1), stTRInfo.nNoEle)*d_theta;

aElePosX = stTRInfo.nRadius*sind(ele_angle);
aElePosZ = stTRInfo.nRadius*cosd(ele_angle);

%% focal points
focal_pointX = (stTRInfo.nRadius + stTRInfo.nLatFocus) * sind(ele_angle);
focal_pointZ = (stTRInfo.nRadius + stTRInfo.nLatFocus) * cosd(ele_angle);

%% Phantom geometry
view_depth = 50e-3;
% view_depth = 92e-3;
margin = 0;
% amp_high_reflect = 20;
amp_point_target = 20;

%% 
numX = 25;
numZ = 25;

totScat = numX * numZ;

metalSize = [20e-3, 5e-3]; % x, z

slope_ = metalSize(2)/metalSize(1); 
intersectZ = 25e-3;

m_point_posTmp = zeros(2*totScat,3);

%% Phantom - sphere object
if(~load_scatter)
    %     axially positioned
%     m_point_pos = ... % x, y, z
%         [0, 0, 10e-3;
%         0, 0, 20e-3;
%         0, 0, 30e-3;
%         0, 0, 40e-3;
%         0, 0, 50e-3;
%         0, 0, 60e-3;
%         0, 0, 70e-3;];
    
    m_point_posTmp(:,1) = (rand(2*totScat,1)-0.5)*metalSize(1); % x
    m_point_posTmp(:,3) = (rand(2*totScat,1)-0.5)*metalSize(2) + intersectZ;% z
    
    idx_ = m_point_posTmp(:,3) < (slope_ * m_point_posTmp(:,1) + intersectZ);
    
    m_point_pos = [m_point_posTmp(idx_,1) m_point_posTmp(idx_,2) m_point_posTmp(idx_,3)];
    
%     amp_point = amp_point_target*ones(size(m_point_pos,1),1);
    amp_point = amp_point_target*rand(size(m_point_pos,1),1);
    
    m_scat_total = m_point_pos;
    m_amp_total = amp_point;
    
    if(save_scatter)
        stScatInfo.scat_pos = m_scat_total;
        stScatInfo.scat_mag = m_amp_total;
        dir_scat = [sCurrentPath '/scatterer']; mkdir(dir_scat);
        save([dir_scat '/' scatter_name], 'stScatInfo');
    end
else % load the scatter file
    [scat_name, dir_scat] = uigetfile('','');
    load([dir_scat '/' scat_name]);
    m_scat_total = stScatInfo.scat_pos;
    m_amp_total = stScatInfo.scat_mag;
    clear stScatInfo
end
% Visualize
% figure(1);
% scatter3(m_scat_total(:,1), m_scat_total(:,2), m_scat_total(:,3),m_amp_total,'filled');
% xlabel('Lateral'); ylabel('Elevational'); zlabel('Axial');
% axis equal; axis tight; grid on;

figure(2);
scatter(m_scat_total(:,1)*1e3, m_scat_total(:,3)*1e3,m_amp_total,'filled');
xlabel('lateral'); ylabel('axial');
axis equal; %axis tight; 
grid on;
ylim([15 35]);
xlim([-30 30]);
set(gca,'Ydir','reverse');
%% Set field
set_field('c',stRFInfo.nC); % Set speed of sound
set_field('Freq_att',stRFInfo.nAttenFactor); % frequency dependency attenuation in [dB/(m*Hz)] around the center frequency
set_field('att', stRFInfo.nFc*stRFInfo.nAttenFactor); % Freqency independent attenuation[dB/m]
set_field('att_f0',stRFInfo.nFc); % Attenuation center frequency[Hz]
set_field('use_att',1); % attenuation on/off
set_field('fs',nSimulFreq); % set the sampling frequency
set_field('use_rectangles',1); % use ractangles for apertures

%% TX: Generate Sources
nFarFieldDepth_x = 0.1e-3; % depth over nFarFieldDepth_x is assumed to be far field on x-z plane
% For accurate directivity pattern, the number of source is more than 5
nFarFieldDepth_y = 0.1e-3;     % depth over nFarFieldDepth_x is assumed to be far field on y-z plane

nMathEleSize_x = sqrt(nFarFieldDepth_x * 4 * stRFInfo.nLambda); % H < sqrt(4*lambda*z)
nMathEleSize_y = sqrt(nFarFieldDepth_y * 4 * stRFInfo.nLambda);

nSubDivNum_x = ceil(stTRInfo.nWidth / nMathEleSize_x);
nSubDivNum_y = ceil(stTRInfo.nHeight/ nMathEleSize_y);

pTxAperture = xdc_convex_focused_array(stTRInfo.nNoEle, stTRInfo.nWidth, stTRInfo.nHeight, stTRInfo.nKerf, stTRInfo.nRadius, stTRInfo.nEleFocus, nSubDivNum_x, nSubDivNum_y, aFocalPoint);
pRxAperture = xdc_convex_focused_array(stTRInfo.nNoEle, stTRInfo.nWidth, stTRInfo.nHeight, stTRInfo.nKerf, stTRInfo.nRadius, stTRInfo.nEleFocus, nSubDivNum_x, nSubDivNum_y, aFocalPoint);

xdc_baffle(pTxAperture,0);
xdc_baffle(pRxAperture,-1);

% Show aperture
% figure(1);
% subplot(2,1,1); show_xdc(pTxAperture); title('Tx Aperture');
% subplot(2,1,2); show_xdc(pRxAperture); title('Rx Aperture');

%% Impulse Response
at2 = 0:dT:2.3*nTc;
aTransImpulsResp = sin(2*pi*stRFInfo.nFc*at2);
aTransImpulsResp = aTransImpulsResp.*(hanning(numel(aTransImpulsResp))'); % Transducer's impulse response
xdc_impulse(pTxAperture,aTransImpulsResp); % Tx impulse response

%% Excitation Pulse
at1 = 0:dT:stTxInfo.cycle*nTc;
aPulseSeq  = sin(2*pi*stRFInfo.nFc*at1);  % excitation pulse
xdc_excitation(pTxAperture,aPulseSeq);

figure(2)
subplot(2,1,1); plot(at1,aPulseSeq), title('Exitation Pulse')
subplot(2,1,2); plot(at2,aTransImpulsResp), title('Transducer`s Impulse response')

nLag = round(((numel(aPulseSeq)+numel(aTransImpulsResp)-1)+numel(aTransImpulsResp)-1)/2);

%% Save parameter
stSaveInfo.stRFInfo = stRFInfo;
stSaveInfo.stTRInfo = stTRInfo;
stSaveInfo.stTxInfo = stTxInfo;

switch stTxInfo.type
    case 'conventional'
        folder_name = ['[Focused]' 'EleFocus' num2str(stTRInfo.nEleFocus*1e3) ...
            'mm_' '_Tx' num2str(stTxInfo.num_tx_scline)];
    case 'plane_wave'
        folder_name = ['[Plane_wave]Volume_Steer_' num2str(stTxInfo.max_steering) 'deg_' num2str(stTxInfo.num_tx_plane) 'compounding' ...
            '_FOV' num2str(FOV_theta) '_Tx' num2str(stTxInfo.num_tx_vol)];
end

save_dir_parameter = [save_dir_master folder_name];
save_dir_data = [save_dir_parameter '/Raw data'];
mkdir(save_dir_parameter);
mkdir(save_dir_data);

save([save_dir_parameter '/stSaveInfo'], 'stSaveInfo');

%% Loop for pressure calculation at each position
mScatXYZPos_Rot = m_scat_total;

vRcvData = zeros(nSample, stTRInfo.nNoEle, stTxInfo.num_tx_scline);

for t_idx = 1:stTxInfo.num_tx_scline
    disp(['>>> scanline ' num2str(t_idx,'%.3d')]); tic;
    
    apod_window = fieldII_get_apod(t_idx, stTRInfo.nLatFocus, stTxInfo.f_num, stTRInfo.nNoEle, stTRInfo.nPitch);
    delay_curve = fieldII_get_delay_convex(apod_window, focal_pointX(t_idx), focal_pointZ(t_idx), aElePosX, aElePosZ, stRFInfo.nC);
    
    xdc_apodization(pTxAperture, 0, apod_window);
    xdc_apodization(pRxAperture, 0, ones(1, stTRInfo.nNoEle));
    xdc_focus_times(pTxAperture, 0, delay_curve);
    xdc_center_focus(pRxAperture, [0 0 0]);% Set the origin for the dynamic focusing line.
    xdc_focus_times(pRxAperture, 0, zeros(1,stTRInfo.nNoEle));
    
    %%% Calculate the received signals from a collection of scatterers for all elements
    %         [mRF, nStartTime] = calc_scat_multi(pTxAperture, pRxAperture, m_point_pos, amp_point);
    [mRF, nStartTime] = calc_scat_multi(pTxAperture, pRxAperture, mScatXYZPos_Rot, m_amp_total);
    
    %%% Zero-padding (so that the first sample is captured at t = 0)
    nPadLen = round(nStartTime/dT);
    mRF_pad = padarray(mRF,[nPadLen 0],'pre');
    
    %%% Truncation or Zero-Padding (so that # of samples = RFInfo.nSample)
    nUpSample = nSample*nUpsampleFactor;
    nFirstSamIdx = 1 + nLag;
    nLastSamIdx = nFirstSamIdx + nUpSample -1;
    
    if(size(mRF_pad,1) > nLastSamIdx)
        mRF_trc = mRF_pad(nFirstSamIdx:nLastSamIdx,:);
    elseif(size(mRF_pad,1) < nLastSamIdx)
        mRF_pad2 = padarray(mRF_pad,[nLastSamIdx-size(mRF_pad,1) 0],'post');
        mRF_trc = mRF_pad2(nFirstSamIdx:nLastSamIdx,:);
    end
    
    %%% Decimate to lower sampling frequency (nFs)
    mRF_deci = mRF_trc(1:nUpsampleFactor:end,:);
    
    %%% Save RF data
    mRcvData = mRF_deci;
    
    vRcvData(:,:,t_idx) = mRcvData;
    
    % save scanline data
    save([save_dir_data '/scanline_' num2str(t_idx, '%.3d')], 'mRcvData');
    toc;
end

%%
field_end;