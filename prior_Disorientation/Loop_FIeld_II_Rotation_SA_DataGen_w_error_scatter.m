% Field II data generate simulation
% Lateral, elevational direction
% Elevation --> element rotating
% Data generation for 3D volumetric image

clc; clear; close all;

sCurrentPath = pwd;
addpath(sCurrentPath);
addpath('Field_II_ver_3_24_mac');
addpath('02.Functions');

save_dir_master = [sCurrentPath '/01.Data/'];

field_init(0);
%% Set parameter
stRFInfo.nFc = 6.5e6;
stRFInfo.nC = 1540;
stRFInfo.nFs = 4*stRFInfo.nFc;
stRFInfo.nAttenFactor = 0.5/(1e-2*1e6); % Attenuation factor[dB/(m*Hz)] % 0.5 dB/(cm*MHz)
stRFInfo.nUnitDis = stRFInfo.nC / stRFInfo.nFs;
stRFInfo.nLambda = stRFInfo.nC / stRFInfo.nFc;

stTRInfo.sType = 'linear'; % Transducer geometry : linear or convex
stTRInfo.nPitch = 430e-6;
stTRInfo.nWidth = 420e-6;
stTRInfo.nHeight = 5e-3;
stTRInfo.nKerf = stTRInfo.nPitch - stTRInfo.nWidth;
stTRInfo.nEleFocus = 25e-3;
stTRInfo.nNoEle = 128;

sDirection = 'elevational'; % lateral or elavational

stTxInfo.nCycle = 1;
stTxInfo.nNoTx_volume = 2048; % volume scan
stTxInfo.nNoTx_plane = stTRInfo.nNoEle; % number of tx at each angle
stTxInfo.nTxFocus = 40e-3;
stTxInfo.nTxFnum = 3;
stTxInfo.sApod = 'boxcar';

aElePosX = stTRInfo.nPitch * linspace(-0.5*(stTRInfo.nNoEle - 1), 0.5*(stTRInfo.nNoEle - 1), stTRInfo.nNoEle);
aElePosZ = zeros(1, stTRInfo.nNoEle);

nTc = 1 / stRFInfo.nFc;
nUpsampleFactor = 5;
nSample = 2450; % 7.26 cm view depth

nSimulFreq = stRFInfo.nFs * nUpsampleFactor;
dT = 1 /nSimulFreq;

aFocalPoint = [0 0 stTxInfo.nTxFocus];
%% Rotate angle
%%%%%%%%%%%%
radius = 5e-3;
FOV_theta = 70;

aPosRng = linspace(0.5* FOV_theta, -0.5*FOV_theta, stTxInfo.nNoTx_volume);
%%%%%%%%%%%%

%% Phantom geometry
view_depth = 70e-3;
margin = 0;
amp_high_reflect = 10;
amp_point_target = 100;

% dimension of the phantom ( square phantom )
size_x = abs(aElePosX(end) - aElePosX(1)) + margin;
size_y = 2*(view_depth + radius)*sind(0.5*FOV_theta) + margin;
size_z = view_depth + margin;
start_z = 0;

volume = size_x * size_y * size_z;

% Scatterers
% density of scatterer (density = num_scatters/volume; [#Sc/mm^3]);
density = 5;
num_scatters = round(density * volume*1e9);

x_pos = (rand(num_scatters, 1) - 0.5) * size_x;
y_pos = (rand(num_scatters, 1) - 0.5) * size_y;
z_pos = rand(num_scatters, 1) *  size_z + start_z;

amp_scatt = rand(num_scatters, 1);
m_scat_pos = cat(2, x_pos, y_pos, z_pos);

% Point targets
% x, y, z
num_axial = 7;
num_ele = 9;
num_lat = 5;
space = 10e-3;

axial_pos = linspace(10e-3, num_axial*space, num_axial);
ele_pos = linspace(-0.5*space*(num_ele-1), 0.5*space*(num_ele-1), num_ele);
lat_pos = linspace(-0.5*space*(num_lat-1), 0.5*space*(num_lat-1), num_lat);
m_point_pos = [];
for a_idx = 1:num_axial
   for e_idx = 1:num_ele
       for l_idx = 1:num_lat
           pos_tmp = [lat_pos(l_idx), ele_pos(e_idx), axial_pos(a_idx)];
           m_point_pos = cat(1, m_point_pos, pos_tmp);
       end
   end
end

amp_point = amp_point_target*ones(size(m_point_pos,1),1);

% Cysts


% Visualize
figure(1);
scatter3(m_point_pos(:,1), m_point_pos(:,2), m_point_pos(:,3),amp_point,'filled');
xlabel('Lateral'); ylabel('Elevational'); zlabel('Axial');
axis equal; axis tight; grid on;

m_scat_total = cat(1, m_scat_pos, m_point_pos);
m_amp_total = cat(1, amp_scatt, amp_point);

%% Set field
set_field('c',stRFInfo.nC); % Set speed of sound
set_field('Freq_att',stRFInfo.nAttenFactor); % frequency dependency attenuation in [dB/(m*Hz)] around the center frequency
set_field('att', stRFInfo.nFc*stRFInfo.nAttenFactor); % Freqency independent attenuation[dB/m]
set_field('att_f0',stRFInfo.nFc); % Attenuation center frequency[Hz]
set_field('use_att',1); % attenuation on/off
set_field('fs',nSimulFreq); % set the sampling frequency
set_field('use_rectangles',1); % use ractangles for apertures

%% TX: Generate Sources
nFarFieldDepth_x = 0.01e-3; % depth over nFarFieldDepth_x is assumed to be far field on x-z plane
% For accurate directivity pattern, the number of source is more than 5
nFarFieldDepth_y = 0.01e-3;     % depth over nFarFieldDepth_x is assumed to be far field on y-z plane

nMathEleSize_x = sqrt(nFarFieldDepth_x * 4 * stRFInfo.nLambda); % H < sqrt(4*lambda*z)
nMathEleSize_y = sqrt(nFarFieldDepth_y * 4 * stRFInfo.nLambda);

nSubDivNum_x = ceil(stTRInfo.nWidth / nMathEleSize_x);
nSubDivNum_y = ceil(stTRInfo.nHeight/ nMathEleSize_y);

pTxAperture = xdc_focused_array (stTRInfo.nNoEle, stTRInfo.nWidth, stTRInfo.nHeight, stTRInfo.nKerf, stTRInfo.nEleFocus, nSubDivNum_x, nSubDivNum_y, aFocalPoint);
pRxAperture = xdc_focused_array (stTRInfo.nNoEle, stTRInfo.nWidth, stTRInfo.nHeight, stTRInfo.nKerf, stTRInfo.nEleFocus, nSubDivNum_x, nSubDivNum_y, aFocalPoint);

xdc_baffle(pTxAperture,0);
xdc_baffle(pRxAperture,-1);
%% Show aperture
% figure(1);
% subplot(2,1,1); show_xdc(pTxAperture); title('Tx Aperture');
% subplot(2,1,2); show_xdc(pRxAperture); title('Rx Aperture');

%% Impulse Response
at2 = 0:dT:2.3*nTc;
aTransImpulsResp = sin(2*pi*stRFInfo.nFc*at2);
aTransImpulsResp = aTransImpulsResp.*(hanning(numel(aTransImpulsResp))'); % Transducer's impulse response
xdc_impulse(pTxAperture,aTransImpulsResp); % Tx impulse response

%% Excitation Pulse
at1 = 0:dT:stTxInfo.nCycle*nTc;
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
stSaveInfo.nRadius = radius;
stSaveInfo.FOV_theta = FOV_theta;
stSaveInfo.aPosRng = aPosRng;

folder_name = ['Volume_LatFocus' num2str(stTxInfo.nTxFocus*1e3) 'mm_EleFocus' num2str(stTRInfo.nEleFocus*1e3) 'mm_FOV' num2str(FOV_theta) '_Tx' num2str(stTxInfo.nNoTx_volume)];

save_dir_parameter = [save_dir_master folder_name];
save_dir_data = [save_dir_parameter '/Raw data'];
mkdir(save_dir_parameter);
mkdir(save_dir_data);

save([save_dir_parameter '/stSaveInfo'], 'stSaveInfo');

%% Search latest progress
raw_list = dir(save_dir_data);
if(numel(raw_list)>=3) % if already processed planes 
    flag = 0;
    if(strcmp(raw_list(3).name, '.DS_Store')), flag = 1; end
    raw_list = raw_list(3+flag:end);
    start_p_idx = numel(raw_list);
else % if it is first
    start_p_idx = 1;
end

last_done_plane = raw_list(start_p_idx).name;
dir_ = [save_dir_data '/' last_done_plane];
scanline_list = dir(dir_);
if(numel(scanline_list)>=3) % if already processed scanline data
    flag = 0;
    if(strcmp(scanline_list(3).name, '.DS_Store')), flag = 1; end
    scanline_list = scanline_list(3+flag:end);
    start_t_idx = numel(scanline_list)+1;
else % if it is first
    start_t_idx = 1;
end

%% Loop for pressure calculation at each position
mScatXYZPos_Rot = m_point_pos;

for p_idx = start_p_idx:stTxInfo.nNoTx_volume
    disp(['>>> Volume scan [' num2str(p_idx) '/' num2str(stTxInfo.nNoTx_volume) ']...']);
    % Field rotation
    nTheta = aPosRng(p_idx);
    aTheta_Zero = atand(m_point_pos(:,2)./(radius+m_point_pos(:,3)));
    aTheta_Prime = aTheta_Zero + nTheta;
    
    aPosZ = ((m_point_pos(:,3)+radius)./cosd(aTheta_Zero)) .* cosd(aTheta_Prime) - radius;
    aPosY = ((m_point_pos(:,3)+radius)./cosd(aTheta_Zero)) .* sind(aTheta_Prime);
    
    mScatXYZPos_Rot(:,3) = aPosZ;
    mScatXYZPos_Rot(:,2) = aPosY;
    % x-axis doesn't change
    
    vRcvData = zeros(nSample, stTRInfo.nNoEle, stTxInfo.nNoTx_plane);
    dir_tmp = [save_dir_data '/RAW_' num2str(p_idx, '%.4d')];
    mkdir(dir_tmp);
    for t_idx = start_t_idx:stTxInfo.nNoTx_plane
        disp(['    Scanline : ' num2str(t_idx)]);tic;
        apod_window = fieldII_get_apod(t_idx,stTxInfo.nTxFocus, stTxInfo.nTxFnum, stTRInfo.nNoEle, stTRInfo.nPitch);
        delay_curve = fieldII_get_delay(t_idx, apod_window, stTxInfo.nTxFocus, aElePosX, stRFInfo.nC);
        
        xdc_apodization(pTxAperture, 0, apod_window);
        xdc_apodization(pRxAperture, 0, ones(1, stTRInfo.nNoEle));
        xdc_focus_times(pTxAperture, 0, delay_curve);
        xdc_center_focus(pRxAperture, [0 0 0]);% Set the origin for the dynamic focusing line.
        xdc_focus_times(pRxAperture, 0, zeros(1,stTRInfo.nNoEle));
        
        %%% Calculate the received signals from a collection of scatterers for all elements
%         [mRF, nStartTime] = calc_scat_multi(pTxAperture, pRxAperture, m_point_pos, amp_point);
        [mRF, nStartTime] = calc_scat_multi(pTxAperture, pRxAperture, mScatXYZPos_Rot, amp_point);
        
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
        save([dir_tmp '/Scanline_' num2str(t_idx, '%.3d')], 'mRcvData');
        toc;
    end
end

%%
field_end;