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
stRFInfo.nFc = 6.5e6;
stRFInfo.nC = 1540;
stRFInfo.nFs = 4*stRFInfo.nFc;
stRFInfo.nAttenFactor = 0.5/(1e-2*1e6); % Attenuation factor[dB/(m*Hz)] % 0.5 dB/(cm*MHz)
stRFInfo.nUnitDis = stRFInfo.nC / stRFInfo.nFs;
stRFInfo.nLambda = stRFInfo.nC / stRFInfo.nFc;

stTRInfo.sType = 'linear'; % Transducer geometry : linear or convex
stTRInfo.nPitch = 430e-6;
stTRInfo.nWidth = 420e-6;
stTRInfo.nHeight = 7e-3;
stTRInfo.nKerf = stTRInfo.nPitch - stTRInfo.nWidth;

stTRInfo.nEleFocus = 5e-3;

stTRInfo.nNoEle = 128;

sDirection = 'elevational'; % lateral or elavational

stTxInfo.type = 'plane_wave';
stTxInfo.cycle = 1;
% stTxInfo.num_tx_vol = 2048; % volume scan
% stTxInfo.num_tx_vol = 128; % volume scan

% stTxInfo.num_tx_vol = 280; % volume scan - 133deg FOV
% stTxInfo.num_tx_vol = 280*1.25; % volume scan - 1/1.25 theta
% stTxInfo.num_tx_vol = 280*1.5; % volume scan - 1/1.5 theta
% stTxInfo.num_tx_vol = 280*1.75; % volume scan - 1/1.75 theta
% stTxInfo.num_tx_vol = 280*2; % volume scan - half theta

% stTxInfo.num_tx_vol = round(280*0.4724/0.1564); % volume scan - 10mm, -20dB


stTxInfo.num_tx_vol = round(280*0.4724/0.0787); % volume scan - 10mm, -20dB

stTxInfo.num_tx_plane = 1; % number of tx at each angle
stTxInfo.max_steering = 0; % [deg]
stTxInfo.apodization = 'boxcar';
stTxInfo.tx_angles = linspace(-stTxInfo.max_steering, stTxInfo.max_steering, stTxInfo.num_tx_plane);

aElePosX = stTRInfo.nPitch * linspace(-0.5*(stTRInfo.nNoEle - 1), 0.5*(stTRInfo.nNoEle - 1), stTRInfo.nNoEle);
aElePosZ = zeros(1, stTRInfo.nNoEle);

nTc = 1 / stRFInfo.nFc;
nUpsampleFactor = 5;
nSample = 2450; % 7.26 cm view depth

nSimulFreq = stRFInfo.nFs * nUpsampleFactor;
dT = 1 /nSimulFreq;

aFocalPoint = [0 0 0];
%% Rotate angle
%%%%%%%%%%%%
radius = 15e-3;
% radius = 10e-3;
% FOV_theta = 60;

% aPosRng = linspace(0.5* FOV_theta, -0.5*FOV_theta, stTxInfo.num_tx_vol);
% aPosRng = linspace(0.5*(stTxInfo.num_tx_vol-1), -0.5*(stTxInfo.num_tx_vol-1), stTxInfo.num_tx_vol) * 0.4724;
% aPosRng = linspace(0.5*(stTxInfo.num_tx_vol-1), -0.5*(stTxInfo.num_tx_vol-1), stTxInfo.num_tx_vol) * 0.4724/1.25;
% aPosRng = linspace(0.5*(stTxInfo.num_tx_vol-1), -0.5*(stTxInfo.num_tx_vol-1), stTxInfo.num_tx_vol) * 0.4724/1.5;
% aPosRng = linspace(0.5*(stTxInfo.num_tx_vol-1), -0.5*(stTxInfo.num_tx_vol-1), stTxInfo.num_tx_vol) * 0.4724/1.75;
% aPosRng = linspace(0.5*(stTxInfo.num_tx_vol-1), -0.5*(stTxInfo.num_tx_vol-1), stTxInfo.num_tx_vol) * 0.4724/2;

% aPosRng = linspace(0.5*(stTxInfo.num_tx_vol-1), -0.5*(stTxInfo.num_tx_vol-1), stTxInfo.num_tx_vol) * 0.1564;

% aPosRng = linspace(0.5*(stTxInfo.num_tx_vol-1), -0.5*(stTxInfo.num_tx_vol-1), stTxInfo.num_tx_vol) * 0.4724/4;
aPosRng = linspace(0.5*(stTxInfo.num_tx_vol-1), -0.5*(stTxInfo.num_tx_vol-1), stTxInfo.num_tx_vol) * 0.0787;

FOV_theta = abs(aPosRng(1)-aPosRng(end));
% aPosRng = 0;
%%%%%%%%%%%%

%% Phantom geometry
view_depth = 72e-3;
% view_depth = 92e-3;
margin = 0;
% amp_high_reflect = 20;
amp_point_target = 20;

%% Phantom - sphere object
if(~load_scatter)
%     sphere_size = 5e-3; % radius or the sphere
%     sphere_vol = 4*pi*sphere_size^3/3;
%     density = 15; % density = num_scatter/mm^3
%     num_scatters = floor(sphere_vol*1e9 * density);
%     
%     sphere_centers = [30e-3, 60e-3];
%     
%     x_pos = []; y_pos = []; z_pos = []; amp_scatt = [];
%     for sp_idx = 1:numel(sphere_centers)
%         rvals = 2*rand(num_scatters,1)-1;
%         elevation = asin(rvals);
%         azimuth = 2*pi*rand(num_scatters,1);
%         radii = sphere_size*(rand(num_scatters,1).^(1/3));
%         [x,y,z] = sph2cart(azimuth,elevation,radii);
%         z = z+sphere_centers(sp_idx);
%         
%         x_pos = cat(1,x_pos,x);
%         y_pos = cat(1,y_pos,y);
%         z_pos = cat(1,z_pos,z);
%         
%         amp_scatt = cat(1,amp_scatt,amp_high_reflect*rand(num_scatters, 1));
%     end
%     
%     m_scat_pos = cat(2, x_pos, y_pos, z_pos);

    % Whole region (prostate mimic)
    radius_prostate = 29e-3;
    center_z = 35e-3;
    center_y = 0;
    
    height_ = 2e-3;
    
    volume_ = (pi*radius_prostate^2)*height_;
    
    scatter_density = 15; % density = num_scatter/mm^3
    num_scatters = floor(volume_*1e9 * scatter_density);
   
    size_y = 70e-3;
    size_z = 70e-3;
    
    offset = 5e-3;
   
    x_pos = (rand(num_scatters,1)-0.5)*height_;
    y_pos = (rand(num_scatters,1)-0.5)*size_y;
    z_pos = rand(num_scatters,1)*(size_z-offset) + offset;
    
    amp = rand(num_scatters,1);
    
    % prostate
    inside_ = ((y_pos-center_y).^2 + (z_pos-center_z).^2 < radius_prostate^2);
    amp = amp .* inside_ + 1e-150;
    
    m_scat_pos = cat(2, x_pos, y_pos, z_pos);
    
    % cyst
    radius_cyst = 4e-3;
    deg_ = [-18 0 18];
    depth_cysts = [20 35 50]*1e-3;
    intensity_cysts = [5.6 2 1e-150]; % at each depth [dB] ( 15dB, 6dB, cyst)
    
    for d_idx = 1:numel(depth_cysts)
        depth_tmp = depth_cysts(d_idx);
        for k = 1:3
            lat_tmp = (radius + depth_tmp) * sind(deg_(k));
            ver_tmp = (radius + depth_tmp) * cosd(deg_(k)) - radius;
            
            inside_tmp = ((y_pos-lat_tmp).^2 + (z_pos-ver_tmp).^2 < radius_cyst^2);
            amp(inside_tmp) = intensity_cysts(k);
        end
    end
    amp_scatt = amp;
    
    % radially positioned
    radial_diff = [6 5 4 3 2 1]*1e-3;
    radial_dist = zeros(1,numel(radial_diff));
    for r_idx = 1:numel(radial_diff)
        radial_dist(r_idx) = sum(radial_diff(1:r_idx));
    end
    
    depth_pos = [12.5 27.5 42.5 57.5]*1e-3;
    
    deg_start = -35;
    deg_tmp = deg_start + (radial_dist/(depth_pos(1)+radius))*180/pi;
    
    lat_ = (radius + depth_pos(1)) * sind(deg_tmp);
    
    m_point_pos = [];
    for d_idx = 1:numel(depth_pos)
        depth_tmp = depth_pos(d_idx);
        
        deg_tmp = asind(lat_/(depth_tmp+radius));
        
        ver_tmp = (radius + depth_tmp) * cosd(deg_tmp) - radius;
        for k = 1:numel(radial_dist)
            m_point_pos = cat(1,m_point_pos, [0 lat_(k) ver_tmp(k)]);
        end
    end

    amp_point = amp_point_target*ones(size(m_point_pos,1),1);
    
    
%     axially positioned
%     m_point_pos = ... % x, y, z
%                 [0, 0, 10e-3;
%                  0, 0, 20e-3;
%                  0, 0, 30e-3;
%                  0, 0, 40e-3;
%                  0, 0, 50e-3;
%                  0, 0, 60e-3;
%                  0, 0, 70e-3;];
%     
%     amp_point = amp_point_target*ones(size(m_point_pos,1),1);
    
    
%     m_scat_total = cat(1, m_scat_pos, m_point_pos);
%     m_amp_total = cat(1, amp_scatt, amp_point);

%     m_scat_total = m_scat_pos;
%     m_amp_total = amp_scatt;
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
scatter(m_scat_total(:,2)*1e3, m_scat_total(:,3)*1e3,m_amp_total,'filled');
xlabel('Elevational'); ylabel('axial');
axis equal; axis tight; grid on;
ylim([0 70]);
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

pTxAperture = xdc_focused_array (stTRInfo.nNoEle, stTRInfo.nWidth, stTRInfo.nHeight, stTRInfo.nKerf, stTRInfo.nEleFocus, nSubDivNum_x, nSubDivNum_y, aFocalPoint);
pRxAperture = xdc_focused_array (stTRInfo.nNoEle, stTRInfo.nWidth, stTRInfo.nHeight, stTRInfo.nKerf, stTRInfo.nEleFocus, nSubDivNum_x, nSubDivNum_y, aFocalPoint);
% pTxAperture = xdc_linear_array (stTRInfo.nNoEle, stTRInfo.nWidth, stTRInfo.nHeight, stTRInfo.nKerf, nSubDivNum_x, nSubDivNum_y, aFocalPoint);
% pRxAperture = xdc_linear_array (stTRInfo.nNoEle, stTRInfo.nWidth, stTRInfo.nHeight, stTRInfo.nKerf, nSubDivNum_x, nSubDivNum_y, aFocalPoint);

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
stSaveInfo.nRadius = radius;
stSaveInfo.FOV_theta = FOV_theta;
stSaveInfo.aPosRng = aPosRng;

switch stTxInfo.type
    case 'conventional'
        folder_name = ['[Focused]Volume_LatFocus' num2str(stTxInfo.nTxFocus*1e3) 'mm_EleFocus' num2str(stTRInfo.nEleFocus*1e3) ...
            'mm_FOV' num2str(FOV_theta) '_Tx' num2str(stTxInfo.num_tx_vol)];
    case 'plane_wave'
        folder_name = ['[Plane_wave]Volume_Steer_' num2str(stTxInfo.max_steering) 'deg_' num2str(stTxInfo.num_tx_plane) 'compounding' ...
            '_FOV' num2str(FOV_theta) '_Tx' num2str(stTxInfo.num_tx_vol)];
end

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
mScatXYZPos_Rot = m_scat_total;

for p_idx = start_p_idx:stTxInfo.num_tx_vol
    disp(['>>> Volume scan [' num2str(p_idx) '/' num2str(stTxInfo.num_tx_vol) ']...']);
    % Field rotation
    nTheta = aPosRng(p_idx);
    aTheta_Zero = atand(m_scat_total(:,2)./(radius+m_scat_total(:,3)));
    aTheta_Prime = aTheta_Zero + nTheta;
    
    aPosZ = ((m_scat_total(:,3)+radius)./cosd(aTheta_Zero)) .* cosd(aTheta_Prime) - radius;
    aPosY = ((m_scat_total(:,3)+radius)./cosd(aTheta_Zero)) .* sind(aTheta_Prime);
    
    mScatXYZPos_Rot(:,3) = aPosZ;
    mScatXYZPos_Rot(:,2) = aPosY;
    % x-axis doesn't change
    
    vRcvData = zeros(nSample, stTRInfo.nNoEle, stTxInfo.num_tx_plane);
    dir_tmp = [save_dir_data '/RAW_' num2str(p_idx, '%.4d')];
    mkdir(dir_tmp);
    for t_idx = start_t_idx:stTxInfo.num_tx_plane
        disp(['    Transmit angle : ' num2str(stTxInfo.tx_angles(t_idx)) ' | [' num2str(t_idx) '/' num2str(stTxInfo.num_tx_plane) ']']);tic;
        
        apod_window = ones(1, stTRInfo.nNoEle); % for plane wave, all element transmit
        delay_curve = fieldII_get_delay_PW(t_idx, aElePosX, stRFInfo.nC, stTxInfo);
        
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
        save([dir_tmp '/Scanline_' num2str(t_idx, '%.3d')], 'mRcvData');
        toc;
    end
end

%%
field_end;