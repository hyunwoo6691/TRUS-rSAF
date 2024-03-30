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
%%
start_idx = input('Enter the start index: ');

end_idx = input('Enter the end index: ');
%% Set parameter
stRFInfo.nFc = 10e6;
stRFInfo.nC = 1540;
stRFInfo.nFs = 4*stRFInfo.nFc;
stRFInfo.nAttenFactor = 0.5/(1e-2*1e6); % Attenuation factor[dB/(m*Hz)] % 0.5 dB/(cm*MHz)
stRFInfo.nUnitDis = stRFInfo.nC / stRFInfo.nFs;
stRFInfo.nLambda = stRFInfo.nC / stRFInfo.nFc;

stTRInfo.nNoEle                 = 1;
stTRInfo.nDiameter = 2e-3;
stTRInfo.nRadius = stTRInfo.nDiameter/2;
stTRInfo.f_num = 1;
stTRInfo.nEleFocus            = 0;

stTxInfo.cycle = 1;
stTxInfo.num_tx_vol = 280*3+1; % volume scan - half theta (0.2362)

aElePosX = 0;
aElePosZ = 0;

nTc = 1 / stRFInfo.nFc;
nUpsampleFactor = 5;
nSample = 3000; % 7.26 cm view depth

nSimulFreq = stRFInfo.nFs * nUpsampleFactor;
dT = 1 /nSimulFreq;

aFocalPoint = [0 0 0];
%% Rotate angle
%%%%%%%%%%%%
radius = 20e-3;

aPosRng = linspace(0.5*(stTxInfo.num_tx_vol-1), -0.5*(stTxInfo.num_tx_vol-1), stTxInfo.num_tx_vol) * 0.1564;
FOV_theta = abs(aPosRng(1)-aPosRng(end));
%%%%%%%%%%%%

%% Phantom geometry
view_depth = 72e-3;

margin = 0;
% amp_high_reflect = 20;
amp_point_target = 20;

%% Phantom - sphere object
if(~load_scatter)
    % radially positioned
    radial_pos = [-1 1; % 10mm
        -1.5 1.5; % 20mm
        -2 2; % 30mm
        -3 3; % 40mm
        -4 4; % 50mm
        -5 5; % 60mm
        -6 6]*1e-3; % 70mm
    %     radial_pos = [0; 0; 0; 0; 0; 0; 0];
    depth_pos = [10 20 30 40 50 60 70]*1e-3;

    %     m_point_pos = []; degs_ = zeros(1,numel(depth_pos));
    %     for d_idx = 1:numel(depth_pos)
    %         depth_tmp = depth_pos(d_idx);
    %
    %         m_point_pos = cat(1,m_point_pos, [-10e-3 0 depth_pos(d_idx)]);
    %         m_point_pos = cat(1,m_point_pos, [0 3e-3 depth_pos(d_idx)]);
    %         m_point_pos = cat(1,m_point_pos, [10e-3 6e-3 depth_pos(d_idx)]);
    %
    %     end

    m_point_pos = [];
    depth_pos = [10 20 30 40 50 60 70]*1e-3;
    for dIdx = 1:numel(depth_pos)
        m_point_pos = cat(1,m_point_pos, [0 0 depth_pos(dIdx)]);
    end
    %
    amp_point = amp_point_target*ones(size(m_point_pos,1),1);
    %
    m_scat_total = m_point_pos;
    m_amp_total = amp_point;


    % Cysts

    % Visualize
    figure(3);
    scatter3(m_point_pos(:,1)*1e3, m_point_pos(:,2)*1e3, m_point_pos(:,3)*1e3,amp_point,'filled');
    xlabel('Lateral'); ylabel('Elevational'); zlabel('Axial');
    axis equal; axis tight; %grid on;

    %     m_scat_total = cat(1, m_scat_pos, m_point_pos);
    %     m_amp_total = cat(1, amp_scatt, amp_point);

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

figure(4);
scatter(m_scat_total(:,2)*1e3, m_scat_total(:,3)*1e3,m_amp_total,'filled');
xlabel('Elevational'); ylabel('axial');
axis equal; axis tight; grid on;
ylim([0 70]); %xlim([-20 20]);
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
nFarFieldDepth_x = 0.01e-3; % depth over nFarFieldDepth_x is assumed to be far field on x-z plane
% For accurate directivity pattern, the number of source is more than 5
nFarFieldDepth_y = 0.01e-3;     % depth over nFarFieldDepth_x is assumed to be far field on y-z plane

nMathEleSize_x = sqrt(nFarFieldDepth_x * 4 * stRFInfo.nLambda); % H < sqrt(4*lambda*z)
nMathEleSize_y = sqrt(nFarFieldDepth_y * 4 * stRFInfo.nLambda);

nSubDivNum_x = ceil(stTRInfo.nDiameter / nMathEleSize_x);
nSubDivNum_y = ceil(stTRInfo.nDiameter / nMathEleSize_y);

pTxAperture = xdc_piston(stTRInfo.nRadius, nMathEleSize_x);
pRxAperture = xdc_piston(stTRInfo.nRadius, nMathEleSize_x);

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


folder_name = ['cis_ii_data'];

save_dir_parameter = [save_dir_master folder_name];
save_dir_data = [save_dir_parameter '/Raw data'];
mkdir(save_dir_parameter);
mkdir(save_dir_data);

save([save_dir_parameter '/stSaveInfo'], 'stSaveInfo');


%% Loop for pressure calculation at each position
mScatXYZPos_Rot = m_scat_total;

mRcvData = zeros(nSample, stTxInfo.num_tx_vol);
for p_idx = start_idx:end_idx
    tic;
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

    dir_tmp = [save_dir_data '/RAW_' num2str(p_idx, '%.4d')];
    mkdir(dir_tmp);

    %%% Calculate the received signals from a collection of scatterers for all elements
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
    mRcvDataTmp = mRF_deci;

    mRcvData(:,p_idx) = mRcvDataTmp;

    % save scanline data
    save([dir_tmp '/raw_data'], 'mRcvDataTmp');
    toc;
end

%%
field_end;

currentTime = clock;
disp(['>>> Finished at ' num2str(currentTime(4)) ':' num2str(currentTime(5),'%.2d')]);