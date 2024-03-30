clc; clear all; close all;

%%
dir_ = uigetdir('','');

load([dir_ '/Parameters.mat']);
stRFInfo = stParam.stRFInfo;
stBFInfo = stParam.stBFInfo;
stTRInfo = stParam.stTRInfo;
aYAxis_DSC = stParam.aYAxis_DSC;
aZAxis_DSC = stParam.aZAxis_DSC;
clear stParam

dir_bf = [dir_ '/Errors_bf'];

dir_gt = [dir_bf '/error_0'];

sigma_ = [0.1 0.2 0.5 1 2 5];

%% ground truth
disp('>>> ground truth data loading');
load([dir_gt '/Sample001/Element_64/[Save_elevational]Conv_Focus25.mat']);
conv_gt_bf = mBFedData;
clear mBFedData
load([dir_gt '/Sample001/Element_64/[Save_elevational]SA_Focus25.mat']);
rsaf_gt_bf = mBFedData;
clear mBFedData
%% error data
disp('>>> error data loading');
err_list = dir(dir_bf);
flag = 0;
if(strcmp(err_list(3).name, '.DS_Store')), flag = 1; end
err_list = err_list(4+flag:end);

num_sample = 100;

CON_ = {}; rSAF_ = {};
for e_idx = 1:numel(err_list)
    disp(['    ' err_list(e_idx).name]);
    sigma_tmp = err_list(e_idx).name;
    dir_err_tmp = [dir_bf '/' sigma_tmp];
    
    con_tmp = []; rsaf_tmp = [];
    for s_idx = 1:num_sample
        dir_err_tmpp = [dir_err_tmp '/Sample' num2str(s_idx,'%.3d') '/Element_64'];
        % conv
        load([dir_err_tmpp '/[Save_elevational]Conv_Focus25.mat']);
        con_tmp = cat(3, con_tmp, mBFedData);
        clear mBFedData
        
        % rSAF
        load([dir_err_tmpp '/[Save_elevational]SA_Focus25.mat']);
        rsaf_tmp = cat(3, rsaf_tmp, mBFedData);
        clear mBFedData
    end
    
    CON_ = cat(3, CON_, con_tmp);
    rSAF_ = cat(3, rSAF_, rsaf_tmp);
end
clear con_tmp rsaf_tmp
%% parameter setting
fov_ = stBFInfo.nFOV;
scline = stBFInfo.nScline;
radius = stBFInfo.nRadius;
bf_depth = 72e-3;
bf_depth_spl = ceil(bf_depth/stRFInfo.nUnitDis*2);

scanline_deg = linspace(-0.5*fov_, 0.5*fov_, scline);
d_theta = abs(scanline_deg(1)-scanline_deg(2));
aDth = linspace(radius, radius+bf_depth, bf_depth_spl);

stMID.nTGC_Atten = 0.5;                % [dB]
stMID.nDCRType = 'high';
stMID.nDCRTap = 128;                   % BPF tap #
stMID.nDCRFcut = 1e6;

da = abs(scanline_deg(1)-scanline_deg(2));
dr = abs(aDth(1)-aDth(2));
nViewDth = stBFInfo.nDth + (stBFInfo.nRadius*(1-cosd(0.5*stBFInfo.nFOV)));
nViewWidth = 2* (stBFInfo.nRadius+stBFInfo.nDth)*sind(0.5*stBFInfo.nFOV);

dz = 1e-4;
dy = 1e-4;
nHeight = round(nViewDth / dz);
nWidth = round(nViewWidth / dy);

%% processing for ground truth data
% CON envelope
[filter_out, spatial_filter] = spatial_filtering(conv_gt_bf, 1, d_theta, 'none'); % 'gauss', 'boxcar', 'none'
[dcr_out, dcr_filter] = DCR(filter_out, stMID, stRFInfo);
[tgc_out, tgc_curve] = fDTGC(dcr_out, stMID, stRFInfo, stBFInfo, size(filter_out,1), stRFInfo.nUnitDis);
conv_gt_env = abs(hilbert(tgc_out));

% rSAF envelope
[filter_out, spatial_filter] = spatial_filtering(rsaf_gt_bf, 1, d_theta, 'none'); % 'gauss', 'boxcar', 'none'
[dcr_out, dcr_filter] = DCR(filter_out, stMID, stRFInfo);
[tgc_out, tgc_curve] = fDTGC(dcr_out, stMID, stRFInfo, stBFInfo, size(filter_out,1), stRFInfo.nUnitDis);
rsaf_gt_env = abs(hilbert(tgc_out));

% CON envelope with filtering
[filter_out, spatial_filter] = spatial_filtering(conv_gt_bf, 1, d_theta, 'boxcar'); % 'gauss', 'boxcar', 'none'
[dcr_out, dcr_filter] = DCR(filter_out, stMID, stRFInfo);
[tgc_out, tgc_curve] = fDTGC(dcr_out, stMID, stRFInfo, stBFInfo, size(filter_out,1), stRFInfo.nUnitDis);
conv_gt_fil_env = abs(hilbert(tgc_out));

%% processing for error data
disp('>>> filtering');
CON_env = {}; CON_env_fil = {}; rSAF_env = {};
spatial_filters = {};

for e_idx = 1:numel(err_list)
    disp(['    ' err_list(e_idx).name]);
    con_tmp = CON_{e_idx};
    rsaf_tmp = rSAF_{e_idx};
    
    con_env_tmp = []; con_env_fil_tmp = []; rsaf_env_tmp = [];
    con_dsc_tmp = []; con_dsc_fil_tmp = []; rsaf_dsc_tmp = [];
    for s_idx = 1:num_sample
        % CON envelope
        [filter_out, spatial_filter] = spatial_filtering(con_tmp(:,:,s_idx), sigma_(e_idx), d_theta, 'none'); % 'gauss', 'boxcar', 'none'
        [dcr_out, dcr_filter] = DCR(filter_out, stMID, stRFInfo);
        [tgc_out, tgc_curve] = fDTGC(dcr_out, stMID, stRFInfo, stBFInfo, size(filter_out,1), stRFInfo.nUnitDis);
        backend_out = abs(hilbert(tgc_out));
        con_env_tmp = cat(3, con_env_tmp, backend_out);
        
        % rSAF envelope
        [filter_out, spatial_filter] = spatial_filtering(rsaf_tmp(:,:,s_idx), sigma_(e_idx), d_theta, 'none'); % 'gauss', 'boxcar', 'none'
        [dcr_out, dcr_filter] = DCR(filter_out, stMID, stRFInfo);
        [tgc_out, tgc_curve] = fDTGC(dcr_out, stMID, stRFInfo, stBFInfo, size(filter_out,1), stRFInfo.nUnitDis);
        backend_out = abs(hilbert(tgc_out));
        rsaf_env_tmp = cat(3, rsaf_env_tmp, backend_out);
        
        % CON envelope with filtering
        [filter_out, spatial_filter] = spatial_filtering(con_tmp(:,:,s_idx), sigma_(e_idx), d_theta, 'gauss'); % 'gauss', 'boxcar', 'none'
        [dcr_out, dcr_filter] = DCR(filter_out, stMID, stRFInfo);
        [tgc_out, tgc_curve] = fDTGC(dcr_out, stMID, stRFInfo, stBFInfo, size(filter_out,1), stRFInfo.nUnitDis);
        backend_out = abs(hilbert(tgc_out));
        con_env_fil_tmp = cat(3, con_env_fil_tmp, backend_out);
    end
    CON_env = cat(3, CON_env, con_env_tmp); rSAF_env = cat(3, rSAF_env, rsaf_env_tmp); CON_env_fil = cat(3, CON_env_fil, con_env_fil_tmp);
    spatial_filters = cat(3, spatial_filters, spatial_filter);
    
end
clear con_env_tmp rsaf_env_tmp con_env_fil_tmp con_tmp rsaf_tmp con_dsc_tmp rsaf_dsc_tmp con_dsc_fil_tmp

%%
figure(1);
imagesc(db(rsaf_gt_env/max(rsaf_gt_env(:))));

%% plot for ground truth
% z_idx = 1015; % 3cm
% z_idx = 1692; % 5cm
z_idx = 2029; % 6cm

% for ground truth
con_gt_1d = conv_gt_env(z_idx,:);
con_gt_fil_1d = conv_gt_fil_env(z_idx,:);
rsaf_gt_1d = rsaf_gt_env(z_idx,:);

axis_ = linspace(-0.5*fov_, 0.5*fov_, numel(con_gt_1d));

figure(111);
plot(axis_, (con_gt_1d/max(rsaf_gt_1d))); hold on;
plot(axis_, (con_gt_fil_1d/max(rsaf_gt_1d)));
plot(axis_, (rsaf_gt_1d/max(rsaf_gt_1d))); hold off;
legend('CON', 'CON fil', 'rSAF');

%% plot for error data
for s_idx = 4:numel(sigma_)
    con_env_tmp = CON_env{s_idx};
    con_env_fil_tmp = CON_env_fil{s_idx};
    rsaf_env_tmp = rSAF_env{s_idx};
    
    for spl_idx = 1:num_sample
        con_env_1d_tmp = con_env_tmp(z_idx,:,spl_idx);
        con_env_fil_1d_tmp = con_env_fil_tmp(z_idx,:,spl_idx);
        rsaf_env_1d_tmp = rsaf_env_tmp(z_idx,:,spl_idx);
        
        figure(s_idx);
        plot(axis_, con_env_1d_tmp,'color', [0.12 0.12 0.12]*8); hold on;
        plot(axis_, con_env_fil_1d_tmp,'color', [0.12 0.12 0.12]*7);
        plot(axis_, rsaf_env_1d_tmp,'color', [0.12 0.12 0.12]*6);
    end
    
    plot(axis_, mean(squeeze(con_env_tmp(z_idx,:,:)),2), 'color', 'k', 'LineWidth', 2);
    plot(axis_, mean(squeeze(con_env_fil_tmp(z_idx,:,:)),2), 'color', 'b', 'LineWidth', 2);
    plot(axis_, mean(squeeze(rsaf_env_tmp(z_idx,:,:)),2), 'color', 'r', 'LineWidth', 2); hold off;
end







