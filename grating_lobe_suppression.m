clc; clear; close all;

addpath('functions');
%% set folder
dir_ = uigetdir('./data','');

folder_list = dir(dir_);
flag = 0;
if(strcmp('.DS_Store',folder_list(3).name)), flag = 1;end
folder_list = folder_list(3+flag:end);

%% load data
stParams = {};
beamforms_ = {};
for f_idx = 1:numel(folder_list)
    disp(['>>> ' folder_list(f_idx).name]);
    dir_tmp = [dir_ '/' folder_list(f_idx).name];
    
    load([dir_tmp '/Parameters.mat']);
    stParams = cat(3, stParams, stParam);
    
        % get optimized N_syn
        fid = fopen([dir_tmp '/FWHM/info.txt'], 'r');
        str = fgetl(fid);
        str = fgetl(fid);
        N_syn_ = split(str,': ');
        N_syn_ = str2double(N_syn_{end});
    
        file_name = ['[Save]SA_VS' num2str(stParam.stTRInfo.nEleFocus*1e3) '_Syn' num2str(N_syn_,'%.3d') '.mat'];
        disp(['    ' file_name]);
        load([dir_tmp '/errors_bf/error_0/Sample001/Element_64/' file_name]);
    
%     file_ = dir([dir_tmp '/errors_bf/error_0/Sample001/Element_64/']);
%     load([dir_tmp '/errors_bf/error_0/Sample001/Element_64/' file_(end).name]);
    
    beamforms_ = cat(3, beamforms_, stSaveInfo.mBFedData);
end
clear stParam stSaveInfo

%% mid & envelope detection
mid_.nTGC_Atten = 0.5;                % [dB]

mid_.nDCRType = 'high';
mid_.nDCRTap = 128;                   % BPF tap #
mid_.nDCRFcut = 1e6;

envs_ = {};
for f_idx = 1:numel(folder_list)
    beamformed_tmp = beamforms_{f_idx};
    param_tmp = stParams{f_idx};
    
    env_tmp = mid_proc(beamformed_tmp, mid_, param_tmp.stRFInfo, param_tmp.stBFInfo);
    envs_ = cat(3, envs_, env_tmp);
end

%%
FOV_ = 130;
radius_ = 10e-3;
view_depth_ = 72e-3;

%% DSC
dsc_data = [];
for f_idx = 1:numel(folder_list)
    env_tmp = envs_{f_idx};
    
    scanline_theta = linspace(-0.5*FOV_, 0.5*FOV_, size(env_tmp,2)); % Ground truth transmitted angle
    depth_ = linspace(radius_, radius_+view_depth_, 2432);
    
    da = abs(scanline_theta(1)-scanline_theta(2));
    dr = abs(depth_(1)-depth_(2));
    view_depth = view_depth_ + (radius_*(1-cosd(0.5*FOV_)));
    view_width = 2* (radius_+view_depth_)*sind(0.5*FOV_);
    
    dz = 1e-4;
    dy = 1e-4;
    height = round(view_depth / dz);
    width = round(view_width / dy);
    
    [axis_y, axis_z, dsc_tmp] = ScanConverter_convex(env_tmp, dr, da, radius_, height, width, dz, dy);
    dsc_data = cat(3, dsc_data, dsc_tmp);
end

dsc_avg = mean(dsc_data,3);

%% plot
for f_idx = 1:(numel(folder_list)+1)
    if(f_idx <= numel(folder_list))
        img_ = dsc_data(:,:,f_idx);
        title_ = folder_list(f_idx).name;
    else
        img_ = dsc_avg;
        title_ = 'averaged';
    end
    
    aROI = find(img_ ~= 50);
    aOutlier = find(img_ == 50);
    mOutput = zeros(size(img_));
    mOutput(aROI) = img_(aROI);
    mOutput_db = db(mOutput/max(mOutput(:)));
    mOutput_db(aOutlier) = -30;
    figure(f_idx);
    imagesc(axis_y*1e3,(axis_z-radius_)*1e3, mOutput_db); caxis([-80 0]);
    axis tight; axis equal;
    title(title_,'Interpreter','none');
    colormap gray; %colorbar;
    set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
    set(gcf,'Position',[300+700*(f_idx-1) 60 623 487]);
end









