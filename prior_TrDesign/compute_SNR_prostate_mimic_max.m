clear; close all;
addpath('functions');
clc;

%%
% N_syn_ = 45:4:155;
N_syn_ = [1];
%%
% gt_posY = [-7.6 -2.74 1.25 4.21 6.18 7.13]*1e-3;
% gt_posZ = [11.2 12.3 12.5 12.1 11.6 11.3;
%            26.7 27.4 27.5 27.3 27.0 26.8;
%            41.9 42.4 42.5 42.3 42.1 42.0;
%            57.1 57.4 57.5 57.4 57.2 57.1]*1e-3;
       
gt_posY = [-7.6 -2.74 1.25]*1e-3;
gt_posZ = [11.2 12.3 12.5;
           26.7 27.4 27.5;
           41.9 42.4 42.5;
           57.1 57.4 57.5]*1e-3;
       
roi_posY = repmat(gt_posY,4,1);
roi_posZ = gt_posZ;

% depth_roi = [20 35 50]*1e-3;
% deg = 0;

% depth_roi = [20 35 50]*1e-3;
roi_size = 2e-3;
%%
dir_ = uigetdir('./data',''); % choose the master folder . e.g., rsaf_optim_fele
folder_list = dir(dir_);
flag = 0;
if(strcmp(folder_list(3).name, '.DS_Store')), flag = 1; end
folder_list = folder_list(3+flag:end);

%%
for folder_idx = 1:numel(folder_list)
    %% set path 
    folder_name = folder_list(folder_idx).name;
    disp(['>>> ' folder_name ]);
    dir_tmp = [dir_ '/' folder_name];
    
    %% set parameter
    load([dir_tmp '/Parameters.mat']);
    
    acoustic_ = stParam.stRFInfo;
    bf_ = stParam.stBFInfo;
    trans_ = stParam.stTRInfo;
    
    imgpos_y = stParam.mImgY;
    imgpos_z = stParam.mImgZ;
    
    clear stParam;
    
    % mid processing parameter
    mid_.nTGC_Atten = 0.5;                % [dB]
    
    mid_.nDCRType = 'high';
    mid_.nDCRTap = 128;                   % BPF tap #
    mid_.nDCRFcut = 1e6;
    
    % dsc parameter
    scanline_theta = linspace(-0.5*bf_.nFOV, 0.5*bf_.nFOV, bf_.nScline); % Ground truth transmitted angle
    depth_ = linspace(bf_.nRadius, bf_.nRadius+bf_.nDth, bf_.nDthSpl);
    
    da = abs(scanline_theta(1)-scanline_theta(2));
    dr = abs(depth_(1)-depth_(2));
    view_depth = bf_.nDth + (bf_.nRadius*(1-cosd(0.5*bf_.nFOV)));
    view_width = 2* (bf_.nRadius+bf_.nDth)*sind(0.5*bf_.nFOV);
    
    dz = 1e-4;
    dy = 1e-4;
    height = round(view_depth / dz);
    width = round(view_width / dy);
    
    %%
%     roi_posZ = (depth_roi+bf_.nRadius)*cosd(deg)-bf_.nRadius;
%     roi_posY = (depth_roi+bf_.nRadius)*sind(deg);
    
    %% get data files & noise files
    dir_save = [dir_tmp '/SNR'];
    mkdir(dir_save);
    
    dir_signal = [dir_tmp '/errors_bf/error_0/'];
    dir_noise = [dir_tmp '/errors_bf_AWGN_noTarget/error_0/'];
    
    data_list_sample = dir(dir_signal);
    flag = 0;
    if(strcmp(data_list_sample(3).name, '.DS_Store')), flag = 1; end
    data_list_sample = data_list_sample(3+flag:end);
    
    noise_list_sample = dir(dir_noise);
    flag = 0;
    if(strcmp(noise_list_sample(3).name,'.DS_Store')), flag = 1; end
    noise_list_sample = noise_list_sample(3+flag:end);
    
    %%
%     num_syn = num_syn_c{folder_idx};
    num_syn = N_syn_(folder_idx);
    %% iterate
    snr_ = zeros(numel(data_list_sample)*numel(noise_list_sample),numel(roi_posZ),numel(num_syn));
    noise_ = zeros(numel(data_list_sample)*numel(noise_list_sample),numel(roi_posZ),numel(num_syn));
    signal_ = zeros(numel(data_list_sample)*numel(noise_list_sample),numel(roi_posZ),numel(num_syn));
    for data_idx = 1:numel(data_list_sample) % SampleXXX
        data_sample_idx = data_list_sample(data_idx).name;
        disp(['    Data name: ' data_sample_idx]);
        
        dir_data_tmpp = [dir_signal '/' data_sample_idx '/Element_64'];
        
        data_list_cases = dir(dir_data_tmpp);
        flag = 0;
        if(strcmp(data_list_cases(3).name,'.DS_Store')), flag = 1; end
        data_list_cases = data_list_cases(3+flag:end);
        
        for case_idx = 1:numel(data_list_cases) % num_syn
            case_ = data_list_cases(case_idx).name;
            disp(['      Case: ' case_]);
            load([dir_data_tmpp '/' case_]);
            
            beamformed_data = stSaveInfo.mBFedData;
            env_data = mid_proc(beamformed_data, mid_, acoustic_, bf_);
            [axis_y, axis_z, dsc_data] = ScanConverter_convex(env_data, dr, da, bf_.nRadius, height, width, dz, dy);
            
            for noise_idx = 1:numel(noise_list_sample) % SampleXXX
                noise_sample_idx = noise_list_sample(noise_idx).name;
                disp(['        Noise name: ' noise_sample_idx]);
                
                dir_noise_tmpp = [dir_noise '/' noise_sample_idx '/Element_64'];
                load([dir_noise_tmpp '/' case_]);
                
                beamformed_noise = stSaveInfo.mBFedData;
                env_noise = mid_proc(beamformed_noise, mid_, acoustic_, bf_);
                [axis_y, axis_z, dsc_noise] = ScanConverter_convex(env_noise, dr, da, bf_.nRadius, height, width, dz, dy);
                
                %% compute SNR at each depth
                for d_idx = 1:numel(roi_posZ)
%                     dth_ = depth_roi(d_idx);
%                     hor_rng = [-roi_size roi_size];
%                     ver_rng = [-roi_size roi_size] + dth_;
                    hor_rng = [-roi_size roi_size] + roi_posY(d_idx);
                    ver_rng = [-roi_size roi_size] + roi_posZ(d_idx);
                    
                    hor_idx = [find(abs(hor_rng(1)-axis_y)==min(abs(hor_rng(1)-axis_y))) ...
                        find(abs(hor_rng(2)-axis_y)==min(abs(hor_rng(2)-axis_y)))];
                    ver_idx = [find(abs(ver_rng(1)-axis_z+bf_.nRadius)==min(abs(ver_rng(1)-axis_z+bf_.nRadius))) ...
                        find(abs(ver_rng(2)-axis_z+bf_.nRadius)==min(abs(ver_rng(2)-axis_z+bf_.nRadius)))];
                    
                    roi_data = dsc_data(ver_idx(1):ver_idx(2), hor_idx(1):hor_idx(2));
                    roi_noise = dsc_noise(ver_idx(1):ver_idx(2), hor_idx(1):hor_idx(2));
                    
                    if(folder_idx == 1), figure(101+d_idx); imagesc(roi_data); pause(0.1); end
                    
                    snr_(numel(data_list_sample)*(data_idx-1)+noise_idx,...
                            d_idx,...
                            case_idx) = ...
                                        mean((roi_data(:)).^2)/var(roi_noise(:));
                                    
                    noise_(numel(data_list_sample)*(data_idx-1)+noise_idx,...
                            d_idx,...
                            case_idx) = var(roi_noise(:));
                        
                    signal_(numel(data_list_sample)*(data_idx-1)+noise_idx,...
                            d_idx,...
                            case_idx) = mean((roi_data(:)).^2);
                end
            end
        end
    end
    save([dir_save '/snr_var.mat'], 'snr_');
    save([dir_save '/signal_mean_noNoise.mat'], 'signal_');
    save([dir_save '/noise_var.mat'], 'noise_');
    
end

disp('all folder processed');







