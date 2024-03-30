% 2021.01.20
% script for compute the spatial resolution of point target at different
% focal depth.
clc; clear; close all;

addpath('./functions');
b_save = 1;
%%
dir_master = uigetdir('./data','');
folder_list = dir(dir_master);
flag = 0;
if(strcmp(folder_list(3).name,'.DS_Store')), flag = 1; end
folder_list = folder_list(3+flag:end);

for folder_idx = 1:numel(folder_list)
    folder_tmp = folder_list(folder_idx).name;
    
    disp(folder_tmp);
    
    %% set directory
    dir_ = [dir_master '/' folder_tmp];
    
    dir_tmp = [dir_ '/errors_bf/error_0/Sample001/Element_64'];
    
    file_list = dir(dir_tmp);
    flag = 0;
    if(strcmp(file_list(3).name,'.DS_Store')), flag = 1; end
    file_list = file_list(3+flag:end);
    
    dir_save = [dir_ '/FWHM'];
    mkdir(dir_save);
    %% cluster beamforming type
    bf_type = {}; num_syn = [];
    for f_idx = 1:numel(file_list)
        file_tmp = split(file_list(f_idx).name,'_');
        file_tmpp = split(file_tmp{1}, ']');
        switch file_tmpp{2}
            case 'Conv'
                bf_type = cat(1, bf_type, 'CON');
            case 'SA'
                bf_type = cat(1, bf_type, 'rSAF');
                
                syn_tmp = split(file_tmp{3},'.');
                syn_tmpp = split(syn_tmp{1}, 'n');
                num_syn = cat(1, num_syn, str2num(syn_tmpp{2}));
        end
    end
    
    %% load parameter
    load([dir_ '/Parameters.mat']);
    
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
    
    %% set region of measurement
    % Coordinate setting
    nNoPT = 7;
    nInitialPosY = 0e-3;
    nInitialPosZ = 10e-3;
    nSpacing = 10e-3;
    
    % ROI position
    aTargetPosY = linspace(nInitialPosY, nInitialPosY, nNoPT);
    aTargetPosZ = linspace(nInitialPosZ, nInitialPosZ+(nNoPT-1)*nSpacing, nNoPT);
    
    nROISize_Hor = 60e-3;
    nROISize_Ver = 8e-3;
    
    mROIPosY = zeros(nNoPT, 2);
    mROIPosZ = zeros(nNoPT,2);
    
    mROIPosY(:,1) = aTargetPosY - 0.5*nROISize_Hor; % left
    mROIPosY(:,2) = aTargetPosY + 0.5*nROISize_Hor; % right
    
    mROIPosZ(:,1) = aTargetPosZ - 0.5*nROISize_Ver; % upper
    mROIPosZ(:,2) = aTargetPosZ + 0.5*nROISize_Ver; % bottom
    
    dB = [-6 -6];
    % dB = [-12 -12];
    
    %% load beamformed data
    rSAF_ = {}; rSAF_syn = {};
    for f_idx = 1:numel(file_list)
        load([dir_tmp '/' file_list(f_idx).name]);
        switch bf_type{f_idx}
            case 'CON'
                CON_ = stSaveInfo.mBFedData;
            case 'rSAF'
                rSAF_ = cat(3, rSAF_, stSaveInfo.mBFedData);
                rSAF_syn = cat(3, rSAF_syn, stSaveInfo.vSynReg);
        end
    end
    
    %% do mid processing
    % for CON
    CON_env = mid_proc(CON_, mid_, acoustic_, bf_);
    % for rSAF
    rSAF_env = {};
    for f_idx = 1:numel(rSAF_)
        rsaf_tmp = mid_proc(rSAF_{f_idx}, mid_, acoustic_, bf_);
        rSAF_env = cat(3, rSAF_env, rsaf_tmp);
    end
    clear rsaf_tmp
    
    %% do dsc
    % for CON
    [axis_y, axis_z, CON_dsc] = ScanConverter_convex(CON_env, dr, da, bf_.nRadius, height, width, dz, dy);
    % for rSAF
    rSAF_dsc = {};
    for f_idx = 1:numel(rSAF_)
        [axis_y, axis_z, rsaf_tmp] = ScanConverter_convex(rSAF_env{f_idx}, dr, da, bf_.nRadius, height, width, dz, dy);
        rSAF_dsc = cat(3, rSAF_dsc, rsaf_tmp);
    end
    clear rsaf_tmp
    
    %% show image
    % CON
    aROI = find(CON_dsc ~= 50);
    aOutlier = find(CON_dsc == 50);
    mOutput = zeros(size(CON_dsc));
    mOutput(aROI) = CON_dsc(aROI);
    mOutput_db = db(mOutput/max(mOutput(:)));
    mOutput_db(aOutlier) = -30;
    figure(100);
    imagesc(axis_y*1e3,(axis_z-bf_.nRadius)*1e3, mOutput_db); caxis([-55 0]);
    axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]');colormap gray; colorbar;
    title('CON');
    set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
    set(gcf,'Position',[818 45 1031 676]);
    
    % rSAF
    for f_idx = 1:numel(rSAF_)
        aROI = find(rSAF_dsc{f_idx} ~= 50);
        aOutlier = find(rSAF_dsc{f_idx} == 50);
        mOutput = zeros(size(rSAF_dsc{f_idx}));
        mOutput(aROI) = rSAF_dsc{f_idx}(aROI);
        mOutput_db = db(mOutput/max(mOutput(:)));
        mOutput_db(aOutlier) = -30;
        figure(100+num_syn(f_idx));
        imagesc(axis_y*1e3,(axis_z-bf_.nRadius)*1e3, mOutput_db); caxis([-55 0]);
        axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]');colormap gray; colorbar;
        title(file_list(f_idx+1).name, 'Interpreter', 'none');
        set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
        set(gcf,'Position',[818 45 1031 676]);
    end
    %%
    close all;
    %% measure the spatial resolution
    % for CON
    FWHM_CON = zeros(1, nNoPT);
    
    for p_idx = 1:nNoPT
        nLft = mROIPosY(p_idx,1);
        nRgt = mROIPosY(p_idx,2);
        nUp = mROIPosZ(p_idx,1);
        nDn = mROIPosZ(p_idx,2);
        
        lIdx = find(abs(axis_y-nLft)==min(abs(axis_y-nLft)));
        rIdx = find(abs(axis_y-nRgt)==min(abs(axis_y-nRgt)));
        uIdx = find(abs(axis_z-bf_.nRadius-nUp)==min(abs(axis_z-bf_.nRadius-nUp)));
        dIdx = find(abs(axis_z-bf_.nRadius-nDn)==min(abs(axis_z-bf_.nRadius-nDn)));
        
        roi_ = CON_dsc(uIdx:dIdx, lIdx:rIdx);
        roi_db = zeros(size(roi_));
        idx_roi = find(roi_ ~= 50);
        idx_reject = find(roi_ == 50);
        roi_db(idx_roi) = db(roi_(idx_roi)/max(roi_(idx_roi)));
        roi_db(idx_reject) = -30;
        
        resol_tmp = measure_spatial_resolution(roi_db,axis_y(lIdx:rIdx)*1e3, (axis_z(uIdx:dIdx)-bf_.nRadius)*1e3, dB, 1, nNoPT, p_idx);
        
        FWHM_CON(p_idx) = resol_tmp;
    end
    figure(1);
    set(gcf, 'Position', [743 82 575 971]);
    % for rSAF
    FWHM_rSAF = zeros(numel(rSAF_), nNoPT);
    
    for f_idx = 1:numel(rSAF_)
        rSAF_tmp = rSAF_dsc{f_idx};
        for p_idx = 1:nNoPT
            nLft = mROIPosY(p_idx,1);
            nRgt = mROIPosY(p_idx,2);
            nUp = mROIPosZ(p_idx,1);
            nDn = mROIPosZ(p_idx,2);
            
            lIdx = find(abs(axis_y-nLft)==min(abs(axis_y-nLft)));
            rIdx = find(abs(axis_y-nRgt)==min(abs(axis_y-nRgt)));
            uIdx = find(abs(axis_z-bf_.nRadius-nUp)==min(abs(axis_z-bf_.nRadius-nUp)));
            dIdx = find(abs(axis_z-bf_.nRadius-nDn)==min(abs(axis_z-bf_.nRadius-nDn)));
            
            roi_ = rSAF_tmp(uIdx:dIdx, lIdx:rIdx);
            roi_db = zeros(size(roi_));
            idx_roi = find(roi_ ~= 50);
            idx_reject = find(roi_ == 50);
            roi_db(idx_roi) = db(roi_(idx_roi)/max(roi_(idx_roi)));
            roi_db(idx_reject) = -30;
            
            resol_tmp = measure_spatial_resolution(roi_db,axis_y(lIdx:rIdx)*1e3, (axis_z(uIdx:dIdx)-bf_.nRadius)*1e3,dB, num_syn(f_idx), nNoPT, p_idx);
            
            FWHM_rSAF(f_idx, p_idx) = resol_tmp;
        end
        figure(num_syn(f_idx));
        set(gcf, 'Position', [743 82 575 971]);
    end
    
    %% pick the best num_syn
    % at each depth, pick the best one
    % get the average num_syn across all depth
    
    indices_ = zeros(1, nNoPT);
    for p_idx = 1:nNoPT
        resol_tmp = FWHM_rSAF(:,p_idx);
        idx_tmp = find(resol_tmp == min(resol_tmp));
        indices_(p_idx) = idx_tmp(1); % it is good to choose the less num_syn
    end
    best_num_syn= max(num_syn(indices_));
    
    %% save
    if(b_save)
        fileID = fopen([dir_save '/info.txt'],'w');
        fprintf(fileID, 'Unit: mm \n');
        fprintf(fileID, ['Best num_syn: ' num2str(best_num_syn) '\n']);
        fprintf(fileID, 'File names \n');
        for f_idx = 1:numel(file_list)
            fprintf(fileID, [file_list(f_idx).name '\n']);
        end
        % CON
        save([dir_save '/FWHM_CON.mat'],'FWHM_CON');
        % rSAF
        for f_idx = 1:numel(rSAF_)
            save([dir_save '/FWHM_rSAF.mat'],'FWHM_rSAF');
        end
    end
end