% clc; clear; close all;

bSave = 1;

%% Coordinate setting
nNoPT = 10;
nInitialPosY = 0e-3;
nInitialPosZ = 10e-3;
nSpacing = 10e-3;

%% ROI position
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

%% set path & load parameter file

% dir_master = uigetdir('/01.Data', '');
dir_master = [dir_master '/Data'];

folder_list = dir(dir_master);

flag = 0;
if(strcmp(folder_list(3).name, '.DS_Store'))
    flag = 1;
end

folder_list = folder_list(3+flag:end);

%% process one by one
% for mf_idx = 1:numel(folder_list)
for mf_idx = folder_start:numel(folder_list)
    folder_name = folder_list(mf_idx).name;
    
    disp(['>>> Folder name: ' folder_name]);
    dir_ = [dir_master '/' folder_name];
    
    dir_save = [dir_ '/XCORR_center_wSA'];
    mkdir(dir_save);
    
    load([dir_ '/Parameters.mat']);
    
    stBFInfo = stParam.stBFInfo;
    aYAxis_DSC = stParam.aYAxis_DSC;
    aZAxis_DSC = stParam.aZAxis_DSC;
    clear stParam;
    
    %%
    dir_error = [dir_ '/Errors'];
    error_list = dir(dir_error);
    
    flag = 0;
    if(strcmp(error_list(3).name, '.DS_Store'))
        flag = 1;
    end
    error_list = error_list(4+flag:end);
    ground_truth = error_list(3+flag).name;
    
    dir_gt = [dir_error '/' ground_truth '/Sample001'];
    gt_list = dir(dir_gt);
    flag = 0;
    if(strcmp(gt_list(3).name, '.DS_Store')), flag = 1; end
    gt_list = gt_list(3+flag:end);
    
    for e_idx = error_start:numel(error_list)
        case_name = error_list(e_idx).name;
        dir_savee = [dir_save '/' case_name];
        mkdir(dir_savee);
        
        dir_tmp = [dir_error '/' case_name];
        sample_list = dir(dir_tmp);
        flag = 0;
        if(strcmp(sample_list(3).name, '.DS_Store'))
            flag = 1;
        end
        sample_list = sample_list(3+flag:end);
        
%         already_done_list = dir(dir_savee);
%         
%         if(numel(already_done_list)>2) % if there exists already processed files
%             if(strcmp(already_done_list(end).name, 'Sample100'))
%                 continue;
%             end
%             last_done = already_done_list(end).name;
%             last_done = split(last_done, '.');
%             last_done = last_done{1};
%             last_done = str2double(last_done(end-2:end));
%         else
%             last_done = 1;
%         end
        
        for s_idx = sample_start:numel(sample_list)
            disp(['    Sample # : ' num2str(s_idx)]); tic;
            
            sample_name = sample_list(s_idx).name;
            
            dir_tmpp = [dir_tmp '/' sample_name];
            file_list = dir(dir_tmpp);
            flag = 0;
            if(strcmp(file_list(3).name, '.DS_Store'))
                flag = 1;
            end
            file_list = file_list(3+flag:end);
            
            xcorr_val = zeros(numel(file_list), nNoPT);
            for f_idx = 1:numel(file_list)
                file_name = file_list(f_idx).name;
%                 file_name_gt = gt_list(f_idx).name;
                file_name_gt = gt_list(end).name; % only SA 
                
                disp(['        File : ' file_name]);
                
                load([dir_tmpp '/' file_name]);
                mImg = stSaveInfo.mImg_DSC;
                
                % load ground truth image
                load([dir_gt '/' file_name_gt]);
                mImg_gt = stSaveInfo.mImg_DSC;
                
                clear stSaveInfo;
                
                for p_idx = 1:nNoPT
                    nLft = mROIPosY(p_idx,1);
                    nRgt = mROIPosY(p_idx,2);
                    nUp = mROIPosZ(p_idx,1);
                    nDn = mROIPosZ(p_idx,2);
                    
                    lIdx = find(abs(aYAxis_DSC-nLft)==min(abs(aYAxis_DSC-nLft)));
                    rIdx = find(abs(aYAxis_DSC-nRgt)==min(abs(aYAxis_DSC-nRgt)));
                    uIdx = find(abs(aZAxis_DSC-stBFInfo.nRadius-nUp)==min(abs(aZAxis_DSC-stBFInfo.nRadius-nUp)));
                    dIdx = find(abs(aZAxis_DSC-stBFInfo.nRadius-nDn)==min(abs(aZAxis_DSC-stBFInfo.nRadius-nDn)));
                    
                    mROI = mImg(uIdx:dIdx, lIdx:rIdx);
                    
                    mROI_gt = mImg_gt(uIdx:dIdx, lIdx:rIdx);
                    
                    correlation_map = normxcorr2(mROI_gt, mROI);
                    
%                     max_xcorr_val(f_idx, p_idx) = max(correlation_map(:));
                    center_coord = ceil(size(correlation_map)/2);
                    xcorr_val(f_idx, p_idx) = correlation_map(center_coord(1), center_coord(2));
                end
            end
            % save
            if(bSave)
                save([dir_savee '/' sample_name],'xcorr_val');
            end
            elapsed_time = toc;
            disp(['    Elapsed time: ' num2str(elapsed_time) ' sec']);
            disp('');
            if(elapsed_time > 5)
                disp('=================== Rebooting script... ===================');
                disp(['Current folder : ' folder_name ', index: ' num2str(folder_start)]);
                disp(['  Error case : ' case_name ', index: ' num2str(e_idx)]);
                stVars.dir_master = dir_master;
                stVars.sample_start = s_idx;
                stVars.error_start = e_idx;
                stVars.folder_start = mf_idx;
                terminate_script
            end
        end
        sample_start = 1;
    end
    error_start = 1;
end
disp('>>> All process completed');

