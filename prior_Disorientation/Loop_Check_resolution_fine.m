clc; clear; close all;

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

dB = [-6 -6];

%% set path & load parameter file

dir_master = uigetdir('/01.Data', '');
folder_list = dir(dir_master);
flag = 0;
if(strcmp(folder_list(3).name, '.DS_Store'))
    flag = 1;
end
folder_list = folder_list(3+flag:end);

%% process one by one
for mf_idx = 1:numel(folder_list)
% for mf_idx = folder_start:numel(folder_list)
    folder_name = folder_list(mf_idx).name;
    
    disp(['>>> Folder name: ' folder_name]);
    dir_ = [dir_master '/' folder_name];
    
    dir_save = [dir_ '/Resolutions_' num2str(dB(1)) 'dB'];
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
    error_list = error_list(3+flag:end);
    
%     for e_idx = 1:numel(error_list)
    for e_idx = numel(error_list)-1:numel(error_list)
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
        
        already_done_list = dir(dir_savee);
        
        if(numel(already_done_list)>2) % if there exists already processed files
            if(strcmp(already_done_list(end).name, 'Sample100'))
                continue;
            end
            last_done = already_done_list(end).name;
            last_done = split(last_done, '.');
            last_done = last_done{1};
            last_done = str2double(last_done(end-2:end));
        else
            last_done = 1;
        end
        
        for s_idx = last_done:numel(sample_list)
            disp(['    Sample # : ' num2str(s_idx)]); tic;
            
            sample_name = sample_list(s_idx).name;
            
            dir_tmpp = [dir_tmp '/' sample_name];
            file_list = dir(dir_tmpp);
            flag = 0;
            if(strcmp(file_list(3).name, '.DS_Store'))
                flag = 1;
            end
            file_list = file_list(3+flag:end);
            
            mLateralResol = zeros(numel(file_list), nNoPT);
            for f_idx = 1:numel(file_list)
                file_name = file_list(f_idx).name;
                
                disp(['        File : ' file_name]);
                
                load([dir_tmpp '/' file_name]);
                mImg = stSaveInfo.mImg_DSC;
                clear stSaveInfo;
                
                figure(f_idx)
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
                    idx_roi = find(mROI ~= -30);
                    mROI(idx_roi) = mROI(idx_roi) + abs(max(mROI(idx_roi)));
                    
                    subplot(2,ceil(nNoPT/2),p_idx); % db scale image
                    imagesc(aYAxis_DSC(lIdx:rIdx)*1e3,(aZAxis_DSC(uIdx:dIdx)-stBFInfo.nRadius)*1e3,mROI); colormap gray; hold on; caxis([-80 0]);
                    [ContourLine, h] = contour(aYAxis_DSC(lIdx:rIdx)*1e3,(aZAxis_DSC(uIdx:dIdx)-stBFInfo.nRadius)*1e3,mROI, dB,'ShowText','off', 'LineColor','r', 'LineWidth', 1);
                    
                    % detect circle
                    nCurrentVal = ContourLine(2,1);
                    nCurrentIdx = 1;
                    nColumn = size(ContourLine,2);
                    
                    k = 1;
                    mCircleIdx = zeros(2,nColumn);
                    mCircleIdx(1,k) = nCurrentIdx;
                    mCircleIdx(2,k) = nCurrentVal;
                    
                    if(nCurrentVal + nCurrentIdx == nColumn)
                        nNoVertices = ContourLine(2,1);
                        aX_Contour = ContourLine(1,2:(2+nNoVertices-1));
                        aY_Contour = ContourLine(2,2:(2+nNoVertices-1));
                        nLatWidth = abs(max(aX_Contour) - min(aX_Contour));
                        nAxlWidth = abs(max(aY_Contour)-min(aY_Contour));
                    else
                        nPrevVal = nCurrentVal;
                        nPreIdx = nCurrentIdx;
                        while(nCurrentIdx+nCurrentVal < nColumn)
                            k = k+1;
                            nNextIdx = nCurrentIdx +1 + nCurrentVal;
                            nNext = ContourLine(2, nNextIdx);
                            nPreVal = nCurrentVal;
                            nPreIdx = nCurrentIdx;
                            nCurrentIdx = nNextIdx;
                            nCurrentVal = nNext;
                            
                            mCircleIdx(1,k) = nCurrentIdx;
                            mCircleIdx(2,k) = nCurrentVal;
                        end
                        nMaxVal = max(mCircleIdx(2,:));
                        nMaxIdx = mCircleIdx(1,find(mCircleIdx(2,:)==nMaxVal));
                        
                        aX_Contour = ContourLine(1, (nMaxIdx+1):(nMaxIdx+nMaxVal));
                        aY_Contour = ContourLine(2, (nMaxIdx+1):(nMaxIdx+nMaxVal));
                        nLatWidth = abs(max(aX_Contour)-min(aX_Contour));
                        nAxlWidth = abs(max(aY_Contour)-min(aY_Contour));
                    end
                    
%                     text(1, (aTargetPosZ(p_idx)+1e-3)*1e3, [num2str(nLatWidth) 'mm'],'Color','y','FontSize',12);
                    
%                     xlabel('Elevational [mm]'); ylabel('Axial [mm]'); title(['Contour image at' num2str(aTargetPosZ(p_idx)*1e2) 'cm']);
%                     set(gca,'Fontsize',18,'FontName','Times New Roman', 'FontWeight','bold');
%                     axis equal; axis tight;
%                     set(gcf,'Position',[818 45 1031 676]);
                    
                    mLateralResol(f_idx,p_idx) = nLatWidth;
                end
            end
            % save
            if(bSave)
                save([dir_savee '/' sample_name],'mLateralResol');
            end
            elapsed_time = toc;
            disp(['    Elapsed time: ' num2str(elapsed_time) ' sec']);
            disp('');
%             if(elapsed_time > 5)
%                 disp('=================== Rebooting script... ===================');
%                 disp(['Current folder : ' folder_name ', index: ' num2str(folder_start)]);
%                 disp(['  Error case : ' case_name ', index: ' num2str(e_idx)]);
%                 stVars.dir_master = dir_master;
% %                 stVars.e_start = e_idx;
%                 stVars.folder_start = mf_idx;
%                 terminate_script
%             end
        end
    end
end
disp('>>> All process completed');

