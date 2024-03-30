clc; clear all; close all;

bSave = 1;
%% Load
sMaster_path = uigetdir('', 'Select path');
aCase_list = dir(sMaster_path);

for cIdx = 3:numel(aCase_list)
    sCase_name = aCase_list(cIdx).name;
    [dir_f, folder_name, ext] = fileparts([sMaster_path '/' sCase_name]);
    if(strcmp(ext, '.DS_Store'))
        disp('>>> pass .DS_Store');
        continue
    end
    tic;
    disp(['>>> ' sCase_name '...']);
    sFolderPath = [sMaster_path '/' sCase_name];
    aFileList = dir(sFolderPath);
    for aIdx = 3:numel(aFileList)
        sFolderName = aFileList(aIdx).name;
        
        [dir_f, folder_name, ext] = fileparts([sFolderPath '/' sFolderName]);
        if(strcmp(ext, '.DS_Store'))
            disp('>>> pass .DS_Store');
            continue
        end
        disp(['>>> ' sFolderName '...']);
        sPath = [sFolderPath '/' sFolderName '/Saved Data'];
        aList = dir(sPath);
        %% Coordinate setting
        nNoPT = 9;
        nInitialPosX = 0e-3;
        nInitialPosZ = 10e-3;
        nSpacing = 10e-3;
        
        %% ROI position
        aTargetPosX = linspace(nInitialPosX, nInitialPosX, nNoPT);
        aTargetPosZ = linspace(nInitialPosZ, nInitialPosZ+(nNoPT-1)*nSpacing, nNoPT);
        
        nROISize_Hor = 60e-3;
        nROISize_Ver = 8e-3;
        
        mROIPosX = zeros(nNoPT, 2);
        mROIPosZ = zeros(nNoPT,2);
        
        mROIPosX(:,1) = aTargetPosX - 0.5*nROISize_Hor; % left
        mROIPosX(:,2) = aTargetPosX + 0.5*nROISize_Hor; % right
        
        mROIPosZ(:,1) = aTargetPosZ - 0.5*nROISize_Ver; % upper
        mROIPosZ(:,2) = aTargetPosZ + 0.5*nROISize_Ver; % bottom
        
        %% load data & resolution calculation
        adB = [-6 -6];
%         mLateralResol = zeros(numel(aList)-2,nNoPT);
%         mAxialResol = zeros(numel(aList)-2,nNoPT);
        flag = 0;
        for fIdx = 3:numel(aList)
            sFileName = aList(fIdx).name;
            [dir_f, folder_name, ext] = fileparts([sPath '/' sFileName]);
            if(strcmp(ext, '.DS_Store'))
                disp('>>> pass .DS_Store');
                flag = flag + 1;
                continue
            end
            disp(['>>> ' sFileName '...'])
            load([sPath '/' sFileName]);
            
%             stRFInfo = stSaveInfo.stRFInfo;
%             stTRInfo = stSaveInfo.stTRInfo;
%             stTxInfo = stSaveInfo.stTxInfo;
%             stBFInfo = stSaveInfo.stBFInfo;
            stRFInfo = stParam.stRFInfo;
            stTRInfo = stParam.stTRInfo;
            stTxInfo = stParam.stTxInfo;
            stBFInfo = stParam.stBFInfo;
            
            switch stBFInfo.sDirection
                case 'elevational'
                    mImgX = stParam.mImgY;
                    mImgZ = stParam.mImgZ;
                    aXAxis = stParam.aYAxis_DSC;
                    aZAxis = stParam.aZAxis_DSC;
                    mImg = stSaveInfo.mImg_DSC;
                case 'lateral'
            end
            
            
            
%             figure(fIdx);
            figure(1);
            for pIdx = 1:nNoPT
                nLft = mROIPosX(pIdx,1);
                nRgt = mROIPosX(pIdx,2);
                nUp = mROIPosZ(pIdx,1);
                nDn = mROIPosZ(pIdx,2);
                
                lIdx = find(abs(aXAxis-nLft)==min(abs(aXAxis-nLft)));
                rIdx = find(abs(aXAxis-nRgt)==min(abs(aXAxis-nRgt)));
                uIdx = find(abs(aZAxis-stBFInfo.nRadius-nUp)==min(abs(aZAxis-stBFInfo.nRadius-nUp)));
                dIdx = find(abs(aZAxis-stBFInfo.nRadius-nDn)==min(abs(aZAxis-stBFInfo.nRadius-nDn)));
                
                mROI_mag = mImg(uIdx:dIdx, lIdx:rIdx);
                
                % magnitude to db
                mROI = mag_to_db(mROI_mag);
                
                subplot(2,ceil(nNoPT/2),pIdx);
                imagesc(aXAxis(lIdx:rIdx)*1e3,(aZAxis(uIdx:dIdx)-stBFInfo.nRadius)*1e3,(mROI+abs(max(mROI(:))))); colormap gray; hold on; caxis([-80 0]);
                [ContourLine, h] = contour(aXAxis(lIdx:rIdx)*1e3,(aZAxis(uIdx:dIdx)-stBFInfo.nRadius)*1e3,(mROI+abs(max(mROI(:)))), adB,'ShowText','off', 'LineColor','r', 'LineWidth', 1);
                
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
%                     aY_Contour = ContourLine(2, (nMaxIdx+1):(nMaxIdx+nMaxVal));
                    nLatWidth = abs(max(aX_Contour)-min(aX_Contour));
%                     nAxlWidth = abs(max(aY_Contour)-min(aY_Contour));
                end
                
                text(1, (aTargetPosZ(pIdx)+1e-3)*1e3, [num2str(nLatWidth) 'mm'],'Color','y','FontSize',12);
                
                xlabel('Elevational [mm]'); ylabel('Axial [mm]'); title(['Contour image at' num2str(aTargetPosZ(pIdx)*1e2) 'cm']);
                set(gca,'Fontsize',18,'FontName','Times New Roman', 'FontWeight','bold');
                axis equal; axis tight;
                set(gcf,'Position',[818 45 1031 676]);
                
%                 mLateralResol(fIdx-2-flag,pIdx) = nLatWidth;
                  mLateralResol(pIdx) = nLatWidth;
%                 mAxialResol(fIdx-2-flag,pIdx) = nAxlWidth;
            end
        end
        
        %% Save
        if(bSave)
            sSaveDir = [sFolderPath '/' sFolderName '/Resolution'];
            mkdir(sSaveDir);
            save([sSaveDir '/Lateral resolution'],'mLateralResol');
%             save([sSaveDir '/Axial resolution'], 'mAxialResol');
        end
    end
    toc;
end
%%
disp('>>> All folder processed');
close all;