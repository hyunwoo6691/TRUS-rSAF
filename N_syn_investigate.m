clc;
%%
acoustic_ = stParam.stRFInfo;
bf_ = stParam.stBFInfo;
trans_ = stParam.stTRInfo;

imgpos_y = stParam.mImgY;
imgpos_z = stParam.mImgZ;

clear stParam;

% mid processing paramete
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
dirTmp = uigetdir('./data');

dir_ = [dirTmp '/Sample001/Element_64/'];
fileList = dir(dir_);
flag = 0;
if(strcmp(fileList(3).name,'.DS_Store')), flag = 1; end
fileList = fileList(3+flag:end);

errorCase = split(dirTmp,'errors_bf/');
errorCase = errorCase{2};
disp(errorCase);

dirSaveTmp = '/Users/songhyunwoo/Desktop/disorientation_paper/figures/prostate_mimic/OPT_NsynCheck';

dirSave = [dirSaveTmp '/' errorCase];
mkdir(dirSave);

for f_idx = 1:numel(fileList)
    close all;
    fileTmp = fileList(f_idx).name;
    
    disp(['    ' fileTmp]);
    
    numSyn_ = split(fileTmp,'Syn');
    numSyn_ = split(numSyn_{end},'.');
    numSyn_ = str2double(numSyn_{1});
    
    figName_ = ['Nsyn' num2str(numSyn_,'%.3d')];
    
    %%
    load([dir_ fileTmp]);
    %%
    beamformed_data = stSaveInfo.mBFedData;
    % beamformed_data = mBFedData;
    % beamformed_data = InterpNan(beamformed_data);
    env_data = mid_proc(beamformed_data, mid_, acoustic_, bf_);
    
    [axis_y, axis_z, dsc_data] = ScanConverter_convex(env_data, dr, da, bf_.nRadius, height, width, dz, dy);
    % [axis_y, axis_z, dsc_data] = ScanConverter_convex(totalNumSyn_, dr, da, bf_.nRadius, height, width, dz, dy);
    
    aROI = find(dsc_data ~= 50);
    aOutlier = find(dsc_data == 50);
    mOutput = zeros(size(dsc_data));
    mOutput(aROI) = dsc_data(aROI);
    mOutput_db = db(mOutput/max(mOutput(:)));
    mOutput_db(aOutlier) = -30;
    
    figure(1321);
    imagesc(axis_y*1e3,(axis_z-bf_.nRadius)*1e3, mOutput_db); caxis([-50 0]);
    % imagesc(axis_y*1e3,(axis_z-bf_.nRadius)*1e3, dsc_data); caxis([0 131]);
    axis tight; axis equal;
    % xlabel('Elevational [mm]'); ylabel('Axial [mm]');
    colormap gray; %colorbar;
    set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
    % set(gcf,'Position',[818 45 1031 676]);
    set(gcf,'Position',[100 100 330 392]);
    ylim([0 70]); xlim([-35 35]);
    
    %% save
    savefig([dirSave '/' figName_ '.fig']);
end

%%
scline = 254;
synDepth_ = zeros(numel(fileList),bf_.nDthSpl);
%%
numSynTot = zeros(1, numel(fileList));

for f_idx = 1:numel(fileList)
    close all;
    fileTmp = fileList(f_idx).name;
    
    disp(['    ' fileTmp]);
    
    numSyn_ = split(fileTmp,'Syn');
    numSyn_ = split(numSyn_{end},'.');
    numSyn_ = str2double(numSyn_{1});
    
    numSynTot(f_idx) = numSyn_;
end
%%
for f_idx = 1:numel(fileList)
    close all;
    fileTmp = fileList(f_idx).name;
    
    disp(['    ' fileTmp]);
    
    numSyn_ = split(fileTmp,'Syn');
    numSyn_ = split(numSyn_{end},'.');
    numSyn_ = str2double(numSyn_{1});
    
    %%
    load([dir_ fileTmp]);
    %%
    vSynReg_ = stSaveInfo.vSynReg;
    
    synReg_ = sum(vSynReg_(:,:,max(scline-0.5*(numSyn_-1),1):min(scline+0.5*(numSyn_-1), size(vSynReg_,3))),3);
    
    numSynDepth_ = synReg_(:,scline);
    
    synDepth_(f_idx,:) = numSynDepth_;
end

%%
axisZ = imgpos_z(:,scline);
figure(1);
subplot(2,1,1);
imagesc((axisZ-bf_.nRadius)*1e3,numSynTot, synDepth_);
xlabel('depth [mm]'); ylabel('N_{syn}');

ratio_ = synDepth_ ./ repmat(synDepth_(end, :), numel(numSynTot), 1);
subplot(2,1,2);
imagesc((axisZ-bf_.nRadius)*1e3,numSynTot, ratio_);
xlabel('depth [mm]'); ylabel('N_{syn}'); title('ratio');
















