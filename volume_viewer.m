clc; clear; close all;

addpath('functions');
%%
dir_ = uigetdir('./data','');

load([dir_ '/Parameters.mat']);
load([dir_ '/stSaveInfo.mat']);

dirData_ = [dir_ '/errors_bf/VolumetricBeamforming'];
planeList = dir(dirData_);
flag = 0;
if(strcmp(planeList(3).name, '.DS_Store')), flag = 1; end
planeList = planeList(3+flag:end);
%%
stRFInfo = stParam.stRFInfo;
stBFInfo = stParam.stBFInfo;
stTRInfo = stParam.stTRInfo;
stTxInfo = stParam.stTxInfo;

vImgX = stParam.vImgX;
vImgY = stParam.vImgY;
vImgZ = stParam.vImgZ;

clear stParam;

mid_.nTGC_Atten = 0.5;                % [dB]
mid_.nDCRType = 'bandpass';
mid_.nDCRTap = 128;                   % BPF tap #
mid_.nDCRFcut = 1e6;
mid_.nDCRF1 = 3.9e6;
mid_.nDCRF2 = 9.1e6;

aDth = linspace(stBFInfo.nRadius, stBFInfo.nRadius+stBFInfo.nDth, stBFInfo.nDthSpl); % axial
aPosRng_fieldII = fliplr(stSaveInfo.aPosRng);
clear stSaveInfo;

%% load beamformed data & mid processing
stBFInfo.nScline = stBFInfo.nLatScline;
vBEedData = [];

disp('>>> Data loading and mid processing');
for scIdx = 1:stTxInfo.num_tx_vol
    if(mod(scIdx, stTxInfo.num_tx_vol/4)==0)
        disp(['    ' num2str(round(100*scIdx/stTxInfo.num_tx_vol)) '%...']);
    end
    
    planeTmp = planeList(scIdx).name;
    load([dirData_ '/' planeTmp]);
    backendTmp = mid_proc(mBFedData, mid_, stRFInfo, stBFInfo);
    vBEedData = cat(3, vBEedData, backendTmp);
end
%%
offset_ = 0.7404*1e-3;
offsetS_ = round(offset_ / stRFInfo.nUnitDis * 2);

method_ = split(dir_,'/data/');
method_ = method_{2};

%% get specific plane - lateral-axial
deg_ = 0; % degree 3, -3, -9
planeIdx = find(abs(aPosRng_fieldII-deg_)==min(abs(aPosRng_fieldII-deg_)));
planeIdx = planeIdx(1);
% planeIdx = 133;

axisX = squeeze(vImgX(1,:,planeIdx));

if(strcmp('REF',method_))
    planeTmp = vBEedData(offsetS_:end,:,planeIdx);
    aDthOffset = aDth(offsetS_:end);
    
else
    planeTmp = vBEedData(:,:,planeIdx);
end
normVal = max(planeTmp(:));

% DSC
disp('>>> dsc');
heightLat = 75e-3;
widthLat = 40e-3;
dz = 1e-4;
dx = 1e-4;
nHeight = round(heightLat / dz);
nWidth = round(widthLat / dx);

[aXAxis, aZAxis, mOutput] = ScanConverter_linear(planeTmp/normVal, stTRInfo.nPitch, stRFInfo.nUnitDis, nHeight, nWidth, dz, dx);

depth_ = 70e-3 - offset_;
lat_ = -10e-3;

depthIdx = find(abs(aZAxis-depth_)==min(abs(aZAxis-depth_)));
latIdx = find(abs(aXAxis-lat_)==min(abs(aXAxis-lat_)));
roiSize = round([10/dx 5/dz]*1e-3);

roiZaxis = depthIdx-round(0.5*roiSize(2)):depthIdx+round(0.5*roiSize(2));
roiXaxis = latIdx-round(0.5*roiSize(1)):latIdx+round(0.5*roiSize(1));
roi_ = mOutput(roiZaxis,roiXaxis);

figure(1);
imagesc(aXAxis*1e3, (aZAxis)*1e3, db(mOutput/max(mOutput(:))));
colormap gray; axis equal; axis tight;
caxis([-50 0]);
%     ylim([0 75]);
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf,'Position',[100 100 330 392]);

figure(3124);
contourf(aXAxis(roiXaxis)*1e3, (aZAxis(roiZaxis))*1e3,db(roi_/max(roi_(:)))); hold on;
[ContourLine, h]=contour(aXAxis(roiXaxis)*1e3, (aZAxis(roiZaxis))*1e3,db(roi_/max(roi_(:))),[-6 -6],'LineWidth',2,'color','k','LineStyle',':'); hold off;
axis equal; axis tight;
set(gca,'Ydir','reverse');
set(gcf,'Position',[462 99 560 420]);
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
[FWHM_Lat, aX_Contour]= get_contour(ContourLine);

%% get 1-d profile
depths_ = [30 50 70]*1e-3;

profile_ = [];
for dIdx = 1:numel(depths_)
    depth_ = depths_(dIdx) - offset_;
    depthIdx = find(abs(aZAxis-depth_)==min(abs(aZAxis-depth_)));
    
    profile_ = cat(1, profile_, mOutput(depthIdx,:));
end

profile_N = profile_ / max(profile_(:));

figure(314);
for dIdx = 1:numel(depths_)
    plot(aXAxis*1e3, (profile_N(dIdx,:)),'LineWidth',2); hold on;
end
hold off;

%% get specific plane - elevation-axial

% dsc parameter for elevational-axial plane
scanline_theta = linspace(-0.5*stBFInfo.nFOV, 0.5*stBFInfo.nFOV, stBFInfo.nRadialScline); % Ground truth transmitted angle
depth_ = linspace(stBFInfo.nRadius, stBFInfo.nRadius+stBFInfo.nDth, stBFInfo.nDthSpl);

da = abs(scanline_theta(1)-scanline_theta(2));
dr = abs(depth_(1)-depth_(2));
view_depth = stBFInfo.nDth + (stBFInfo.nRadius*(1-cosd(0.5*stBFInfo.nFOV)));
view_width = 2* (stBFInfo.nRadius+stBFInfo.nDth)*sind(0.5*stBFInfo.nFOV);

dz = 1e-4;
dy = 1e-4;
heightRad = round(view_depth / dz);
widthRad = round(view_width / dy);

% ele_ = 64;
ele_ = 42;
% ele_ = 41;
if(strcmp('REF',method_))
    planeTmp = squeeze(vBEedData(offsetS_:end,ele_,:));
    aDthOffset = aDth(offsetS_:end);
    
else
    planeTmp = squeeze(vBEedData(:,ele_,:));
end

[axis_y, axis_z, dsc_data] = ScanConverter_convex(planeTmp, dr, da, stBFInfo.nRadius, heightRad, widthRad, dz, dy);

aROI = find(dsc_data ~= 50);
aOutlier = find(dsc_data == 50);
mOutput = zeros(size(dsc_data));
mOutput(aROI) = dsc_data(aROI);

normVal = max(mOutput(:));

mOutput_db = db(mOutput/normVal);
mOutput_db(aOutlier) = -30;

depth_ = 30e-3 - offset_;
depthIdx = find(abs(axis_z-stBFInfo.nRadius-depth_)==min(abs(axis_z-stBFInfo.nRadius-depth_)));
roiSize = round([10/dy 5/dz]*1e-3);

roiZaxis = depthIdx-round(0.5*roiSize(2)):depthIdx+round(0.5*roiSize(2));
roiYaxis = round(0.5*(numel(axis_y)-roiSize(1))):round(0.5*(numel(axis_y)+roiSize(1)));
roi_ = dsc_data(roiZaxis,roiYaxis);

figure(1321);
imagesc(axis_y*1e3,(axis_z-stBFInfo.nRadius)*1e3, mOutput_db); 
caxis([-50 0]);
axis tight; axis equal; 
colormap gray; %colorbar;
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf,'Position',[100 100 330 392]);
xlim([-35 35]); ylim([0 75]);

figure(3125);
contourf(axis_y(roiYaxis)*1e3, (axis_z(roiZaxis)-stBFInfo.nRadius)*1e3,db(roi_/max(roi_(:)))); hold on;
[ContourLine, h]=contour(axis_y(roiYaxis)*1e3, (axis_z(roiZaxis)-stBFInfo.nRadius)*1e3,db(roi_/max(roi_(:))),[-6 -6],'LineWidth',2,'color','b','LineStyle',':'); hold off;
set(gca,'Ydir','reverse'); axis tight; axis equal;
set(gcf,'Position',[462 99 560 420]);
[FWHM_rad, aX_Contour]= get_contour(ContourLine);

%% Volume rendering (stack of 2d (elevational-axial) interpolation)
vDSCVolume = zeros(heightRad,widthRad,stTRInfo.nNoEle);
vBEedDataNorm = vBEedData/max(vBEedData(:));
for eIdx = 1:stTRInfo.nNoEle
    planeTmp = squeeze(vBEedDataNorm(:,eIdx,:));
    
    [axis_y, axis_z, dsc_data] = ScanConverter_convex(planeTmp, dr, da, stBFInfo.nRadius, heightRad, widthRad, dz, dy);
    
    vDSCVolume(:,:,eIdx) = dsc_data;
end

aROI = find(vDSCVolume ~= 50);
aOutlier = find(vDSCVolume == 50);
vOutput = zeros(size(vDSCVolume));
vOutput(aROI) = vDSCVolume(aROI);

vOutput_db = db(vOutput/max(vOutput(:)));
vOutput_db(aOutlier) = -30;

%%
%elevation - axial 
mipYZ = squeeze(max(vOutput,[],[3]));

figure;
imagesc(axis_y*1e3, (axis_z-stBFInfo.nRadius)*1e3, db(mipYZ/max(mipYZ(:)))); colormap gray;
axis equal; axis tight; caxis([-50 0]);
set(gcf,'Position',[100 100 330 392]);
xlim([-35 35]); ylim([0 75]);

%% Volume rendering (3d interpolation)
% 
% if(strcmp('REF',method_))
%     grid_x = vImgX(offsetS_:end,:,1);
% else
%     grid_x = vImgX(:,:,1);
% end
% grid_y = zeros(size(grid_x));
% grid_z = repmat(aDthOffset',1,stTRInfo.nNoEle);
% 
% % coordinate for acquired data
% volume_grid_x = repmat(grid_x, 1, 1, numel(aPosRng_fieldII));
% volume_grid_y = zeros(size(volume_grid_x));
% volume_grid_z = zeros(size(volume_grid_x));
% 
% for o_idx = 1:numel(aPosRng_fieldII)
%     theta_ = aPosRng_fieldII(o_idx);
%     
%     theta_zero = atand(grid_y./(stBFInfo.nRadius+grid_z));
%     theta_prime = theta_zero + theta_;
%     
%     volume_grid_z(:,:,o_idx) = ((grid_z+stBFInfo.nRadius)./cosd(theta_zero)) .* cosd(theta_prime) - stBFInfo.nRadius;
%     volume_grid_y(:,:,o_idx) = ((grid_z+stBFInfo.nRadius)./cosd(theta_zero)) .* sind(theta_prime);
% end
% 
% % dsc parameter
% scanline_theta = aPosRng_fieldII; % Ground truth transmitted angle
% 
% viewDepth = stBFInfo.nDth + (stBFInfo.nRadius*(1-cosd(0.5*stBFInfo.nFOV)));
% viewEleWidth = 2* (stBFInfo.nRadius+stBFInfo.nDth)*sind(0.5*stBFInfo.nFOV);
% viewLatWidth = abs(vImgX(1,1,1) - vImgX(1,end,1));
% 
% dx = 1e-4;
% dy = 1e-4;
% dz = 1e-4;
% 
% volume_size = [viewLatWidth, 40e-3, viewDepth]; % x,y,z
% voxel_size = [dx, dy, dz]; % x,y,z
% 
% [dscOut,axisX, axisY, axisZ,rejectedPixels] = dsc_volume(vBEedData/max(vBEedData(:)), volume_size, voxel_size, aPosRng_fieldII, volume_grid_x, volume_grid_z,stBFInfo.nRadius);
% 
% axisZ_ = axisZ - stBFInfo.nRadius;

%% MIP
% elevation - axial 
% mipYZ = squeeze(max(dscOut,[],[1]))';
% 
% figure;
% imagesc(axisY*1e3, (axisZ-stBFInfo.nRadius)*1e3, db(mipYZ/max(mipYZ(:)))); colormap gray;
% axis equal; axis tight; caxis([-50 0]);
% set(gcf,'Position',[100 100 330 392]);
% xlim([-20 20]); ylim([0 75]);

%% PSF
% roiSize = [10e-3 10e-3 10e-3];
% roiSizeSample = round(roiSize ./ voxel_size);
% 
% depth_ = 30e-3;
% rangeZ = depth_ + [-roiSize(3) roiSize(3)]*0.5;
% rangeZIdx = (find(abs(axisZ_-rangeZ(1))==min(abs(axisZ_-rangeZ(1)))) : ... 
%             find(abs(axisZ_-rangeZ(2))==min(abs(axisZ_-rangeZ(2)))));
%         
% lat_ = -10e-3;
% rangeX = lat_ + [-roiSize(1) roiSize(1)]*0.5;
% rangeXIdx = (find(abs(axisX-rangeX(1))==min(abs(axisX-rangeX(1)))) : ...
%             find(abs(axisX-rangeX(2))==min(abs(axisX-rangeX(2)))));
% 
% 
% ele_ = 0;
% rangeY = ele_ + [-roiSize(2) roiSize(2)]*0.5;
% rangeYIdx = (find(abs(axisY-rangeY(1))==min(abs(axisY-rangeY(1)))) : ... 
%             find(abs(axisY-rangeY(2))==min(abs(axisY-rangeY(2)))));
% 
% 
% axisXROI = axisX(rangeXIdx);
% axisYROI = axisY(rangeYIdx);
% axisZROI = axisZ_(rangeZIdx);
% volROI = dscOut(rangeXIdx, rangeYIdx, rangeZIdx);
% volROI = volROI/max(volROI(:));
% volROIdB = db(volROI);
% 
% load('config_Volshow.mat');
% figure(3210);
% v = volshow(volROIdB,config);
% 
% % contour cutoff
% dB_ = -80;
% volROIdB_cutoff = volROIdB;
% volROIdB_cutoff(volROIdB_cutoff<dB_) = -120;
% figure(3322);
% v = volshow(volROIdB_cutoff,config);
