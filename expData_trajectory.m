clc; clear; close all;

addpath('functions');
%%
dirMaster_ = uigetdir('./data','');

rList = dir(dirMaster_);
flag = 0;
if(strcmp(rList(3).name, '.DS_Store')), flag = 1; end
rList = rList(3+flag:end);
%%
posGT_ = [];
posGTStage_ = [];
scanningAngle_ = [];
for rIdx = 1:numel(rList)
    %%
    dir_ = [dirMaster_ '/' rList(rIdx).name];
    %%
    load([dir_ '/rRotate.mat']);
    load([dir_ '/rTR.mat']);
    
    load([dir_ '/origin.mat']);
    origin_(1:2) = origin_(1:2)*1e-3;
    
    load([dir_ '/scanningAngle.mat']);
    %% ground truth trajectory
    posGTTmp = [0; -(rRotate-rTR)] + (rRotate-rTR)*[sind(scanningAngle(1,:)); cosd(scanningAngle(1,:))];
    posGTStageTmp = [0; -(rRotate-rTR)] + (rRotate-rTR+rTR)*[sind(scanningAngle(1,:)); cosd(scanningAngle(1,:))];
    
    %%
    posGT_ = cat(3, posGT_, posGTTmp);
    posGTStage_ = cat(3, posGTStage_, posGTStageTmp);
    scanningAngle_ = cat(3, scanningAngle_, scanningAngle);
    %%
    figure(rIdx);
    subplot(3,1,1);
    plot(posGTTmp(1,:),'LineWidth',2,'color','k');
    
    subplot(3,1,2);
    plot(posGTTmp(2,:),'LineWidth',2,'color','k');
    
    subplot(3,1,3);
    plot(scanningAngle(1,:),'LineWidth',2,'color','k');
end

%%
for rIdx = 1:numel(rList)
    figure(510);
    plot(posGT_(1,:,rIdx)*1e3, (posGT_(2,:,rIdx)-rTR)*1e3,'LineWidth',2,'color',[0.2 0.2 0.2]*rIdx, 'LineStyle',':'); hold on;
    plot(posGTStage_(1,:,rIdx)*1e3, (posGTStage_(2,:,rIdx)-rTR)*1e3,'LineWidth',2,'color',[0.2 0.2 0.2]*rIdx);
end
line([-30 30], [0 0],'LineWidth',1,'LineStyle','--','color','k');
line([0 0], [-30 30],'LineWidth',1,'LineStyle','--','color','k'); hold off;
axis equal; grid on;
xlim([-10 10]); ylim([-15 15]);
xlabel('Frontal (mm)'); ylabel('Sagittal (mm)');

figure(511);
plot(scanningAngle_(1,:,rIdx),'LineWidth',2,'color','k', 'LineStyle',':');
xlabel('Frontal (mm)'); ylabel('Angle (deg)'); grid on; hold off;

