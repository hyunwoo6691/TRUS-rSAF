%%% Plot for same delta-theta, different number of synthesize
clc; clear all; close all;

bCutoff = 0;
nCutoffVal = 20;
%% Load resolution data
sFolderPath = uigetdir('','Choose path');

load([sFolderPath '/Resolution/Lateral resolution.mat']);
%% Plot
nNoFile = size(mLateralResol,1);
nNoPT = size(mLateralResol,2);

nInitialDth = 10e-3;
nStep = 10e-3;
aZAxis = linspace(nInitialDth, nInitialDth+(nNoPT-1)*nStep, nNoPT);

figure(1);
plot(aZAxis*1e3, mLateralResol(1,:),'LineWidth', 2, 'Marker', 'o', 'LineStyle','--'); hold on;
for fIdx = 2:nNoFile
%     plot(aZAxis*1e3, mLateralResol(fIdx,:),'LineWidth',2,'LineStyle','-.','Marker','*');
    plot(aZAxis*1e3, mLateralResol(fIdx,:),'LineWidth',2,'Marker','o'); hold on;
end
hold off; %ylim([0 20]);
xlabel('Depth [mm]'); ylabel('Lateral resolution [mm]'); grid on;
title({'FWHM at different number of synthesize'});
if(bCutoff), ylim([0 nCutoffVal]);end
% leg = legend('Conv', 'Syn4', 'Syn8', 'Syn16', 'Syn32', 'Syn64','Location','northwest');
leg = legend('Conv', 'Syn4', 'Syn8', 'Syn16', 'Syn32', 'Location','northwest');
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf,'Position',[379 53 1500 676]);
