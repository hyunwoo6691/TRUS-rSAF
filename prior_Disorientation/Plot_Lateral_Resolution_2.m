%%% Plot for different delta-theta, same number of synthesize
clc; clear all; close all;

bCutoff = 1;
nCutoffVal = 20;

aNoSyn = [4 8 16 32 64 128];

aCheckSyn = [16];
%% Load resolution data
sFolderPath = uigetdir('','Choose path');
aList = dir(sFolderPath);

for fIdx = 3:numel(aList)
    [dir_f, folder_name, ext] = fileparts([sFolderPath '/' aList(fIdx).name]);
    if(strcmp(ext, '.DS_Store'))
        disp('>>> pass .DS_Store');
        continue
    end
    
    load([sFolderPath '/' aList(fIdx).name '/Resolution/Lateral resolution.mat']);
    for cIdx = 1:numel(aCheckSyn)
        Idx = find(aCheckSyn(cIdx) == aNoSyn);
        aLateralResol = mLateralResol(Idx, :);
        %% Plot
        nNoFile = size(mLateralResol,1);
        nNoPT = size(mLateralResol,2);
        
        nInitialDth = 10e-3;
        nStep = 10e-3;
        aZAxis = linspace(nInitialDth, nInitialDth+(nNoPT-1)*nStep, nNoPT);
        
        figure(cIdx);
        plot(aZAxis*1e3, aLateralResol,'LineWidth', 2, 'Marker', 'o'); hold on;
        xlabel('Depth [mm]'); ylabel('Lateral resolution [mm]'); grid on;
        title({'FWHM at different \Delta\theta',['Maximum number of synthesize: ' num2str(aCheckSyn(cIdx))]});
        if(bCutoff), ylim([0 nCutoffVal]);end
        leg = legend('Tx128', 'Tx256', 'Tx512','Location','northwest');
        set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
        set(gcf,'Position',[379 53 1500 676]);
    end
end
hold off;