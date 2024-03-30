clc; close all;

waveLength = 1540 / 6.5e6;

dTheta = 0.4724 * pi / 180;

radius_ = [5 10 15] * 1e-3;

z_vs = 5e-3;

depth_ = linspace(0, 70e-3, 100);

scanlineInterval = dTheta*(radius_ + z_vs) / waveLength; % unit: lambda
figure;
plot(radius_*1e3, scanlineInterval,'LineWidth',2,'color','k');

scanlineInterval = dTheta*depth_ / waveLength; % unit: lambda
figure;
plot(depth_*1e3, scanlineInterval,'LineWidth',2,'color','k');


scanlineInterval  = [];
for rIdx = 1:numel(radius_)
    scIntTmp = dTheta*(radius_(rIdx) + depth_) / waveLength;
    
    scanlineInterval = cat(1, scanlineInterval, scIntTmp);
    
    figure(1314);
    plot(depth_*1e3, scIntTmp,'LineWidth',2, 'color',[0.2 0.2 0.2]*rIdx); hold on;
end
figure(1314);
line([5 5], [0 5],'LineStyle','--','color','r','LineWidth',2);
hold off;
ylim([0 ceil(max(scanlineInterval(:)))]);
xlabel('Imaging depth (mm)'); ylabel('Radial scanline Interval (\lambda)'); 
legend('r = 5mm', 'r = 10mm', 'r = 15mm','location','northwest');