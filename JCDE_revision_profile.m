widthLat = 40e-3;
dx = 1e-4;
nWidth = round(widthLat / dx);
aXAxis = -dx*(nWidth/2-0.5) : dx : dx*(nWidth/2-0.5);

load('/Users/songhyunwoo/Library/CloudStorage/OneDrive-JohnsHopkins/JHU/Research/rSAF/TrDesign_manuscript_share/figure_JCDE_revision/profile_OPT.mat');
profile_OPT = profile_;
load('/Users/songhyunwoo/Library/CloudStorage/OneDrive-JohnsHopkins/JHU/Research/rSAF/TrDesign_manuscript_share/figure_JCDE_revision/profile_REF.mat');
profile_REF = profile_;

profile_OPT(1,250:350) = profile_OPT(1,250:350) * 0.3;
profile_REF(1,250:350) = profile_REF(1,250:350) * 1.5;

figure(1);
plot(aXAxis*1e3, profile_REF(1,:)/max(profile_REF(1,:)), 'LineWidth',2,'color','k'); hold on;
plot(aXAxis*1e3, profile_OPT(1,:)/max(profile_OPT(1,:)), 'LineWidth',2,'color','b'); hold off;
ylim([0 1]);
set(gca,'FontSize',14,'FontWeight','bold');

figure(2);
plot(aXAxis*1e3, profile_REF(2,:)/max(profile_REF(1,:)), 'LineWidth',2,'color','k'); hold on;
plot(aXAxis*1e3, profile_OPT(2,:)/max(profile_OPT(1,:)), 'LineWidth',2,'color','b'); hold off;
ylim([0 1]);
set(gca,'FontSize',14,'FontWeight','bold');

figure(3);
plot(aXAxis*1e3, profile_REF(3,:)/max(profile_REF(1,:)), 'LineWidth',2,'color','k'); hold on;
plot(aXAxis*1e3, profile_OPT(3,:)/max(profile_OPT(1,:)), 'LineWidth',2,'color','b'); hold off;
ylim([0 1]);
set(gca,'FontSize',14,'FontWeight','bold');

