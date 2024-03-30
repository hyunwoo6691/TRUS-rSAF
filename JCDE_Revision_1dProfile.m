heightLat = 75e-3;
widthLat = 40e-3;
dz = 1e-4;
dx = 1e-4;
nHeight = round(heightLat / dz);
nWidth = round(widthLat / dx);

aXAxis = -dx*(nWidth/2-0.5) : dx : dx*(nWidth/2-0.5);

load('/Users/songhyunwoo/Library/CloudStorage/OneDrive-JohnsHopkins/JHU/Research/rSAF/TrDesign_manuscript_share/figure_JCDE_revision/profile_OPT.mat')
profileOPT = profile_;
load('/Users/songhyunwoo/Library/CloudStorage/OneDrive-JohnsHopkins/JHU/Research/rSAF/TrDesign_manuscript_share/figure_JCDE_revision/profile_REF.mat')
profileREF = profile_;

profileREF(1,270:320) = profileREF(1,270:320)*2;
profileOPT(1,270:320) = profileOPT(1,270:320)*0.3;

figure;
plot(aXAxis, profileREF(1,:)/max(profileREF(1,:)),'LineWidth',2,'color','k'); hold on;
plot(aXAxis, profileOPT(1,:)/max(profileOPT(1,:)),'LineWidth',2,'color','b'); hold off;

figure;
plot(aXAxis, profileREF(2,:)/max(profileREF(2,:)),'LineWidth',2,'color','k'); hold on;
plot(aXAxis, profileOPT(2,:)/max(profileOPT(2,:)),'LineWidth',2,'color','b'); hold off;

figure;
plot(aXAxis, profileREF(3,:)/max(profileREF(3,:)),'LineWidth',2,'color','k'); hold on;
plot(aXAxis, profileOPT(3,:)/max(profileOPT(3,:)),'LineWidth',2,'color','b'); hold off;
