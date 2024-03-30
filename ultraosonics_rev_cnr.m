clc; clear; close all;

%%
addpath( '/Users/songhyunwoo/Desktop');

rsaf_gt = double(im2gray(imread('rsaf_gt.png')));
rsaf0_1 = double(im2gray(imread('rsaf0_1.png')));
rsaf0_2 = double(im2gray(imread('rsaf0_2.png')));
rsaf0_5 = double(im2gray(imread('rsaf0_5.png')));
rsaf1 = double(im2gray(imread('rsaf1.png')));
rsaf2 = double(im2gray(imread('rsaf2.png')));
rsaf5 = double(im2gray(imread('rsaf5.png')));

con_gt = double(im2gray(imread('con_gt.png')));
con0_1 = double(im2gray(imread('con0_1.png')));
con0_2 = double(im2gray(imread('con0_2.png')));
con0_5 = double(im2gray(imread('con0_5.png')));
con1 = double(im2gray(imread('con1.png')));
con2 = double(im2gray(imread('con2.png')));
con5 = double(im2gray(imread('con5.png')));

%%
cnr_rSAF = zeros(1,7);
cnr_CON = zeros(1,7);
%% rSAF
figure(1);
imagesc(rsaf_gt);
roi_tmp = drawrectangle();
roiIdx = round(roi_tmp.Position);

roi = rsaf_gt(roiIdx(2):roiIdx(2)+roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
bg = rsaf_gt(roiIdx(2):roiIdx(2)+1.5*roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
meanROI = mean(roi(:)); meanBG = mean(bg(:));
stdROI = std(roi(:)); stdBG = std(bg(:));
cnr_rSAF(1) = abs(meanROI-meanBG)/sqrt(stdROI^2+stdBG^2);

roi = rsaf0_1(roiIdx(2):roiIdx(2)+roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
bg = rsaf0_1(roiIdx(2):roiIdx(2)+1.5*roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
meanROI = mean(roi(:)); meanBG = mean(bg(:));
stdROI = std(roi(:)); stdBG = std(bg(:));
cnr_rSAF(2) = abs(meanROI-meanBG)/sqrt(stdROI^2+stdBG^2);

roi = rsaf0_2(roiIdx(2):roiIdx(2)+roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
bg = rsaf0_2(roiIdx(2):roiIdx(2)+1.5*roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
meanROI = mean(roi(:)); meanBG = mean(bg(:));
stdROI = std(roi(:)); stdBG = std(bg(:));
cnr_rSAF(3) = abs(meanROI-meanBG)/sqrt(stdROI^2+stdBG^2);

roi = rsaf0_5(roiIdx(2):roiIdx(2)+roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
bg = rsaf0_5(roiIdx(2):roiIdx(2)+1.5*roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
meanROI = mean(roi(:)); meanBG = mean(bg(:));
stdROI = std(roi(:)); stdBG = std(bg(:));
cnr_rSAF(4) = abs(meanROI-meanBG)/sqrt(stdROI^2+stdBG^2);

roi = rsaf1(roiIdx(2):roiIdx(2)+roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
bg = rsaf1(roiIdx(2):roiIdx(2)+1.5*roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
meanROI = mean(roi(:)); meanBG = mean(bg(:));
stdROI = std(roi(:)); stdBG = std(bg(:));
cnr_rSAF(5) = abs(meanROI-meanBG)/sqrt(stdROI^2+stdBG^2);

roi = rsaf2(roiIdx(2):roiIdx(2)+roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
bg = rsaf2(roiIdx(2):roiIdx(2)+1.5*roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
meanROI = mean(roi(:)); meanBG = mean(bg(:));
stdROI = std(roi(:)); stdBG = std(bg(:));
cnr_rSAF(6) = abs(meanROI-meanBG)/sqrt(stdROI^2+stdBG^2);

roi = rsaf5(roiIdx(2):roiIdx(2)+roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
bg = rsaf5(roiIdx(2):roiIdx(2)+1.5*roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
meanROI = mean(roi(:)); meanBG = mean(bg(:));
stdROI = std(roi(:)); stdBG = std(bg(:));
cnr_rSAF(7) = abs(meanROI-meanBG)/sqrt(stdROI^2+stdBG^2);

%% CON
figure(2);
imagesc(con_gt);
roi_tmp = drawrectangle();
roiIdx = round(roi_tmp.Position);

roi = con_gt(roiIdx(2):roiIdx(2)+roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
bg = con_gt(roiIdx(2):roiIdx(2)+1.5*roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
meanROI = mean(roi(:)); meanBG = mean(bg(:));
stdROI = std(roi(:)); stdBG = std(bg(:));
cnr_CON(1) = abs(meanROI-meanBG)/sqrt(stdROI^2+stdBG^2);

roi = con0_1(roiIdx(2):roiIdx(2)+roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
bg = con0_1(roiIdx(2):roiIdx(2)+1.5*roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
meanROI = mean(roi(:)); meanBG = mean(bg(:));
stdROI = std(roi(:)); stdBG = std(bg(:));
cnr_CON(2) = abs(meanROI-meanBG)/sqrt(stdROI^2+stdBG^2);

roi = con0_2(roiIdx(2):roiIdx(2)+roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
bg = con0_2(roiIdx(2):roiIdx(2)+1.5*roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
meanROI = mean(roi(:)); meanBG = mean(bg(:));
stdROI = std(roi(:)); stdBG = std(bg(:));
cnr_CON(3) = abs(meanROI-meanBG)/sqrt(stdROI^2+stdBG^2);

roi = con0_5(roiIdx(2):roiIdx(2)+roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
bg = con0_5(roiIdx(2):roiIdx(2)+1.5*roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
meanROI = mean(roi(:)); meanBG = mean(bg(:));
stdROI = std(roi(:)); stdBG = std(bg(:));
cnr_CON(4) = abs(meanROI-meanBG)/sqrt(stdROI^2+stdBG^2);

roi = con1(roiIdx(2):roiIdx(2)+roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
bg = con1(roiIdx(2):roiIdx(2)+1.5*roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
meanROI = mean(roi(:)); meanBG = mean(bg(:));
stdROI = std(roi(:)); stdBG = std(bg(:));
cnr_CON(5) = abs(meanROI-meanBG)/sqrt(stdROI^2+stdBG^2);

roi = con2(roiIdx(2):roiIdx(2)+roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
bg = con2(roiIdx(2):roiIdx(2)+1.5*roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
meanROI = mean(roi(:)); meanBG = mean(bg(:));
stdROI = std(roi(:)); stdBG = std(bg(:));
cnr_CON(6) = abs(meanROI-meanBG)/sqrt(stdROI^2+stdBG^2);

roi = con5(roiIdx(2):roiIdx(2)+roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
bg = con5(roiIdx(2):roiIdx(2)+1.5*roiIdx(4), roiIdx(1):roiIdx(1)+roiIdx(3));
meanROI = mean(roi(:)); meanBG = mean(bg(:));
stdROI = std(roi(:)); stdBG = std(bg(:));
cnr_CON(7) = abs(meanROI-meanBG)/sqrt(stdROI^2+stdBG^2);

%%
axis_ = [0 0.1 0.2 0.5 1 2 5];
figure(533);
plot(axis_,cnr_CON,'LineWidth',4,'color','k'); hold on;
plot(axis_,cnr_rSAF,'LineWidth',4,'color','b'); hold off;
grid on; xlabel('\sigma'); ylabel('CNR');
