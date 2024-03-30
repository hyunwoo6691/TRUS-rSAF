sigma_ = {'\sigma = 0.1\circ', '\sigma = 0.2\circ', '\sigma = 0.5\circ', ...
    '\sigma = 1\circ', '\sigma = 2\circ', '\sigma = 5\circ', };

figure(1);
subplot(2,7,1);
imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, img_con_gt); hold on;
text(-17, 5, 'GT', 'color', 'w', 'FontWeight', 'bold', 'FontSize', 20); hold off;
colormap gray; caxis([-50 0]);
axis tight; axis equal; xlim([-20 20]);
set(gca, 'FontWeight', 'bold', 'FontSize', 12);
subplot(2,7,8);
imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, img_rsaf_gt); hold on;
text(-17, 5, 'GT', 'color', 'w', 'FontWeight', 'bold', 'FontSize', 20); hold off;
colormap gray; caxis([-50 0]);
axis tight; axis equal; xlim([-20 20]);
set(gca, 'FontWeight', 'bold', 'FontSize', 12);

for e_idx = 1:numel(err_list)
    con_tmp = CON{e_idx};
    con_tmp = con_tmp(:,:,1);
    subplot(2,7,e_idx+1);
    imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, con_tmp); hold on;
    text(-17, 5, sigma_{e_idx}, 'color', 'w', 'FontWeight', 'bold', 'FontSize', 20); hold off;
    colormap gray; caxis([-50 0]);
    axis tight; axis equal; xlim([-20 20]);
    set(gca, 'FontWeight', 'bold', 'FontSize', 12);
    
    rsaf_tmp = rSAF{e_idx};
    rsaf_tmp = rsaf_tmp(:,:,1);
    subplot(2,7,e_idx+8);
    imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, rsaf_tmp); hold on;
    text(-17, 5, sigma_{e_idx}, 'color', 'w', 'FontWeight', 'bold', 'FontSize', 20); hold off;
    colormap gray; caxis([-50 0]);
    axis tight; axis equal; xlim([-20 20]);
    set(gca, 'FontWeight', 'bold', 'FontSize', 12);
end
set(gcf, 'Position', [-1599 572 1600 709]);

%% single frame
output_video = VideoWriter(fullfile(dir_save, '/video_averaged'), 'MPEG-4');
output_video.FrameRate = 30;
open(output_video)

for s_idx = 1:numel(sample_list)
    disp(['>>> sample: ' num2str(s_idx)]);
    figure(1);
    subplot(2,7,1);
    imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, img_con_gt); hold on;
    text(-17, 5, 'GT', 'color', 'w', 'FontWeight', 'bold', 'FontSize', 20); hold off;
    colormap gray; caxis([-50 0]);
    axis tight; axis equal; xlim([-20 20]);
    set(gca, 'FontWeight', 'bold', 'FontSize', 12);
    subplot(2,7,8);
    imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, img_rsaf_gt); hold on;
    text(-17, 5, 'GT', 'color', 'w', 'FontWeight', 'bold', 'FontSize', 20); hold off;
    colormap gray; caxis([-50 0]);
    axis tight; axis equal; xlim([-20 20]);
    set(gca, 'FontWeight', 'bold', 'FontSize', 12);
    
    for e_idx = 1:numel(err_list)
        con_tmp = CON{e_idx};
        con_tmp = con_tmp(:,:,s_idx);
        subplot(2,7,e_idx+1);
        imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, con_tmp); hold on;
        text(-17, 5, sigma_{e_idx}, 'color', 'w', 'FontWeight', 'bold', 'FontSize', 20); hold off;
        colormap gray; caxis([-50 0]);
        axis tight; axis equal; xlim([-20 20]);
        set(gca, 'FontWeight', 'bold', 'FontSize', 12);
        
        rsaf_tmp = rSAF{e_idx};
        rsaf_tmp = rsaf_tmp(:,:,s_idx);
        subplot(2,7,e_idx+8);
        imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, rsaf_tmp); hold on;
        text(-17, 5, sigma_{e_idx}, 'color', 'w', 'FontWeight', 'bold', 'FontSize', 20); hold off;
        colormap gray; caxis([-50 0]);
        axis tight; axis equal; xlim([-20 20]);
        set(gca, 'FontWeight', 'bold', 'FontSize', 12);
    end
    set(gcf, 'Position', [-1599 572 1600 709]);
    f = getframe(gcf);
    
    writeVideo(output_video, f);
end

close(output_video);

%% moving average
average_num = 50;

output_video = VideoWriter(fullfile(dir_save, ['/video_averaged_' num2str(average_num)]), 'MPEG-4');
output_video.FrameRate = 30;
open(output_video)

for s_idx = 1:(numel(sample_list)-average_num+1)
    disp(['>>> sample: ' num2str(s_idx)]);
    figure(1);
    subplot(2,7,1);
    imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, img_con_gt); hold on;
    text(-17, 5, 'GT', 'color', 'w', 'FontWeight', 'bold', 'FontSize', 20); hold off;
    colormap gray; caxis([-50 0]);
    axis tight; axis equal; xlim([-20 20]);
    set(gca, 'FontWeight', 'bold', 'FontSize', 12);
    subplot(2,7,8);
    imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, img_rsaf_gt); hold on;
    text(-17, 5, 'GT', 'color', 'w', 'FontWeight', 'bold', 'FontSize', 20); hold off;
    colormap gray; caxis([-50 0]);
    axis tight; axis equal; xlim([-20 20]);
    set(gca, 'FontWeight', 'bold', 'FontSize', 12);
    
    for e_idx = 1:numel(err_list)
        con_tmp = CON{e_idx};
        con_tmp = mean(con_tmp(:,:,s_idx:s_idx+average_num-1),3);
        subplot(2,7,e_idx+1);
        imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, con_tmp); hold on;
        text(-17, 5, sigma_{e_idx}, 'color', 'w', 'FontWeight', 'bold', 'FontSize', 20); hold off;
        colormap gray; caxis([-50 0]);
        axis tight; axis equal; xlim([-20 20]);
        set(gca, 'FontWeight', 'bold', 'FontSize', 12);
        
        rsaf_tmp = rSAF{e_idx};
        rsaf_tmp = mean(rsaf_tmp(:,:,s_idx:s_idx+average_num-1),3);
        subplot(2,7,e_idx+8);
        imagesc(aYAxis_DSC*1e3, (aZAxis_DSC-stBFInfo.nRadius)*1e3, rsaf_tmp); hold on;
        text(-17, 5, sigma_{e_idx}, 'color', 'w', 'FontWeight', 'bold', 'FontSize', 20); hold off;
        colormap gray; caxis([-50 0]);
        axis tight; axis equal; xlim([-20 20]);
        set(gca, 'FontWeight', 'bold', 'FontSize', 12);
    end
    set(gcf, 'Position', [-1599 572 1600 709]);
    f = getframe(gcf);
    
    writeVideo(output_video, f);
end

close(output_video);