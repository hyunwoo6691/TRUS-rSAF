depth = 30.3e-3;

scanline = 64;

no_syn = 16;

beams_used = max(scanline-0.5*no_syn, 1):min(scanline+0.5*no_syn-1, 128);

v_beamformed = vBFedData(:, beams_used, :);

depth_p = depth+stBFInfo.nRadius;
roi_ = depth_p * cosd(aSclineTheta);
roi_ = repmat(roi_,numel(aDth),1);
find_ = abs(roi_ - mImgZ);

idx_ = zeros(1, 128);
for k = 1:128
    tmp = find(find_(:,k) == min(find_(:,k)));
%     idx_(k) = tmp(1);
    idx_(k) = 2045;
end

m_psf = zeros(no_syn, 128);
for k = 1:no_syn
    beamformed_tmp = squeeze(v_beamformed(:,k,:));
    for e = 1:128
        m_psf(k, e) = beamformed_tmp(idx_(e), e);
    end
end

psf_sum = sum(m_psf)/sqrt(no_syn);

%% plot
close all;
figure(1); hold on;
for k = 1:no_syn
    subplot(2, no_syn/2, k); 
    plot(m_psf(k,:), 'LineWidth', 2); grid minor;
    line([64, 64], [min(m_psf(:)), max(m_psf(:))],'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
%     line([64, 64], [0, 1],'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    title(['Beam idx: ', num2str(beams_used(k))]);
    ylim([min(m_psf(:)) max(m_psf(:))]);
end
hold off;
set(gcf, 'Position', [957 486 1630 720]);

figure(2); hold on;
for k = 1:no_syn
    plot(m_psf(k,:)/max(psf_sum), 'LineWidth', 1);
end
% plot(psf_sum/max(psf_sum), 'LineWidth', 2.5, 'LineStyle', '-.', 'Color', 'k');
plot(mBFedData(idx_(1),:)/max(mBFedData(idx_(1),:)), 'LineWidth', 2.5, 'LineStyle', '-.', 'Color', 'r');
grid minor;

%% for conventional
depth = 30.49e-3;

depth_p = depth+stBFInfo.nRadius;

scanline = 64;

aZ = mImgZ(:,scanline);
idx_ = find(abs(aZ - depth_p) == min(abs(aZ - depth_p)));
idx_ = idx_(1);

psf = mBFedData(idx_, :);
% psf(1:127) = psf(2:end);
% psf(end) = 0;

figure(2); 
plot(psf/max(psf), 'LineWidth', 2.5, 'LineStyle', '-.', 'Color','r'); 
grid minor;
% ylim([-0.2 1]);

%%
psf_db = mBFedData(idx_(1),:)/max(mBFedData(idx_(1),:));
figure(3);hold on;
plot(db(psf_db), 'LineWidth', 2.5, 'LineStyle', '-.');
ylim([-80 0]);
set(gcf, 'Position', [1 45 1680 616]);

%%
figure(3);
line([40 88], [-6 -6], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r');
text(80, -7, '-6dB', 'FontSize', 14);

grid minor;


