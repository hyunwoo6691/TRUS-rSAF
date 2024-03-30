%% Image output
[aYAxis, aZAxis, mOutputTmp] = ScanConverter_convex(syn_16, dr, da, nRadius, stDSCInfo.nHeight, stDSCInfo.nWidth, stDSCInfo.dz, stDSCInfo.dy);

aROI = find(mOutputTmp ~= 50);
aOutlier = find(mOutputTmp == 50);
mOutput = zeros(size(mOutputTmp));
mOutput(aROI) = mOutputTmp(aROI);
mOutput(aOutlier) = 1e-2*min(mOutput(aROI));
mOutput_dB = db(mOutput/max(mOutput(:)));
% mOutput_dB = db(mOutput/max_syn);

%%%%%%%% Normailze at each depth
% mOutput_norm = zeros(size(mOutput));
% for z_idx = 1:numel(aZAxis)
%     mOutput_norm(z_idx,:) = mOutput(z_idx,:)/max(mOutput(z_idx,:));
% end
% mOutput_dB = db(mOutput_norm/max(mOutput_norm(:)));

figure; imagesc(aYAxis*1e3,(aZAxis-nRadius)*1e3, mOutput_dB); hold on;
hold off; colormap jet;
axis tight; axis equal; xlabel('Elevational [mm]'); ylabel('Axial [mm]');
set(gca,'Fontsize',20,'FontName','Times New Roman', 'FontWeight','bold');
set(gcf,'Position',[-2935 104 351 676]);
xlim([-10 10]); %ylim([0 90]);
caxis([-100 0]);

%% get mask for each radial depth
depth_ = (3:7) * 1e-2; % radial direction
offset = 0.1e-3;
% offset = 0.05e-3;
% depth = depth_ + offset;
depth = depth_;
depth(2:5) = depth(2:5) + 0.1e-3;
% depth = depth_;

y_axis = aYAxis_DSC;
z_axis = aZAxis_DSC - stBFInfo.nRadius;

[z_grid, y_grid] = ndgrid(z_axis, y_axis);

dist_map = sqrt(z_grid.^2 + y_grid.^2);

mask_map = logical(zeros(size(dist_map)));
for d_idx = 1:numel(depth)
    for k = 1:numel(y_axis)
        tmp = dist_map(:,k);
        idx = find(abs(tmp - depth(d_idx)) == min(abs(tmp-depth(d_idx))));
        mask_map(idx(1), k) = true;
    end
end

reject_area = (img_con_gt == -30);
mask_map(reject_area) = false;
figure(3);imagesc(mask_map); title('mask map');
%% take -6dB width
num_syn = [4 8 16 32];
depths = [3 4 5 6 7 8 9]*10e-3;
depths = 50e-3;
width_ = zeros(numel(depths), numel(num_syn));

for d_idx = 1:numel(depths)
    depth_ = depths(d_idx);
    for s_idx = 1:numel(num_syn)
        syn_ = num_syn(s_idx);
        beam_pattern_tmp = eval(['syn_' num2str(syn_)]);
        [aYAxis, aZAxis, mOutputTmp] = ScanConverter_convex(beam_pattern_tmp, dr, da, nRadius, stDSCInfo.nHeight, stDSCInfo.nWidth, stDSCInfo.dz, stDSCInfo.dy);
        
        aROI = find(mOutputTmp ~= 50);
        aOutlier = find(mOutputTmp == 50);
        mOutput = zeros(size(mOutputTmp));
        mOutput(aROI) = mOutputTmp(aROI);
        mOutput(aOutlier) = 1e-2*min(mOutput(aROI));
        
        mOutput_norm = zeros(size(mOutput));
        for z_idx = 1:numel(aZAxis)
            mOutput_norm(z_idx,:) = mOutput(z_idx,:)/max(mOutput(z_idx,:));
        end
        mOutput_dB = db(mOutput_norm/max(mOutput_norm(:)));
        
        % mask
        y_axis = aYAxis;
        z_axis = aZAxis - nRadius;
        
        [z_grid, y_grid] = ndgrid(z_axis, y_axis);
        
        dist_map = sqrt(z_grid.^2 + y_grid.^2);
        
        mask_map = logical(zeros(size(dist_map)));
        for k = 1:numel(y_axis)
            tmp = dist_map(:,k);
            idx = find(abs(tmp - depth_) == min(abs(tmp-depth_)));
            mask_map(idx(1), k) = true;
        end
        
        %     depth_idx = find(abs(aZAxis-nRadius-depth_)==min(abs(aZAxis-nRadius-depth_)));
        
        %     psf_ = mOutput_dB(depth_idx,:); % elevational
        psf_ = mOutput_dB(mask_map)'; % radial
        y_interpn = linspace(aYAxis(1), aYAxis(end), 10*numel(aYAxis));
        psf_ = interpn(aYAxis, psf_, y_interpn, 'linear');
        
        peak_idx = find(psf_==max(psf_),1,'first');
        psf_left = psf_(1:peak_idx-1); y_left = y_interpn(1:peak_idx-1);
        psf_right = psf_(peak_idx+1:end); y_right = y_interpn(peak_idx+1:end);
        
        half_left_idx = find(psf_left>=-6, 1, 'first'); left_idx = half_left_idx;
        half_right_idx = find(psf_right>=-6, 1, 'last'); right_idx = peak_idx + half_right_idx;
        
        width_(d_idx,s_idx) = y_right(half_right_idx) - y_left(half_left_idx);
        
        figure(d_idx); subplot(4,1, s_idx);
        plot(y_interpn, psf_); hold on;
        scatter(y_interpn(left_idx), psf_(left_idx));
        scatter(y_interpn(right_idx), psf_(right_idx)); hold off;
        ylim([-10 0]);
    end
    figure(d_idx*10);
    plot(num_syn, width_(d_idx,:)*1e3);
    ylabel('FWHM [mm]'); xlabel('N_syn');
end

%% Overlay outermost beams
v_syn_side = v_syn16(:,:,[1 end]);
for f_idx = 1:size(v_syn_side,3)
    beam_pattern_tmp = v_syn_side(:,:,f_idx)';
    [aYAxis, aZAxis, mOutputTmp] = ScanConverter_convex(beam_pattern_tmp, dr, da, nRadius, stDSCInfo.nHeight, stDSCInfo.nWidth, stDSCInfo.dz, stDSCInfo.dy);
    
    aROI = find(mOutputTmp ~= 50);
    aOutlier = find(mOutputTmp == 50);
    mOutput = zeros(size(mOutputTmp));
    mOutput(aROI) = mOutputTmp(aROI);
    mOutput(aOutlier) = 1e-2*min(mOutput(aROI));
    
    mOutput_dB = db(mOutput/max(mOutput(:)));
    
    if(f_idx == 1), A = mOutput_dB; else B = mOutput_dB; end
    
    figure(10131+f_idx);
    imagesc(aYAxis*1e3,(aZAxis-nRadius)*1e3, mOutput_dB);
    axis tight; axis equal;
    set(gcf, 'Position', [1857 541 325 561]);
    xlim([-5 5]); ylim([0 90]);
    hold off; colormap default; caxis([-100 0]);
end

figure(413221);
a1 = axes;
h1 = imagesc(a1, A); colormap(a1, 'default'); c1 = colorbar; axis image; caxis([-100 0]);

a2 = axes;
h2 = imagesc(a2, B); colormap(a2, 'default'); c2 = colorbar; axis image; caxis([-100 0]);

alpha(h2, isfinite(B)*0.5);

a1.Visible = 'off';
a2.Visible = 'off';

c1.Visible = 'off';

set(gcf, 'Position', [1885 -66 942 717]);


%% contour
colors = ['k', 'b'];
v_syn_side = v_syn16(:,:,[1 end]);
depth = 50e-3;
roi_size = [5e-3 10e-3];
% normalize in each depth
for f_idx = 1:size(v_syn_side,3)
    beam_pattern_tmp = v_syn_side(:,:,f_idx)';
    [aYAxis, aZAxis, mOutputTmp] = ScanConverter_convex(beam_pattern_tmp, dr, da, nRadius, stDSCInfo.nHeight, stDSCInfo.nWidth, stDSCInfo.dz, stDSCInfo.dy);
    
    aROI = find(mOutputTmp ~= 50);
    aOutlier = find(mOutputTmp == 50);
    mOutput = zeros(size(mOutputTmp));
    mOutput(aROI) = mOutputTmp(aROI);
    mOutput(aOutlier) = 1e-2*min(mOutput(aROI));
    
    mOutput_norm = zeros(size(mOutput));
    for z_idx = 1:numel(aZAxis)
        mOutput_norm(z_idx,:) = mOutput(z_idx,:)/max(mOutput(z_idx,:));
    end
    mOutput_dB = db(mOutput_norm/max(mOutput_norm(:)));
    
    % roi_position
    up_idx = find(abs(aZAxis-(depth-roi_size(2)))==min(abs(aZAxis-(depth-roi_size(2)))));
    down_idx = find(abs(aZAxis-(depth+roi_size(2)))==min(abs(aZAxis-(depth+roi_size(2)))));
    left_idx = find(abs(aYAxis+roi_size(1))==min(abs(aYAxis+roi_size(1))));
    right_idx = find(abs(aYAxis-roi_size(1))==min(abs(aYAxis-roi_size(1))));
    
    roi_ = mOutput_dB(up_idx:down_idx,left_idx:right_idx);
    roi_ = roi_ + max(roi_(:)); % normalize
    
    if(f_idx == 1), A = roi_; else B = roi_; end
    
    figure(10131+f_idx);
    imagesc(aYAxis(left_idx:right_idx)*1e3,(aZAxis(up_idx:down_idx))*1e3, roi_); hold on;
    [ContourLine, h] = contour(aYAxis(left_idx:right_idx)*1e3,(aZAxis(up_idx:down_idx))*1e3, roi_, ...
                                [-6 -6],'ShowText','on', 'LineColor', colors(f_idx), 'LineWidth', 2.5);
    axis tight; axis equal;
    set(gca, 'Ydir', 'reverse');
    set(gcf, 'Position', [1857 541 325 561]);
%     xlim([-5 5]); ylim([depth-10e-3 depth+10e-3]*1e3);
    hold off; colormap default;
end
% overlay
figure(41341);
a1 = axes;
h1 = imagesc(a1, A); colormap(a1, 'default'); c1 = colorbar; axis image; caxis([-50 0]);

a2 = axes;
h2 = imagesc(a2, B); colormap(a2, 'default'); c2 = colorbar; axis image; caxis([-50 0]);

alpha(h2, isfinite(B)*0.5);

a1.Visible = 'off';
a2.Visible = 'off';

c1.Visible = 'off';

hold on;
[ContourLine, h] = contour(A, [-6 -6],'ShowText','off', 'LineColor', 'r', 'LineWidth', 3);
[ContourLine, h] = contour(B, [-6 -6],'ShowText','off', 'LineColor', 'b', 'LineWidth', 3);

set(gcf, 'Position', [1857 -55 325 561]);

hold off;








