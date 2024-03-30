function nLatWidth = measure_spatial_resolution(input_img,axis_y, axis_z, dB, fig_ , num_pt , p_idx)

    figure(fig_);
    subplot(num_pt, 1, p_idx);
    imagesc(axis_y,axis_z,input_img); colormap gray; hold on; caxis([-80 0]); axis tight; axis equal;
    [ContourLine, h] = contour(axis_y,axis_z,input_img, dB,'ShowText','off', 'LineColor','r', 'LineWidth', 1);
    
    % detect circle
    nCurrentVal = ContourLine(2,1);
    nCurrentIdx = 1;
    nColumn = size(ContourLine,2);
    
    k = 1;
    mCircleIdx = zeros(2,nColumn);
    mCircleIdx(1,k) = nCurrentIdx;
    mCircleIdx(2,k) = nCurrentVal;
    
    if(nCurrentVal + nCurrentIdx == nColumn)
        nNoVertices = ContourLine(2,1);
        aX_Contour = ContourLine(1,2:(2+nNoVertices-1));
        aY_Contour = ContourLine(2,2:(2+nNoVertices-1));
        nLatWidth = abs(max(aX_Contour) - min(aX_Contour));
        nAxlWidth = abs(max(aY_Contour)-min(aY_Contour));
    else
        nPrevVal = nCurrentVal;
        nPreIdx = nCurrentIdx;
        while(nCurrentIdx+nCurrentVal < nColumn)
            k = k+1;
            nNextIdx = nCurrentIdx +1 + nCurrentVal;
            nNext = ContourLine(2, nNextIdx);
            nPreVal = nCurrentVal;
            nPreIdx = nCurrentIdx;
            nCurrentIdx = nNextIdx;
            nCurrentVal = nNext;
            
            mCircleIdx(1,k) = nCurrentIdx;
            mCircleIdx(2,k) = nCurrentVal;
        end
        nMaxVal = max(mCircleIdx(2,:));
        nMaxIdx = mCircleIdx(1,find(mCircleIdx(2,:)==nMaxVal));
        
        aX_Contour = ContourLine(1, (nMaxIdx+1):(nMaxIdx+nMaxVal));
        aY_Contour = ContourLine(2, (nMaxIdx+1):(nMaxIdx+nMaxVal));
        nLatWidth = abs(max(aX_Contour)-min(aX_Contour));
    end

end

