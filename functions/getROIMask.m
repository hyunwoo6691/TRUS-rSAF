function mask_map = getROIMask(y_axis, z_axis, depth, radius, dor, reject_idx,bShow)

[z_grid, y_grid] = ndgrid(z_axis, y_axis);

dist_map = sqrt(z_grid.^2 + y_grid.^2);

mask_map = logical(zeros(size(dist_map)));

mask_map(dist_map>depth-dor & dist_map<depth+dor) = true;

mask_map(reject_idx) = false;
if(bShow)
    figure(3131);imagesc(y_axis, z_axis-radius,mask_map);
    title('mask map');
    axis tight; axis equal;
end
end

