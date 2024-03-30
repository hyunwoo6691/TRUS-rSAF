function [dsc_out,axis_x, axis_y, axis_z] = dsc_volume(input,volume_size, voxel_size, orientation, volume_grid_x, volume_grid_z)
tic;
volume_sample = round(volume_size./voxel_size);
num_sample = volume_sample(1)*volume_sample(2)*volume_sample(3);

axis_x = linspace(-0.5*volume_size(1), 0.5*volume_size(1), volume_sample(1));
axis_y = linspace(-0.5*volume_size(2), 0.5*volume_size(2), volume_sample(2));
axis_z = linspace(0, volume_size(3), volume_sample(3));

[dsc_grid_x, dsc_grid_y, dsc_grid_z] = ndgrid(axis_x, axis_y, axis_z);

dsc_out = zeros(numel(axis_x), numel(axis_y), numel(axis_z));

d_theta = abs(orientation(1)-orientation(2));

for v_idx = 1:num_sample
    
    if(mod(v_idx,round(num_sample/5))==0), disp(['>>> DSC: ' num2str(v_idx) '/' num2str(num_sample)]); end
    
    x_ = dsc_grid_x(v_idx);
    y_ = dsc_grid_y(v_idx);
    z_ = dsc_grid_z(v_idx);
    
    % find closest coordinate
    % orientation
    theta_ = atand(y_/z_);
    if(theta_ < orientation(1) || theta_>orientation(end)), continue; end
    th_idx = (theta_ - orientation(1))/d_theta + 1;
    th_idx_int = floor(th_idx);
    th_idx_fr = th_idx - th_idx_int;
    
    % x
    x_axis_tmp = volume_grid_x(1,:,th_idx_int);
    d_x = abs(x_axis_tmp(1) - x_axis_tmp(2));
    x_idx = (x_ - x_axis_tmp(1))/d_x + 1;
    if(x_idx < 1 || x_idx > numel(x_axis_tmp)), continue; end
    x_idx_int = floor(x_idx);
    x_idx_fr = x_idx - x_idx_int;
    
    % z
    z_axis_tmp = volume_grid_z(:,1,th_idx_int);
    d_z = abs(z_axis_tmp(1) - z_axis_tmp(2));
    z_idx = (z_ - z_axis_tmp(1))/d_z + 1;
    if(z_idx< 1 || z_idx > numel(z_axis_tmp)), continue; end
    z_idx_int = floor(z_idx);
    z_idx_fr = z_idx - z_idx_int;
    
    % only when input sample exists
    if((th_idx_int < numel(orientation)) && (x_idx_int < numel(x_axis_tmp)) && (z_idx_int < numel(z_axis_tmp)))
        % get the vertices of the local cube
        c000 = input(z_idx_int, x_idx_int, th_idx_int);
        c001 = input(z_idx_int, x_idx_int, th_idx_int+1);
        c010 = input(z_idx_int, x_idx_int+1, th_idx_int);
        c011 = input(z_idx_int, x_idx_int+1, th_idx_int+1);
        c100 = input(z_idx_int+1, x_idx_int, th_idx_int);
        c101 = input(z_idx_int+1, x_idx_int, th_idx_int+1);
        c110 = input(z_idx_int+1, x_idx_int+1, th_idx_int);
        c111 = input(z_idx_int+1, x_idx_int+1, th_idx_int+1);
        
        % interpolate along x-axis
        c00 = c000*(1-x_idx_fr) + c010*x_idx_fr;
        c01 = c100*(1-x_idx_fr) + c110*x_idx_fr;
        c10 = c001*(1-x_idx_fr) + c011*x_idx_fr;
        c11 = c101*(1-x_idx_fr) + c111*x_idx_fr;
        
        % interpolate along z-axis
        c0 = c00*(1-z_idx_fr) + c01*z_idx_fr;
        c1 = c10*(1-z_idx_fr) + c11*z_idx_fr;
        
        % interpolate along y-axis (between slices)
        dsc_out(v_idx) = c0*(1-th_idx_fr) * c1*th_idx_fr;
    end
end
toc;
end

