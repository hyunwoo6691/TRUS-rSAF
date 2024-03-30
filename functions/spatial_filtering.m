function [output, fil] = spatial_filtering(input, sigma_, d_theta, shape)

if(~strcmp(shape,'none'))
%     disp(['spatial filtering: ' shape]);
end

% filter length
% L = round((2*sigma_/d_theta)/2)*2 + 1; % make it odd number
% L = round((2*sigma_/d_theta)/2)*2; % make it even number
L = 131; % fix to the length same as N_syn
% L = 5; % fix to the length same as N_syn

switch shape
    case 'gauss'
        alpha = (L-1)/(2*sigma_);
        fil_tmp = gausswin(L, alpha);
        fil = fil_tmp / sum(fil_tmp);
    case 'boxcar'
        fil_tmp = ones(1,L);
        fil = fil_tmp / sum(fil_tmp);
end

if(strcmp(shape, 'none'))
    output = input;
    fil = 'none';
else
    output = zeros(size(input));
    for z_idx = 1:size(input,1)
        tmp = input(z_idx,:);
        output(z_idx,:) = conv(tmp, fil, 'same');
    end
end

end

