function output = InterpNan(input)
[sample_, channel] = size(input);

output = zeros(size(input));
for c_idx = 1:channel
    aline_tmp = input(:,c_idx);
    output(:,c_idx) = interpn(0:1:sample_-1, aline_tmp, 0:1:sample_-1,'spline')';
    clc;
end

end

