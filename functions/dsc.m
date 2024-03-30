function [axis_y, axis_z, mOutput, aOutlier] = dsc(env_data, dr, da, bf_, height, width, dz, dy)

[axis_y, axis_z, dsc_data] = ScanConverter_convex(env_data, dr, da, bf_.nRadius, height, width, dz, dy);

aROI = find(dsc_data ~= 50);
aOutlier = find(dsc_data == 50);
mOutput = ones(size(dsc_data))*1e-25;
mOutput(aROI) = dsc_data(aROI);

end

