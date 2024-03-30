function [mOutput,axis_y, axis_z] = load_bfedData(dir_)
load([dir_ '/Parameters.mat']);
dir_file = [dir_ '/errors_bf/error_0/Sample001/Element_64/'];
file_ = dir(dir_file);
load([dir_file '/' file_(end).name]);

acoustic_ = stParam.stRFInfo;
bf_ = stParam.stBFInfo;
trans_ = stParam.stTRInfo;


% mid processing parameter
mid_.nTGC_Atten = 0.5;                % [dB]

mid_.nDCRType = 'high';
mid_.nDCRTap = 128;                   % BPF tap #
mid_.nDCRFcut = 1e6;

% dsc parameter
scanline_theta = linspace(-0.5*bf_.nFOV, 0.5*bf_.nFOV, bf_.nScline); % Ground truth transmitted angle
depth_ = linspace(bf_.nRadius, bf_.nRadius+bf_.nDth, bf_.nDthSpl);

da = abs(scanline_theta(1)-scanline_theta(2));
dr = abs(depth_(1)-depth_(2));
view_depth = bf_.nDth + (bf_.nRadius*(1-cosd(0.5*bf_.nFOV)));
view_width = 2* (bf_.nRadius+bf_.nDth)*sind(0.5*bf_.nFOV);

dz = 1e-4;
dy = 1e-4;
height = round(view_depth / dz);
width = round(view_width / dy);

%%
beamformed_data = stSaveInfo.mBFedData;
% remove NaN
for cols = 1:size(beamformed_data,2)
    col_tmp = beamformed_data(:,cols);
    nans = find(isnan(col_tmp));
    interped = interpn(0:1:size(beamformed_data,1)-1, col_tmp, nans,'spline');
    beamformed_data(nans,cols) = interped;
end
env_data = mid_proc(beamformed_data, mid_, acoustic_, bf_);

[axis_y, axis_z, dsc_data] = ScanConverter_convex(env_data, dr, da, bf_.nRadius, height, width, dz, dy);

aROI = find(dsc_data ~= 50);
aOutlier = find(dsc_data == 50);
mOutput = zeros(size(dsc_data));
mOutput(aROI) = dsc_data(aROI);
mOutput_db = db(mOutput/max(mOutput(:)));
mOutput_db(aOutlier) = -30;

end

