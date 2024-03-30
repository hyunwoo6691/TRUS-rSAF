clc; clear all; close all;

%%
dir_ = uigetdir('./data','');
case_list = dir(dir_);
flag = 0;
if(strcmp(case_list(3).name,'.DS_Store')), flag = 1; end
case_list = case_list(3+flag:end);

%%
numSyn = 45:4:153;
%%
synReg = {};
for c_idx = 1:numel(case_list)
    case_tmp = case_list(c_idx).name;
    disp(case_tmp);
    dir_tmp = [dir_ '/' case_tmp '/errors_bf/error_0/Sample001/Element_64/'];
    file_ = dir(dir_tmp);
    load([dir_tmp file_(end).name]);
    
    synReg_tmp = stSaveInfo.vSynReg;
    synReg = cat(3,synReg, synReg_tmp);
end
%%
fc = 6.5e6;
fs = 4*fc;
c = 1540;
unit_dis = c/fs;
axis_ = linspace(0,2432,2432)*unit_dis/2;

doi = [10 20 30 40 50 60 70]*1e-3; % depth of interest

doi_idx = zeros(1,numel(doi));
for d_idx = 1:numel(doi)
    doi_tmp = doi(d_idx);
    doi_idx(d_idx) = find(abs(axis_-doi_tmp) == min(abs(axis_-doi_tmp)));
end

%%
numSynMap = zeros(numel(case_list),numel(doi));
for c_idx = 1:numel(case_list)
    syn_tmp = numSyn(c_idx);
    synReg_tmp = synReg{c_idx};
    
    scline_tmp = round(0.5*size(synReg_tmp,2));
    
    beam_idx = max(scline_tmp-0.5*(syn_tmp-1),1):min(scline_tmp+0.5*(syn_tmp-1),size(synReg_tmp,3));
    numSyn_scline = sum(squeeze(synReg_tmp(:,scline_tmp,beam_idx)),2);
    
    for d_idx = 1:numel(doi)
        numSynMap(c_idx,d_idx) = numSyn_scline(doi_idx(d_idx));
    end
end
%%
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
for c_idx = 1
    syn_tmp = numSyn(c_idx);
    synReg_tmp = synReg{c_idx};
    
    scline_tmp = round(0.5*size(synReg_tmp,2));
    
    beam_idx = max(scline_tmp-0.5*(syn_tmp-1),1):min(scline_tmp+0.5*(syn_tmp-1),size(synReg_tmp,3));
    mOutput_sum = zeros(height,width);
    for b_idx = 1:numel(beam_idx)
        tmp = squeeze(synReg_tmp(:,:,beam_idx(b_idx)));
        
        [axis_y, axis_z, dsc_data] = ScanConverter_convex(tmp, dr, da, bf_.nRadius, height, width, dz, dy);
        aROI = find(dsc_data ~= 50);
        aOutlier = find(dsc_data == 50);
        mOutput = zeros(size(dsc_data));
        mOutput(aROI) = dsc_data(aROI);
        
        mOutput_sum = mOutput_sum + mOutput;
        
        figure(b_idx);
        imagesc(axis_y, axis_z-bf_.nRadius,mOutput); axis tight; axis equal;
        set(gcf,'Position',[126 718 783 598]);
        pause(0.5);
%         waitforbuttonpress;
    end
    
end


%%
