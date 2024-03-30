clc; close all; clear;
%%
dir_ = uigetdir('./data', '');

dir_tmp = [dir_ '/r5mm_lens_focus_5mm/errors_bf/error_0/Sample001/Element_64'];

file_list = dir(dir_tmp);
flag = 0;
if(strcmp(file_list(3).name,'.DS_Store')), flag = 1; end
file_list = file_list(3+flag:end);

clear flag
%%
load([dir_ '/r5mm_lens_focus_5mm/Parameters.mat']);

acoustic_ = stParam.stRFInfo;
bf_ = stParam.stBFInfo;
trans_ = stParam.stTRInfo;

%% cluster beamforming type
bf_type = {}; num_syn = [];
for f_idx = 1:numel(file_list)
    file_tmp = split(file_list(f_idx).name,'_');
    file_tmpp = split(file_tmp{1}, ']');
    switch file_tmpp{2}
        case 'Conv'
            bf_type = cat(1, bf_type, 'CON');
        case 'SA'
            bf_type = cat(1, bf_type, 'rSAF');
            
            syn_tmp = split(file_tmp{3},'.');
            syn_tmpp = split(syn_tmp{1}, 'n');
            num_syn = cat(1, num_syn, str2num(syn_tmpp{2}));
    end
end

%%
radius_ = 5;
con_ = ['r12_' num2str(radius_) 'mm_lens_focus_20mm'];
rsaf_foc5 = ['r12_' num2str(radius_) 'mm_lens_focus_5mm'];
rsaf_foc10 = ['r12_' num2str(radius_) 'mm_lens_focus_10mm'];
rsaf_foc15 = ['r12_' num2str(radius_) 'mm_lens_focus_15mm'];
rsaf_foc20 = ['r12_' num2str(radius_) 'mm_lens_focus_20mm'];

dir_con = [dir_ '/' con_ '/sidelobe/sidelobe_CON.mat']; load(dir_con);
dir_rsaf_foc5 = [dir_ '/' rsaf_foc5 '/sidelobe/sidelobe_rSAF.mat']; load(dir_rsaf_foc5);
rsaf_1d_foc5 = rsaf_1d;
dir_rsaf_foc10 = [dir_ '/' rsaf_foc10 '/sidelobe/sidelobe_rSAF.mat']; load(dir_rsaf_foc10);
rsaf_1d_foc10 = rsaf_1d;
dir_rsaf_foc15 = [dir_ '/' rsaf_foc15 '/sidelobe/sidelobe_rSAF.mat']; load(dir_rsaf_foc15);
rsaf_1d_foc15 = rsaf_1d;
dir_rsaf_foc20 = [dir_ '/' rsaf_foc20 '/sidelobe/sidelobe_rSAF.mat']; load(dir_rsaf_foc20);
rsaf_1d_foc20 = rsaf_1d;

clear rsaf_1d

%%
% syns = [41 19 13 5]; % focal point: 5mm, 10mm, 15mm, 20mm       r = 5mm
% syns = [41 17 7 5]; % focal point: 5mm, 10mm, 15mm, 20mm       r = 7.5mm
% syns = [41 15 5 3]; % focal point: 5mm, 10mm, 15mm, 20mm         r = 10mm
syns = [41 15 5 3]; % focal point: 5mm, 10mm, 15mm, 20mm       r = 12.5mm
%%
depths_ = 1:7;

% figure(1);
for d_idx = depths_
%     d_idx = 3; % number indicates depth [cm]
    
    con_1d_ = con_1d(:,d_idx);
    con_1d_ = con_1d_(con_1d_ ~= 50);
    
%     subplot(numel(depths_), 1, d_idx);
    figure(d_idx);
    % for CON
    radial_axis = linspace(-0.5*bf_.nFOV, 0.5*bf_.nFOV, numel(con_1d_));
    plot(radial_axis, con_1d_/max(con_1d_), 'LineWidth', 2); hold on;
    
    % for rSAF
    idx_tmp = syns(1); idx_tmpp = find(idx_tmp == num_syn); rsaf_tmp = rsaf_1d_foc5{idx_tmpp};
    rsaf_tmpp = rsaf_tmp(:,d_idx); rsaf_tmpp = rsaf_tmpp(rsaf_tmpp ~= 50);
    radial_axis = linspace(-0.5*bf_.nFOV, 0.5*bf_.nFOV, numel(rsaf_tmpp));
    plot(radial_axis, rsaf_tmpp/max(rsaf_tmpp), 'LineWidth', 2);
    
    idx_tmp = syns(2); idx_tmpp = find(idx_tmp == num_syn); rsaf_tmp = rsaf_1d_foc10{idx_tmpp};
    rsaf_tmpp = rsaf_tmp(:,d_idx); rsaf_tmpp = rsaf_tmpp(rsaf_tmpp ~= 50);
    radial_axis = linspace(-0.5*bf_.nFOV, 0.5*bf_.nFOV, numel(rsaf_tmpp));
    plot(radial_axis, rsaf_tmpp/max(rsaf_tmpp), 'LineWidth', 2);
    
    idx_tmp = syns(3); idx_tmpp = find(idx_tmp == num_syn); rsaf_tmp = rsaf_1d_foc15{idx_tmpp};
    rsaf_tmpp = rsaf_tmp(:,d_idx); rsaf_tmpp = rsaf_tmpp(rsaf_tmpp ~= 50);
    radial_axis = linspace(-0.5*bf_.nFOV, 0.5*bf_.nFOV, numel(rsaf_tmpp));
    plot(radial_axis, rsaf_tmpp/max(rsaf_tmpp), 'LineWidth', 2);
    
    idx_tmp = syns(4); idx_tmpp = find(idx_tmp == num_syn); rsaf_tmp = rsaf_1d_foc20{idx_tmpp};
    rsaf_tmpp = rsaf_tmp(:,d_idx); rsaf_tmpp = rsaf_tmpp(rsaf_tmpp ~= 50);
    radial_axis = linspace(-0.5*bf_.nFOV, 0.5*bf_.nFOV, numel(rsaf_tmpp));
    plot(radial_axis, rsaf_tmpp/max(rsaf_tmpp), 'LineWidth', 2);
    
    hold off;
    xlabel('axial [mm]'); ylabel('1d psf');
    xlim([-10 10]);
    legend('CON', 'rSAF-foc5', 'rSAF-foc10', 'rSAF-foc15', 'rSAF-foc20', 'Location', 'northwest');
    set(gcf, 'Position', [696 320 612 369]);
end
