clc; close all;
%%
dir_ = uigetdir('./data','');
%%
error_case = [0 0.1 0.2 0.5 1 2 5];
% scline = 554;
scline = round(130.9680/0.4724);

disorientations_ = zeros(numel(error_case), scline);
for e_idx = 1:numel(error_case)
    %%%%%%% RF data extract for beamforming %%%%%%%
    nError_range = error_case(e_idx);
    aError = normrnd(0, nError_range, [1, scline]); % Gaussian distribution, mean: 0, std: nError_range
    disorientations_(e_idx,:) = aError;
end

save([dir_ '/disorientations.mat'],'disorientations_');