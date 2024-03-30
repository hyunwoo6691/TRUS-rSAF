% %% loop - Resolution
% close all; clc;
% 
% clearvars -except stVars
% 
% dir_master = '/Users/songhyunwoo/Documents/JohnsHopkins/Research/Elevational SA/01.Data/fine_focal_point/Data';
% 
% if(exist('stVars'))
%     sample_start = stVars.sample_start;
%     error_start = stVars.error_start;
%     folder_start = stVars.folder_start;
% else % if excute first time
%     sample_start = 1;
%     error_start = 1;
%     folder_start = 1;
% end
% 
% Loop_Check_resolution_fine_takeall

%% loop - Resolution
% close all; clc;
% 
% clearvars -except stVars
% 
% dir_master = '/Users/songhyunwoo/Documents/JohnsHopkins/Research/Elevational SA/01.Data/fine_focal_point/Data';
% 
% if(exist('stVars'))
%     sample_start = stVars.sample_start;
%     error_start = stVars.error_start;
%     folder_start = stVars.folder_start;
% else % if excute first time
%     sample_start = 1;
%     error_start = 1;
%     folder_start = 1;
% end
% 
% Loop_Check_resolution_fine_specific_depth
%% loop - XCORR
% close all; clc;
% 
% clearvars -except stVars
% 
% % dir_master = '/Users/songhyunwoo/Documents/JohnsHopkins/Research/Elevational SA/01.Data/fine_focal_point';
% dir_master = '/Volumes/Research/eSAF/Data_psf/fine_focal_point';
% 
% if(exist('stVars'))
%     sample_start = stVars.sample_start;
%     error_start = stVars.error_start;
%     folder_start = stVars.folder_start;
% else % if excute first time
%     sample_start = 1;
%     error_start = 1;
%     folder_start = 1;
% end
% 
% Loop_Check_xcorr_fine
%% single
% close all; clc; clear all;
% 
% dir_ = '/Users/songhyunwoo/Documents/JohnsHopkins/Research/Elevational SA/01.Data/fine_focal_point/[elevational]Focus20mm_Tx2048';
% 
% Check_resolution_fine
