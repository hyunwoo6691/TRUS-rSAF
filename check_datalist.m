clc; clear; close all;

%%
dir_ = uigetdir('./data','');
dir_ = [dir_ '/Raw data'];

%%
folder_list = dir(dir_);
flag = 0;
if(strcmp(folder_list(3).name,'.DS_Store')), flag = 1; end
folder_list = folder_list(3+flag:end);

%%
list = 1:1:4095;

%%
list_mask = zeros(numel(list),1);
file_mask = zeros(numel(list),1);
for k = 1:numel(folder_list)
    name_tmp = split(folder_list(k).name,'_');
    list_mask(str2double(name_tmp{end})) = 1;
   
    % check if the file exists
    if(list_mask(str2double(name_tmp{end})))
        file_list = dir([dir_ '/' folder_list(k).name]);
        if(numel(file_list)==2), continue; end
        flag = 0;
        if(strcmp(file_list(3).name,'.DS_Store')), flag = 1; end
        file_list = file_list(3+flag:end);
        if(strcmp(file_list(end).name,'Scanline_001.mat'))
            file_mask(str2double(name_tmp{end})) = 1;
        end
    end
end

%%
find(list_mask == 0)
find(file_mask == 0)