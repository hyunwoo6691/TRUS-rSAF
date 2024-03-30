for k = 1:numel(scanline_data)
    name_tmp = scanline_data(k).name;
    components_ = split(name_tmp,'_');
    idx_ = str2double(components_{end});
    
    new_name = [components_{1} '_' num2str(idx_+1024,'%.4d')];
%     
%     if(strcmp(name_tmp, new_name))
%         continue;
%     else
        movefile([dir_data '/' name_tmp], [dir_data_p '/' new_name]);
%     end
end