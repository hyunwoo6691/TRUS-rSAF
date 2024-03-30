function delay_curve = fieldII_get_delay(ele_idx, apod_window, tx_focus, aElePosX, nC)
delay_curve = zeros(1, numel(apod_window));
mask_ = logical(apod_window);
active_ele = aElePosX(mask_);

dist_ = sqrt((active_ele-aElePosX(ele_idx)).^2 + tx_focus^2);
delay_tmp = max(dist_/nC) - (dist_ /nC);

delay_curve(mask_) = delay_tmp;
end

