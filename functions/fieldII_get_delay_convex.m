function delay_curve = fieldII_get_delay_convex(apod_window, tx_focusX, tx_focusZ, aElePosX, aElePosZ, nC)
delay_curve = zeros(1, numel(apod_window));
mask_ = logical(apod_window);
active_eleX = aElePosX(mask_);
active_eleZ = aElePosZ(mask_);

dist_ = sqrt((active_eleX - tx_focusX).^2 + (active_eleZ - tx_focusZ).^2);

delay_tmp = max(dist_/nC) - (dist_ /nC);

delay_curve(mask_) = delay_tmp;
end

