function apod_window = fieldII_get_apod(ele_idx,tx_focus, tx_fnum, num_ele, pitch)
apod_window = zeros(1, num_ele);

num_tx_ele = ceil((tx_focus / tx_fnum)/pitch);
num_tx_ele = num_tx_ele + (1-mod(num_tx_ele,2)); % make it odd number
element_use = max(ele_idx-0.5*(num_tx_ele-1), 1):min(ele_idx+0.5*(num_tx_ele-1), num_ele);

apod_window(element_use) = 1;
end

