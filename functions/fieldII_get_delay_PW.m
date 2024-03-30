function delay_curve = fieldII_get_delay_PW(t_idx, aElePosX, sound_speed, stTxInfo)
tx_angle = stTxInfo.tx_angles(t_idx);

dist_tmp = aElePosX * tand(tx_angle);
dist_tmp = dist_tmp + abs(min(dist_tmp));

delay_curve = dist_tmp / sound_speed;

end

