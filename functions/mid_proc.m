function env_out = mid_proc(input, mid_param, acoustic_param, bf_param)
[dcr_out, Fil] = DCR(input, mid_param, acoustic_param);
[tgc_out, aTGCCurve] = fDTGC(dcr_out, mid_param, acoustic_param, bf_param, size(dcr_out,1), acoustic_param.nUnitDis);
env_out = abs(hilbert(tgc_out));
end

