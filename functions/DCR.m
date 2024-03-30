function [mDCROut, Fil] = DCR( mBFedData , stMID, stRFInfo)
switch stMID.nDCRType
    case 'bandpass'
        Passband = [stMID.nDCRF1 stMID.nDCRF2];
        PassbandW = Passband/stRFInfo.nFs*2;
    case 'high'
        Passband = stMID.nDCRFcut;
        PassbandW = Passband/stRFInfo.nFs*2;
end

Fil = fir1(stMID.nDCRTap,PassbandW, stMID.nDCRType);

mDCROut = convn(mBFedData, Fil','same');

end
