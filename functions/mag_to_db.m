function img_db = mag_to_db(img_mag)
    aROI = find(img_mag ~= 50);
    aOutlier = find(img_mag== 50);
    mOutput = zeros(size(img_mag));
    mOutput(aROI) = img_mag(aROI);
    img_db = mOutput;
    img_db = db(mOutput/max(mOutput(:)));
    img_db(aOutlier) = -30;
end

