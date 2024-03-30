% mInput            : Input data, size: nSample x nScanline
% dr                : Sample distance of mInput [m]
% da                : Scanline angle difference of mInput [deg]
% nTransRadius      : Radius of convex array transducer[m]
% nHeight           : Height of Output Image [pixel]
% nWidth            : Width of Output Image [pixel]
% dz                : Grid interval of Output Image [m]
% dx                : Grid interval of Output Image [m]
%
function [aXAxis, aZAxis, mOutput] = ScanConverter_convex(mInput, dr, da, nTransRadius, nHeight, nWidth, dz, dx)

    [nSample, nScanline] = size(mInput); % nInXlen: sample num, nInYlen:scanline num

    nHalfViewAngle = (nScanline-1)*da/2; % [degree] half of view angle

    aXAxis = -dx*(nWidth/2-0.5) : dx : dx*(nWidth/2-0.5);
    aZAxis = (0 : dz : dz*(nHeight-1)) + nTransRadius*cosd(nHalfViewAngle);

    [mZ, mX] = ndgrid(aZAxis, aXAxis); % (x, z) coordinates

    mR = sqrt(mZ.^2 + mX.^2); % (r, a) coordinates
    mA = atand(mX./mZ);

    % !!!!!!!!!!!!!!!!!!!!!!!!!!!! View depth is adjusted !!!!!!!!!!!!!!
    nEndDepth = dr*( (nSample)-1) + nTransRadius; % z-coordinate of end of image
    mFov = (mR>nTransRadius)&(mR<nEndDepth)&(mA>-nHalfViewAngle)&(mA<nHalfViewAngle); % Field of view (Image region)

%     figure;imagesc(mFov);
    
    mOutput = ones(nHeight, nWidth)*50; % Converted image (resultant image)
    for zidx = 1:nHeight
       for xidx = 1:nWidth
           % Only within imaging region
           if(mFov(zidx,xidx)==1)
               
               sidx     = (mR(zidx,xidx)-nTransRadius)/dr + 1;% sample index  %  (mR(zidx,xidx))/dr + 1; 
               sidx_int = floor(sidx); % integer component of sample index
               sidx_fr  = sidx - floor(sidx); % fractional component of sample index
               
               cidx 	= (mA(zidx,xidx)-(-nHalfViewAngle))/da + 1; % scanline index
               cidx_int = floor(cidx); % integer component of scanline index
               cidx_fr  = cidx - cidx_int; % fractional component of scanline index
               
               % Only when input sample exists
               if( (sidx_int < nSample) && (cidx_int < nScanline) )                   
                   nTmp1 = mInput(sidx_int  ,cidx_int  )*(1-sidx_fr)  +  mInput(sidx_int+1,cidx_int  )*sidx_fr;
                   nTmp2 = mInput(sidx_int  ,cidx_int+1)*(1-sidx_fr)  +  mInput(sidx_int+1,cidx_int+1)*sidx_fr;

                   mOutput(zidx,xidx) = nTmp1*(1-cidx_fr) + nTmp2*cidx_fr; 
               end

           end
       end
    end

end
