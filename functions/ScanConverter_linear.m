function [aXAxis, aZAxis, mOutput] = ScanConverter_linear(mInput, pitch, unit_dis, nHeight, nWidth, dz, dx)

    [nSample, nScanline] = size(mInput); % nInXlen: sample num, nInYlen:scanline num
    
    input_width = pitch*(nScanline-1)*0.5; % one side
    input_depth = (nSample-1)*unit_dis;

    aXAxis = -dx*(nWidth/2-0.5) : dx : dx*(nWidth/2-0.5);
    aZAxis = (0 : dz : dz*(nHeight-1));

    [mZ, mX] = ndgrid(aZAxis, aXAxis); % (x, z) coordinates
    
    mFov = (mX>(-input_width))&(mX<input_width)&(mZ>0)&(mZ<input_depth);

    mOutput = ones(nHeight, nWidth)*1e-3; % Converted image (resultant image)
    for zidx = 1:nHeight
       for xidx = 1:nWidth
           % Only within imaging region
           if(mFov(zidx,xidx)==1)
               
               sidx     = mZ(zidx,xidx)/unit_dis*2+ 1;% sample index  %  (mR(zidx,xidx))/dr + 1; 
               sidx_int = floor(sidx); % integer component of sample index
               sidx_fr  = sidx - floor(sidx); % fractional component of sample index
               
               cidx 	= (mX(zidx,xidx)+input_width)/pitch + 1; % scanline index
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