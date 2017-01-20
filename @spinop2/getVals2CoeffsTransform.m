function F = getVals2CoeffsTransform(~)
%GETVALS2COEFFSTRANSFORM   Returns the FFT2 in 2D.
%
% See also SPINOP2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = @(u) fft2(u);
    
end