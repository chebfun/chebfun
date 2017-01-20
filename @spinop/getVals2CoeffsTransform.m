function F = getVals2CoeffsTransform(~)
%GETVALS2COEFFSTRANSFORM   Returns the FFT in 1D.
%
% See also SPINOP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = @(u) fft(u);
    
end