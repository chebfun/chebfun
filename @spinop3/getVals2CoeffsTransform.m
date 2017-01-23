function F = getVals2CoeffsTransform(~)
%GETVALS2COEFFSTRANSFORM   Returns the FFTN in 3D.
%
% See also SPINOP3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = @(u) fftn(u);
    
end