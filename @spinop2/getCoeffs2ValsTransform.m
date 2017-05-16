function F = getCoeffs2ValsTransform(~)
%GETCOEFFS2VALSTRANSFORM   Returns the IFFT2 in 2D.
%
% See also SPINOP2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = @(u) ifft2(u);
    
end