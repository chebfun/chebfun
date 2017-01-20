function F = getCoeffs2ValsTransform(~)
%GETCOEFFS2VALSTRANSFORM   Returns the IFFTN in 3D.
%
% See also SPINOP3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = @(u) ifftn(u);
    
end