function F = getCoeffs2ValsTransform(~)
%GETCOEFFS2VALSTRANSFORM   Returns the IFFT in 1D.
%
% See also SPINOP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = @(u) ifft(u);
    
end