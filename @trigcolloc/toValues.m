function fx = toValues(disc, f, varargin)
%TOVALUES   Convert CHEBFUN to a TRIGCOLLOC discretization.
%   TOVALUES(DISC, F) converts a chebfun F to values at equi-spaced points 
%   in the DISC.DOMAIN.
%
% See also TOFUNCTION.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isnumeric(f) )
    fx = f;
    return
end

x = functionPoints(disc);
fx = f(x);

end