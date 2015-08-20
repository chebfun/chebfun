function fx = toValues(disc, f, varargin)
%TOVALUES   Convert CHEBFUN to a CHEBCOLLOC2 discretization.
%   TOVALUES(DISC,F) converts a chebfun F to values at 2nd kind points in the
%   DISC.DOMAIN.
%
% See also TOFUNCTION.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isnumeric(f) )
    fx = f;
    return
end

% The easy part.
x = functionPoints(disc);
fx = f(x);

% Evaluate left- and right-sided limits at breaks, which are part of the
% discretization.
n = disc.dimension;
csn = [0, cumsum(n)]; % Offsets into breakpoints
dxloc = csn(2:end-1);
fx(dxloc) = feval(f, x(dxloc), 'left');
fx(dxloc+1) = feval(f, x(dxloc), 'right');

end
