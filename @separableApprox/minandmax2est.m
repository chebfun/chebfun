function mM = minandmax2est(f, N)
%MINANDMAX2EST   Estimates the minimum and maximum of a SEPARABLEAPPROX.
%   mM = MINANDMAX2EST(F) returns estimates for the minimum and maximum of the
%   SEPARABLEAPPROX F over its domain.  mM is a vector of length 2 such that
%   mM(1) is the estimated minimum and mM(2) is the estimated maximum.
%
%   mM = MINANDMAX2EST(F, N) returns estimates for the minimum and maximum of
%   the SEPARABLEAPPROX F over its domain, based on samples on an N by N grid
%   (N = 33 by default).
%
% See also MINANDMAX2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    mM = [0, 0];
    return
end

if ( ( nargin < 2 ) || isempty(N) )
    % Default to N = 33:
    N = 33;
end

% Sample f on an appropriate grid:
vals = sample(f, N, N);

% Make result a column vector:
vals = vals(:);

% Get min and max:
mM = [ min(vals), max(vals) ];

end
