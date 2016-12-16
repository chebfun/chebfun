function mM = minandmax3est(f, N)
%MINANDMAX2EST   Estimates the minimum and maximum of a CHEBFUN3.
%   mM = MINANDMAX3EST(F) returns estimates for the minimum and maximum of the
%   CHEBFUN3 F over its domain.  mM is a vector of length 2 such that
%   mM(1) is the estimated minimum and mM(2) is the estimated maximum.
%
%   mM = MINANDMAX3EST(F, N) returns estimates for the minimum and maximum of
%   the CHEBFUN3 F over its domain, based on samples on an N by N by N grid
%   (N = 25 by default).
%
% See also CHEBFUN3/MINANDMAX3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    mM = [0, 0];
    return
end

if ( ( nargin < 2 ) || isempty(N) )
    % Default to N = 25:
    N = 25;
end

% Sample f on an appropriate grid:
vals = sample(f, N, N, N);

% Make result a column vector:
vals = vals(:);

% Get min and max:
mM = [ min(vals), max(vals) ];

end
