function mM = minandmax2est(f, N)
%MINANDMAX2EST   Estimates the minimum and maximum of a SEPARABLEAPPROX.
%   mM = MINANDMAX2EST(F) returns estimates for the minimum and maximum of the
%   SEPARABLEAPPROX F over its domain.  mM is a vector of length 2 such that
%   mM(1) is the estimated minimum and mM(2) is the estimated maximum.
%
%   mM = MINANDMAX2EST(F, N) returns estimates for the minimum and maximum of
%   the SEPARABLEAPPROX F over its domain, based on the evaluation on an N by N
%   Chebyshev grid in the domain of F (N = 32 by default).
%
% See also MINANDMAX2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    mM = [0, 0];
    return
end

if ( ( nargin < 2 ) || isempty(N) )
    % Default to N = 32:
    N = 32;
end

% Create N by N Chebyshev grid:
[XX, YY] = chebpts2(N, N, f.domain);
% Evaluate f on the grid:
vals = feval(f, XX, YY);
% Make result a column vector:
vals = vals(:);

mM = [ min(vals), max(vals) ];

end