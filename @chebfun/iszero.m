function out = iszero(f, tol)
%ISZERO   True for zero CHEBFUN objects.
%   ISZERO(F) returns logical TRUE if F is approximately zero everywhere and
%   logical FALSE otherwise.  "Approximately zero" means that each of the FUNs
%   of F is zero (per their ISZERO() function), the values of F at each of its
%   breakpoints is smaller than VSCALE(F)*EPSLEVEL(F) in magnitude, and F has
%   no higher-order impulses.
%
%   ISZERO(F, TOL) does the same but uses TOL instead of VSCALE(F)*EPSLEVEL(F)
%   as the tolerance for checking the size of the breakpoint values of F.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Use default tolerance if none supplied.
if ( nargin < 2 )
    tol = vscale(f)*epslevel(f)
end

% Impulses:
out = all(all(abs(f.impulses(:,:,1)) < tol)) && (size(f.impulses, 3) <= 1);

% Loop over each of the funs:
% [TODO]:  We don't use tol here?
k = 0;
while ( (k < numel(f.funs)) && out )
    k = k + 1;
    out = iszero(f.funs{k});
end

end
