function out = iszero(f, tol)
%ISZERO   Check if a CHEBFUN is identically zero on its domain.
%   ISZERO(F) returns true if F is identically zero or empty on F.domain and
%   false otherwise. If F is an array-valued CHEBFUN, the a true/false value is
%   returned for each column.

% TODO:  Document the TOL input.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% An empty CHEBFUN is zero:
if ( isempty(f) )
    out = true;
    return
end

% Choose a tolerance:
if ( nargin < 2 )
    tol = vscale(f)*epslevel(f);
    % TODO: Remove this once epslevels of zero CHEBFUNs have been improved.
    if ( isnan(tol) )
        tol = chebfun.pref('eps');
    end
end

% Impulses:
% TODO:  What about higher-order impulses?
out = all(f.impulses(:,:,1) <= tol);

% Loop over each of the FUNs:
% TODO:  We don't use TOL here?
k = 0;
while ( k < numel(f.funs) && any(out) )
    k = k + 1;
    out = out & iszero(f.funs{k});
end

end
