function f = power(f, b)
% .^   BNDFUN power.
%   F.^G returns a BNDFUN F to the scalar power G, a scalar F to the BNDFUN
%   power G, or a BNDFUN F to the BNDFUN power G. F and or G may be complex. 
%   Note, it is assumed that F is non-zero on its domain. If F has zeros, then
%   the output is garbage without throwing a warning.
%
%   H = POWER(F, G) is called for the syntax 'F .^ G'.
%
% See also SQRT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% If there are roots at the end of the domain, then make the f.onefun a singfun:
lval = get(f, 'lval');                         % Value at left of domain.
rval = get(f, 'rval');                         % Value at right of domain.
tol = 1e2*get(f, 'epslevel')*get(f, 'vscale'); % Tolerance for a root.
if ( any(abs(lval) < tol) || any(abs(rval) < tol) ) && ...
        ( ~isa(f.onefun, 'singfun') )
    f.onefun = singfun(f.onefun);              % Cast f.onefun to a SINGFUN.
    f.onefun = extractBoundaryRoots(f.onefun);
end

% Call POWER() of the ONEFUN:
f.onefun = power(f.onefun, b);

end