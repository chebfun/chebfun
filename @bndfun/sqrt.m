function f = sqrt(f)
%SQRT   Square root of a BNDFUN.
%   SQRT(F) returns the square root of a BNDFUN F. Note, it is assumed that the
%   only roots of F are located at the endpoints of F.domain.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% If there are roots at the end of the domain, then make the f.onefun a singfun:
lval = get(f, 'lval');                         % Value at left of domain.
rval = get(f, 'rval');                         % Value at right of domain.
tol = 100*get(f, 'epslevel')*get(f, 'vscale'); % Tolerance for a root.
if ( any(abs(lval) < tol) || any(abs(rval) < tol) )
    f.onefun = singfun(f.onefun);              % Cast f.onefun to a SINGFUN.
end

% Call SQRT() of the ONEFUN:
f.onefun = sqrt(f.onefun);

end