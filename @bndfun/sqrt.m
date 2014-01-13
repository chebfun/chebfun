function f = sqrt(f)
%SQRT   Square root of a BNDFUN.
%   SQRT(F) returns the square root of a BNDFUN F.
%
%   This function assumes that the curve traced out by F in the complex plane
%   both (1) does not come too close to zero except at the domain boundaries 
%   +/- 1 and (2) does not cross over the branch cut in SQRT along the negative
%   real axis.  That is, F should not vanish at any point of the interior of its
%   domain, and the imaginary part of F should not vanish at any point of the
%   interior of its domain where the real part of F is negative.  If any of
%   these assumptions are violated, garbage may be returned with no warning.
%
% See also POWER.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If there are roots at the end of the domain, then make the f.onefun a singfun:
lval = get(f, 'lval');                          % Value at left of domain.
rval = get(f, 'rval');                          % Value at right of domain.
tol = 100*get(f, 'epslevel').*get(f, 'vscale'); % Tolerance for a root.

% Whether F has a vanishing value at any end point determines which SQRT() we'll
% call, the SMOOTHFUN one or the SINGFUN one.

if ( any(abs(lval) < tol) || any(abs(rval) < tol) ) && ...
        ( ~isa(f.onefun, 'singfun') )
    f.onefun = singfun(f.onefun);               % Cast f.onefun to a SINGFUN.
end

% Call SQRT of the ONEFUN:
f.onefun = sqrt(f.onefun);

end