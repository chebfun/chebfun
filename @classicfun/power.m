function f = power(f, b)
%.^   CLASSICFUN power.
%   F.^G returns a CLASSICFUN F to the scalar power G, a scalar F to the
%   CLASSICFUN power G, or a CLASSICFUN F to the CLASSICFUN power G. F and or G
%   may be complex.
%
%   This function assumes that the curve traced out by F in the complex plane
%   both (1) does not come too close to zero except at the domain boundaries 
%   +/- 1 and (2) does not cross over the branch cut in POWER along the negative
%   real axis.  That is, F should not vanish at any point of the interior of its
%   domain, and the imaginary part of F should not vanish at any point of the
%   interior of its domain where the real part of F is negative.  If any of
%   these assumptions are violated, garbage may be returned with no warning.
%
%   H = POWER(F, G) is called for the syntax 'F .^ G'.
%
% See also SQRT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% If there are roots at the end of the domain, then make the f.onefun a singfun:
lval = get(f, 'lval');           % Value at left of domain.
rval = get(f, 'rval');           % Value at right of domain.
tol = 1e3*eps*get(f, 'vscale'); % Tolerance for a root.
if ( any(abs(lval) < tol) || any(abs(rval) < tol) ) && ...
        ( ~isa(f.onefun, 'singfun') )
    f.onefun = singfun(f.onefun);               % Cast f.onefun to a SINGFUN.
end

% Call POWER() of the ONEFUN:
f.onefun = power(f.onefun, b);

end
