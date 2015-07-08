function f = simplify(f, tol)
%SIMPLIFY  Remove small trailing Chebyshev coeffs of a happy CHEBTECH object.
%  G = SIMPLIFY(F) attempts to compute a 'simplified' version G of the happy
%  CHEBTECH object F such that LENGTH(G) <= LENGTH(F) but ||G - F|| is small in
%  a relative sense: ||G - F|| < G.EPSLEVEL*G.VSCALE. It does this by removing
%  trailing coefficients of F that are relatively small; more precisely, those
%  that are smaller in magnitude than the product of F.VSCALE and F.EPSLEVEL.
%  G.EPSLEVEL is set to F.EPSLEVEL.
%
%  If F is not happy, F is returned unchanged.
%
%  G = SIMPLIFY(F, TOL) does the same as above but uses TOL instead of
%  F.EPSLEVEL as the relative threshold level for deciding whether a coefficient
%  is small enough to be removed. Here, G.EPSLEVEL is set to the maximum of
%  F.EPSLEVEL and TOL.
%
% See also HAPPINESSCHECK.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    return
end

% Do nothing to an unhappy CHEBTECH:
if ( ~f.ishappy )
    return
end

% Grab coefficients
coeffs = f.coeffs;
[n,m] = size(coeffs);

% Use the chebfunpref.eps if no tolerance was supplied:
if ( nargin < 2 )
    p = chebfunpref;
    tol = p.eps;
end
if length(tol) ~= m
    tol = max(tol)*ones(1,m);
end

% extend to 17 if necessary
if ( n < 17 )
    coeffs = [coeffs;ones(17-n,1)*min(min(abs(coeffs), [],1),tol)];
end

% Loop through columns to compute cutoff
cutoff = 1;
for k = 1:m
    cutoff = max(cutoff,standardChop(coeffs(:,k),tol(k)));
end
cutoff = min(cutoff,n);

% Chop coefficients
f.coeffs = coeffs(1:cutoff,:);

% Make sure epslevel is eps
f.epslevel = eps + 0*f.epslevel;

end
