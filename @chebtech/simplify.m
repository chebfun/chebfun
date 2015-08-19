function f = simplify(f, tol)
%SIMPLIFY  Remove small trailing Chebyshev coeffs of a happy CHEBTECH object.
%  G = SIMPLIFY(F) attempts to compute a 'simplified' version G of the happy
%  CHEBTECH object F such that LENGTH(G) <= LENGTH(F) but ||G - F|| is small in
%  a relative sense. It does this by calling the routine STANDARDCHOP.
%
%  If F is not happy, F is returned unchanged.
%
%  G = SIMPLIFY(F, TOL) does the same as above but uses TOL instead of
%  EPS. 
%
% See also STANDARDCHOP.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case.
if ( isempty(f) )
    return
end

% Set F.EPSLEVEL to be MATLAB EPS.
f.epslevel = eps + 0*f.epslevel;

% Do nothing to an unhappy CHEBTECH.
if ( ~f.ishappy )
    return
end

% STANDARDCHOP requires at least 17 coefficients, so for F such that LENGTH(F) <
% 17, the coefficients are padded using prolong. The following
% parameters are chosen explicitly to work with STANDARDCHOP.
% See STANDARDCHOP for details.
nold = length(f);
N = max(17, round(nold*1.25 + 5));
f = prolong(f,N);

% Grab coefficients of F.
coeffs = f.coeffs;
[n, m] = size(coeffs);
coeffs = chebtech2.vals2coeffs(chebtech2.coeffs2vals(coeffs));

% Use CHEBFUNPREF.EPS if no tolerance was supplied.
p = chebfunpref;
if ( nargin < 2 )
    tol = p.eps;
end

% Recast TOL as a row vector.
if ( size(tol, 2) ~= m )
    tol = max(max(tol), p.eps)*ones(1, m);
end

% Loop through columns to compute CUTOFF.
cutoff = 1;
for k = 1:m
    cutoff = max(cutoff, standardChop(coeffs(:,k), tol(k)));
end

% Take the minimum of CUTOFF and LENGTH(F). This is necessary when padding was
% required.
cutoff = min(cutoff, nold);

% Chop coefficients using CUTOFF.
f.coeffs = coeffs(1:cutoff,:);

end
