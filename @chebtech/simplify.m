function f = simplify(f, tol)
%SIMPLIFY   Zero small Chebyshev coefficients of a happy CHEBTECH object.
%   G = SIMPLIFY(F) attempts to compute a 'simplified' version G of the happy
%   CHEBTECH object F such that LENGTH(G) <= LENGTH(F) but ||G - F|| <
%   F.EPSLEVEL. It does this by zeroing out all coefficients of F smaller in
%   magnitude than the default CHEBTECH EPS preference and then removing all
%   trailing zero coefficients. G.EPSLEVEL is set to the maximum of F.EPSLEVEL
%   and the default CHEBTECH EPS.
%
%   If F is not happy, F is returned unchanged.
%
%   G = SIMPLIFY(F, TOL) does the same as above but uses TOL as the threshold
%   level for deciding whether a coefficient is small enough to be zeroed.
%   Here, G.EPSLEVEL is set to the maximum of F.EPSLEVEL and TOL.
%
% See also HAPPINESSCHECK.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    return
end

% Do nothing to an unhappy CHEBTECH.  ([TODO]:  Is this the right thing to do?)
if ( ~f.ishappy )
    return;
end

% Use the default tolerance if none was supplied:
if ( nargin < 2 )
    pref = chebtech.pref();
    tol = pref.chebtech.eps;
end

% Zero all coefficients smaller than the tolerance:
f.coeffs(abs(f.coeffs) < tol) = 0;

% Remove zero coefficients from the tail:
[ignored, numCoeffsToRemove] = find(f.coeffs.' ~= 0, 1);
if ( numCoeffsToRemove > 1 )
    f = prolong(f, length(f) - numCoeffsToRemove + 1);
end

% Update epslevel:
f.epslevel = max(f.epslevel, tol);

end
