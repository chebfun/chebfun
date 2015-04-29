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

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    return
end

% Do nothing to an unhappy CHEBTECH:
if ( ~f.ishappy )
    return
end

% Use the f.epslevel if no tolerance was supplied:
if ( nargin < 2 )
    tol = f.epslevel;
end

% Until July 2014 we used to zero interior coefficients as well
% as trailing ones with the following code:
% f.coeffs(bsxfun(@minus, abs(f.coeffs), tol.*f.vscale) < 0) = 0;
% Check for trailing zero coefficients:
% [ignored, firstNonZeroRow] = find(f.coeffs.' ~= 0, 1);

% Check for trailing coefficients smaller than the tolerance relative
% to F.VSCALE:
largeCoeffs = (bsxfun(@minus, abs(f.coeffs), tol.*f.vscale) > 0);
[ignored, lastNonZeroRow] = find(largeCoeffs.' == 1, 1, 'last');

% If the whole thing is now zero, leave just one coefficient:
if ( isempty(lastNonZeroRow) )
    lastNonZeroRow = 1;
    f.coeffs = 0*f.coeffs;
end

% Remove trailing zeros:
if ( lastNonZeroRow > 0 )
    f.coeffs = f.coeffs(1:lastNonZeroRow, :);
end

% Update epslevel:
f.epslevel = max(f.epslevel, tol);

end
