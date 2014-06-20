function f = simplify(f, tol)
%SIMPLIFY  Zero small Chebyshev coefficients of a happy CHEBTECH object.
%  G = SIMPLIFY(F) attempts to compute a 'simplified' version G of the happy
%  CHEBTECH object F such that LENGTH(G) <= LENGTH(F) but ||G - F|| is small in
%  a relative sense: ||G - F|| < G.EPSLEVEL*G.VSCALE. It does this by zeroing
%  out all coefficients of F that are relatively small; more precisely, it sets
%  to zero all coefficients smaller in magnitude than the product of F.VSCALE
%  and F.EPSLEVEL. It then removes all trailing zero coefficients from F if
%  there are any. G.EPSLEVEL is set to F.EPSLEVEL.
%
%  If F is not happy, F is returned unchanged.
%
%  G = SIMPLIFY(F, TOL) does the same as above but uses TOL instead of the
%  F.EPSLEVEL as the relative threshold level for deciding whether a coefficient
%  is small enough to be zeroed. Here, G.EPSLEVEL is set to the maximum of
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

% Zero all coefficients smaller than the tolerance relative to F.VSCALE:
f.coeffs(bsxfun(@minus, abs(f.coeffs), tol.*f.vscale) < 0) = 0;

% Check for trailing zero coefficients:
[ignored, firstNonZeroRow] = find(f.coeffs.' ~= 0, 1);

% If the whole thing's now zero, leave just one coefficient:
if ( isempty(firstNonZeroRow) )
    firstNonZeroRow = size(f, 1);
end

% Remove trailing zeros:
if ( firstNonZeroRow > 0 )
    f.coeffs = f.coeffs(firstNonZeroRow:end, :);
end

% Update values and epslevel:
f.epslevel = max(f.epslevel, tol);

end
