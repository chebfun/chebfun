function f = simplify(f, tol)
%SIMPLIFY  Zero small Fourier coefficients of a happy FOURTECH object.
%  G = SIMPLIFY(F) attempts to compute a 'simplified' version G of the happy
%  FOURTECH object F such that LENGTH(G) <= LENGTH(F) but ||G - F|| is small in
%  a relative sense: ||G - F|| < G.EPSLEVEL*G.VSCALE. It does this by zeroing
%  out all coefficients of F that are relatively small; more precisely, it sets
%  to zero all coefficients smaller in magnitude than the product of F.VSCALE
%  and the default FOURTECH EPS preference. It then removes all trailing zero
%  coefficients from F if there are any. G.EPSLEVEL is set to the maximum of
%  F.EPSLEVEL and the default FOURTECH EPS.
%
%  If F is not happy, F is returned unchanged.
%
%  G = SIMPLIFY(F, TOL) does the same as above but uses TOL instead of the
%  default FOURTECH EPS preference as the relative threshold level for deciding
%  whether a coefficient is small enough to be zeroed. Here, G.EPSLEVEL is set
%  to the maximum of F.EPSLEVEL and TOL.
%
% See also HAPPINESSCHECK.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    return
end

% Do nothing to an unhappy FOURTECH. ([TODO]: Is this the right thing to do?)
if ( ~f.ishappy )
    return;
end

% Use the default tolerance if none was supplied:
if ( nargin < 2 )
    tol = f.epslevel/2;
end

c = f.coeffs;  % Obtain Fourier coefficients {c_k}
numCoeffs = size(c, 1);
fIsEven = ( mod(numCoeffs, 2) == 0 );

% Split the coefficients into the positive and negative Fourier modes.
if ( fIsEven )
    % In this case the negative cofficients have an additional term
    % corresponding to the cos(N/2*x) coefficient. We account for this by
    % making the positive coefficients symmetric.
    cn = c(numCoeffs:-1:numCoeffs/2,:);
    cn(1,:) = 0.5*c(numCoeffs,:);
    cp = [cn(1,:); c(1:numCoeffs/2,:)];
else
    cn = c(numCoeffs:-1:(numCoeffs+1)/2,:);
    cp = c(1:(numCoeffs+1)/2,:);
end

% Need to check both the positive and negative coeficients in the Fourier
% expansion.

% Zero all coefficients smaller than the tolerance relative to F.VSCALE:
idp = bsxfun(@minus, abs(cp), tol.*f.vscale) < 0;
idn = bsxfun(@minus, abs(cn), tol.*f.vscale) < 0;
cp(idp) = 0;
cn(idn) = 0;

% Check for trailing zero coefficients:
[ignored, firstNonZeroRowP] = find(cp.' ~= 0, 1);
[ignored, firstNonZeroRowN] = find(cn.' ~= 0, 1);

% If the whole thing's now zero, leave just one coefficient:
if ( isempty(firstNonZeroRowP) && isempty(firstNonZeroRowN) )
    firstNonZeroRow = size(cp, 1);
% The negative and positive cofficient vectors need to be the same length
% So, we remove the smaller of the tails from both.
else
    firstNonZeroRow = min(firstNonZeroRowP, firstNonZeroRowN);
end

% Remove trailing zeros:
if ( firstNonZeroRow > 0 )
    cp = cp(firstNonZeroRow:end,:);
    cn = cn(firstNonZeroRow:end,:);
end

% Now put the coefficients vector back together.
f.coeffs = [cp(1:end,:); cn(end-1:-1:1,:)];

% Update values and epslevel:
f.values = f.coeffs2vals(f.coeffs);
f.vscale = max(abs(f.values), [], 1);
f.epslevel = max(f.epslevel, tol);

end