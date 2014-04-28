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
    pref = fourtech.techPref();
    tol = f.epslevel.*f.vscale;
    % TODO: Document this.
%     vscale = f.vscale;
%     vscale(vscale < f.epslevel) = 1;
%     tol = max(tol)./vscale;
end

c = f.coeffs;  % Obtain Fourier coefficients {c_k}
numCoeffs = size(c,1);
fIsEven = mod(numCoeffs,2) == 0;
% Split the coefficients into the positive and negative Fourier modes
if fIsEven
    % In this case the positive cofficients have an additional term
    % corresponding to the cos(N/2*x) coefficient. We account for this by
    % making the positive coefficients symmetric 
    cp = c(numCoeffs:-1:numCoeffs/2,:);
    cp(1,:) = 0.5*c(numCoeffs,:);
    cn = [cp(1,:);c(1:numCoeffs/2,:)];
else
    cp = c(numCoeffs:-1:(numCoeffs+1)/2,:);
    cn = c(1:(numCoeffs+1)/2,:);
end

% Only need to work with the coefficients corresponding to the positve
% modes of the Fourier expansion;

% Zero all coefficients smaller than the tolerance relative to F.VSCALE:
id = bsxfun(@minus, abs(cp), tol.*f.vscale) < 0;
cp(id) = 0;
cn(id) = 0;

% Check for trailing zero coefficients:
[ignored, firstNonZeroRow] = find(cp.' ~= 0, 1);

% If the whole thing's now zero, leave just one coefficient:
if ( isempty(firstNonZeroRow) )
    firstNonZeroRow = length(cp);
end

% Remove trailing zeros:
if ( firstNonZeroRow > 0 )
    cp = cp(firstNonZeroRow:end,:);
    cn = cn(firstNonZeroRow:end,:);
end

% % Now put the coefficients vector back together.
% if fIsEven
%     c = [cn(end-1:-1:2,:);cp(1:end,:)];
%     c(end,:) = 2*c(end,:);
% else
%     c = [cn(end:-1:2,:);cp(1:end,:)];
% end
f.coeffs = [cn(1:end-1,:);cp(end:-1:1)];

% Update values and epslevel:
f.values = f.coeffs2vals(f.coeffs);
f.vscale = max(abs(f.values), [], 1);
f.epslevel = max(f.epslevel, tol);

end
