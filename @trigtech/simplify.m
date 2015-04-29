function f = simplify(f, tol)
%SIMPLIFY  Remove small trailing Fourier coefficients of a happy TRIGTECH object.
%  G = SIMPLIFY(F) attempts to compute a 'simplified' version G of the happy
%  TRIGTECH object F such that LENGTH(G) <= LENGTH(F) but ||G - F|| is small in
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

% Do nothing to an unhappy TRIGTECH. ([TODO]: Is this the right thing to do?)
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
    % In this case the negative coefficients have an additional term
    % corresponding to the cos(N/2*x) coefficient. We account for this by
    % making the positive coefficients symmetric.
    numModes = numCoeffs/2+1;
    cn = c(numModes:-1:1,:);
    cn(numModes,:) = 0.5*cn(numModes,:);
    cp = [c(numModes:numCoeffs,:); cn(numModes,:)];
else
    numModes = (numCoeffs+1)/2;
    cp = c(numModes:numCoeffs,:);
    cn = c(numModes:-1:1,:);
end

% Need to check both the positive and negative coefficients in the Fourier
% expansion.

% Check for trailing coefficients smaller than the tolerance relative
% to F.VSCALE:
idp = bsxfun(@minus, abs(cp), tol.*f.vscale) > 0;
idn = bsxfun(@minus, abs(cn), tol.*f.vscale) > 0;

% Before July 2014 we used to zero all small coefficients:
% cp(idp) = 0;
% cn(idn) = 0;
% Check for trailing zero coefficients:
% [ignored, firstNonZeroRowP] = find(cp.' ~= 0, 1);
% [ignored, firstNonZeroRowN] = find(cn.' ~= 0, 1);

% Check for trailing small coefficients:
[ignored, lastNonZeroRowP] = find(idp.' == 1, 1, 'last');
[ignored, lastNonZeroRowN] = find(idn.' == 1, 1, 'last');

% If the whole thing's now zero, leave just one coefficient:
if ( isempty(lastNonZeroRowP) && isempty(lastNonZeroRowN) )
    lastNonZeroRowP = 1;
    lastNonZeroRowN = 1;
    cp = 0*cp; 
    cn = 0*cn;
end

lastNonZeroRow = max(lastNonZeroRowP, lastNonZeroRowN);

% Remove trailing zeros:
if ( lastNonZeroRow > 0 )
    cp = cp(1:lastNonZeroRow,:);
    cn = cn(1:lastNonZeroRow,:);
end

% Now put the coefficients vector back together.
f.coeffs = [cn(end:-1:2,:); cp(1:end,:)];

% Update values and epslevel:
f.values = f.coeffs2vals(f.coeffs);
f.vscale = max(abs(f.values), [], 1);
f.epslevel = max(f.epslevel, tol);

end
