function [ishappy, epslevel, cutoff] = standardCheck(f, values, vscl, pref)
%STANDARDCHECK   Attempt to trim trailing Chebyshev coefficients in a CHEBTECH.
%   [ISHAPPY, EPSLEVEL, CUTOFF] = STANDARDCHECK(F) uses the routine STANDARDCHOP
%   to compute a positive integer CUTOFF which represents the number of
%   coefficients of F that are deemed accurate enough to keep. ISHAPPY is TRUE
%   if the CUTOFF value returned by STANDARDCHOP is less than LENGTH(F) and
%   FALSE otherwise. EPSLEVEL is always returned as MATLAB EPS.
%
%   [ISHAPPY, EPSLEVEL, CUTOFF] = STANDARDCHECK(F, VALUES, VSCL, PREF) allows
%   additional preferences to be passed. VALUES is a matrix of the function
%   values of F at the corresponding interpolation points. VSCL is an
%   approximation of the maximum function value of F on a possibly larger
%   approximation interval. PREF is a data structure used to pass in additional
%   information, e.g. a target accuracy tolerance could be passed using
%   PREF.EPS.
%
% See also CLASSICCHECK, STRICTCHECK, LOOSECHECK.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Grab the coefficients of F.
coeffs = f.coeffs;
[n,m] = size(coeffs);

% Initialize ISHAPPY to be FALSE.
ishappy = false;

% Initialize EPSLEVEL to be MATLAB EPS.
epslevel = eps*ones(1,m); 

% Initialize CUTOFF to be LENGTH(F).
cutoff = length(f);

% Check for NaNs and exit if any are found.
if ( any(isnan(coeffs)) )
    error('CHEBFUN:CHEBTECH:standardCheck:nanEval', ...
        'Function returned NaN when evaluated.')
end

% If a PREF object is passed use the supplied value of EPS.
% Otherwise, use the CHEBTECH default EPS.
if ( nargin == 1 )
    pref = f.techPref();
    tol = pref.eps;
elseif ( isnumeric(pref) )
    tol = pref;
else
    tol = pref.eps;
end

% Convert TOL to a row vector.
if ( size(tol,2) ~= m )
  tol = ones(1,m)*max(tol);
end

% Compute function values of F if none were given.
if ( nargin < 2 || isempty(values) )
    values = f.coeffs2vals(f.coeffs);
end

% Scale TOL by the MAX(||F||*F.HSCALE, VSCL);
nrmf = max(abs(values), [], 1);
if ( isempty(vscl) )
    vscl = nrmf;
end
%tol = tol.*max(nrmf.*f.hscale, vscl);
tol = tol.*vscl./nrmf;
%ratio = vscl./nrmf

% Loop through columns of coeffs.
ishappy = false(1,m);
cutoff = zeros(1,m);
for k = 1:m

    % Call STANDARDCHOP.
    cutoff(k) = standardChop(coeffs(:,k), tol(k));

    % Check for happiness.
    ishappy(k) = ( cutoff(k) < n );

    % Exit if any column is unhappy.
    if ( ~ishappy(k) )
        break
    end

end

% Set outputs.
ishappy = all(ishappy); 
cutoff = max(cutoff);

end


