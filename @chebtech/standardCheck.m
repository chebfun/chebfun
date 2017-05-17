function [ishappy, cutoff] = standardCheck(f, values, data, pref)
%STANDARDCHECK   Attempt to trim trailing Chebyshev coefficients in a CHEBTECH.
%   [ISHAPPY, CUTOFF] = STANDARDCHECK(F) uses the routine STANDARDCHOP to
%   compute a positive integer CUTOFF which represents the number of
%   coefficients of F that are deemed accurate enough to keep. ISHAPPY is TRUE
%   if the CUTOFF value returned by STANDARDCHOP is less than LENGTH(F) and
%   FALSE otherwise.
%
%   [ISHAPPY, CUTOFF] = STANDARDCHECK(F, VALUES, DATA, PREF) allows additional
%   preferences to be passed. VALUES is a matrix of the function values of F at
%   the corresponding interpolation points. DATA.VSCALE is an approximation of
%   the maximum function value of F on a possibly larger approximation
%   interval.  PREF is a data structure used to pass in additional information,
%   e.g. a target accuracy tolerance could be passed using PREF.CHEBFUNEPS.
%
% See also CLASSICCHECK, STRICTCHECK, LOOSECHECK.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Grab the coefficients of F.
coeffs = f.coeffs;
[n, m] = size(coeffs);

% Check for NaNs and exit if any are found.
if ( any(isnan(coeffs)) )
    error('CHEBFUN:CHEBTECH:standardCheck:nanEval', ...
        'Function returned NaN when evaluated.')
end

% If a PREF object is passed use the supplied value of EPS.
% Otherwise, use the CHEBTECH default EPS.
if ( nargin == 1 )
    pref = f.techPref();
    tol = pref.chebfuneps;
elseif ( isnumeric(pref) )
    tol = pref;
else
    tol = pref.chebfuneps;
end

% Convert TOL to a row vector.
if ( size(tol, 2) ~= m )
  tol = ones(1, m)*max(tol);
end

% Compute function values of F if none were given.
if ( nargin < 2 || isempty(values) )
    values = f.coeffs2vals(f.coeffs);
end

% Scale TOL by the MAX(DATA.HSCALE, DATA.VSCALE/VSCALE(F)).
% This choice of scaling is the result of undesirable behavior when using
% standardCheck to construct the function f(x) = sqrt(1-x) on the interval [0,1]
% with splitting turned on. Due to the way standardChop checks for plateaus, the
% approximations on the subdomains were chopped incorrectly leading to poor
% quality results. This choice of scaling corrects this by giving less weight to
% subintervals that are much smaller than the global approximation domain, i.e.
% HSCALE >> 1. For functions on a single domain with no breaks, this scaling has
% no effect since HSCALE = 1.
vscaleF = max(abs(values), [], 1);
tol = tol.*max(data.hscale, data.vscale./vscaleF);

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


