function [ishappy, cutoff] = standardCheck(f, values, data, pref)
%STANDARDCHECK   Attempt to trim trailing Fourier coefficients in a TRIGTECH.
%   [ISHAPPY, CUTOFF] = STANDARDCHECK(F) uses the routine STANDARDCHOP to
%   compute a positive integer CUTOFF which represents the number of
%   coefficients of F that are deemed accurate enough to keep.  ISHAPPY is TRUE
%   if the CUTOFF value returned by STANDARDCHOP is less than LENGTH(F) and
%   FALSE otherwise.
%
%   [ISHAPPY, CUTOFF] = STANDARDCHECK(F, VALUES, DATA, PREF) allows additional
%   preferences to be passed. VALUES is a matrix of the function values of F at
%   the corresponding interpolation points. DATA.VSCALE is an approximation of
%   the maximum function value of F on a possibly larger approximation
%   interval.  PREF is a data structure used to pass in additional information,
%   e.g. a target accuracy tolerance could be passed using PREF.EPS.
%
% See also CLASSICCHECK, STRICTCHECK, LOOSECHECK.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Grab the coefficients of F.
coeffs = abs(f.coeffs(end:-1:1,:));
[n, m] = size(coeffs);

% Compute some values if none were given.
if ( nargin < 2 || isempty(values) )
    values = f.coeffs2vals(f.coeffs);
end

% In order to work with STANDARDCHOP, the coefficients of F are modified so that
% the entries corresponding to wave numbers k and -k appear sequentially in the
% new matrix of coefficients. These entries are also replaced by the sum of the
% absolute values of the k and -k coefficients.

% Need to handle odd/even cases separately.
isEven = mod(n, 2) == 0;
if ( isEven )
    coeffs = [coeffs(n,:) ; coeffs(n-1:-1:n/2+1,:) + coeffs(1:n/2-1,:) ; coeffs(n/2,:)];
else
    coeffs = [coeffs(n:-1:(n+1)/2+1,:) + coeffs(1:(n+1)/2-1,:) ; coeffs((n+1)/2,:)];
end
coeffs = flipud(coeffs);
coeffs = [coeffs(1,:) ; kron(coeffs(2:end,:), [1 ; 1])];

% Initialize ISHAPPY.
ishappy = false;

% Initialize CUTOFF.
cutoff = n; 

% NaNs are not allowed.
if ( any(isnan(coeffs)) )
    error('CHEBFUN:TRIGTECH:standardCheck:nanEval', ...
        'Function returned NaN when evaluated.')
end

% Grab some preferences.
if ( nargin == 1 )
    pref = f.techPref();
    tol = pref.eps;
elseif ( isnumeric(pref) )
    tol = pref;
else
    tol = pref.eps;
end

% Reshape TOL.
if ( size(tol, 2) ~= m )
  tol = ones(1, m)*max(tol);
end

% Scale TOL by VSCL/||F||;
nrmf = max(abs(values), [], 1);
tol = tol.*data.vscale./nrmf;

% Loop through columns of coeffs
ishappy = false(1, m);
cutoff = zeros(1, m);
for k = 1:m

    % Call STANDARDCHOP.
    cutoff(k) = standardChop(coeffs(:,k), tol(k));

    % Check for happiness.
    ishappy(k) = ( cutoff(k) < n );

    % Divide CUTOFF by 2.
    if ( mod(cutoff(k), 2) == 0 )
        cutoff(k) = cutoff(k)/2;
    else
        cutoff(k) = (cutoff(k) - 1)/2;
    end

    % Break if unhappy.
    if ( ~ishappy(k) )
        break
    end

end

% Set outputs.
ishappy = all(ishappy); 

% CUTOFF is always odd.
cutoff = 2*max(cutoff) + 1;

end
