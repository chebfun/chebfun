function [ishappy, epslevel, cutOff] = standardCheck(f, values, vscl, pref)
%STANDARDCHECK  This function is a wrapper for Nick's standardChop
%  routine to chop a series of Chebyshev coefficients, see below.

% Grab the coefficients:
coeffs = f.coeffs;
[n,m] = size(coeffs);

% initialize ishappy
ishappy = false;

% initialize epslevel
epslevel = eps*ones(1,m); 

% NaNs are not allowed.
if ( any(isnan(coeffs)) )
    error('CHEBFUN:CHEBTECH:standardCheck:nanEval', ...
        'Function returned NaN when evaluated.')
end

% Grab some preferences:
if ( nargin == 1 )
    pref = f.techPref();
    tol = pref.eps;
elseif ( isnumeric(pref) )
    tol = pref;
else
    tol = pref.eps;
end
if size(tol,2) ~= m
  tol = ones(1,m)*max(tol);
end

% Compute some values if none were given:
if ( nargin < 2 || isempty(values) )
    values = f.coeffs2vals(f.coeffs);
end
if ( isempty(vscl) || isempty(vscl) )
    vscl = max(abs(values), [],  1);
end

% get hscale
hscale = f.hscale;

% Loop through columns of coeffs
ishappy = false(1,m);
cutOff = zeros(1,m);
for k = 1:m

    % call standardChop
    [cutOff(k)] = standardChop(coeffs(:,k), tol(k), vscl(k)*hscale);
    ishappy(k) = ( cutOff(k) < n );
    if ( ~ishappy(k) )
        % No need to continue if it fails on any column.
        break
    end
end

% set outputs
ishappy = all(ishappy); 
cutOff = max(cutOff);

end


