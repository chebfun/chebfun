function [ishappy, epslevel, cutoff] = standardCheck(f, values, vscl, pref)
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

% Compute shift
if ( isempty(vscl) )
    shift = ones(1,m);
else
    shift = max(abs(values), [],  1)./max(vscl,1);
end
shift = max(shift.*f.hscale,1);
shift = 1./shift;

% Loop through columns of coeffs
ishappy = false(1,m);
cutoff = zeros(1,m);
for k = 1:m

    % call standardChop
    [cutoff(k)] = standardChop(coeffs(:,k), tol(k), shift(k));
    ishappy(k) = ( cutoff(k) < n );
    if ( ~ishappy(k) )
        % No need to continue if it fails on any column.
        break
    end
end

% set outputs
ishappy = all(ishappy); 
cutoff = max(cutoff);

end


