function [ishappy, epslevel, cutOff] = standardCheck(f, values, vscl, pref)
%STANDARDCHECK  This function is a wrapper for Nick's standardChop
%  routine to chop a series of Chebyshev coefficients, see below.

% Grab the coefficients:
coeffs = f.coeffs;
[n,m] = size(coeffs);

%% initialize ishappy
ishappy = false;

%% initialize epslevel
epslevel = eps*ones(1,m); 

% NaNs are not allowed.
if ( any(isnan(coeffs)) )
    error('CHEBFUN:CHEBTECH:standardCheck:nanEval', ...
        'Function returned NaN when evaluated.')
end

% Compute some values if none were given:
%if ( nargin < 2 || isempty(values) )
%    values = f.coeffs2vals(f.coeffs);
%end
%if ( isempty(vscl) || isempty(vscl) )
%    vscl = max(abs(values), [],  1);
%end
%vscl = max(vscl ,max(abs(values), [],  1));

% Set the function scaling for each vector of values.
maxvals = max(abs(values), [], 1);

% set tolerance
tol = eps;
tol = pref.eps;

%% Loop through columns of coeffs
ishappy = false(1,m);
cutOff = zeros(1,m);
for k = 1:m
    [ishappy(k), cutOff(k)] = standardChop(coeffs(:,k), tol, vscl(k));
    if ( ~ishappy(k) )
        % No need to continue if it fails on any column.
        break
    end
end

ishappy = all(ishappy); 
cutOff = max(cutOff);

end


