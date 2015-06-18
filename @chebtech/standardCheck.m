function [ishappy, epslevel, cutOff] = standardCheck(f, op, values, vscl, pref)
%STANDARDCHECK  This function is a wrapper for Nick's standardChop
%  routine to chop a series of Chebyshev coefficients, see below.

% Grab the coefficients:
coeffs = f.coeffs;
[n,m] = size(coeffs);

%% initialize ishappy
ishappy = false;

%% initialize epslevel
epslevel = eps*ones(1,m); 

%% initialize tol
tol = eps;

% NaNs are not allowed.
if ( any(isnan(coeffs)) )
    error('CHEBFUN:CHEBTECH:standardCheck:nanEval', ...
        'Function returned NaN when evaluated.')
end

% Set the function scaling for each vector of values.
maxvals = max(abs(values), [], 1);

% Check the vertical scale:
if ( max(maxvals) == 0 )
    % This is the zero function, so we must be happy!
    ishappy = true;
    cutOff = 1;
    return
elseif ( any(isinf(maxvals)) )
    % Inf located. No cutoff.
    ishappy = false;
    cutOff = n;
    return
end

%% Loop through columns of coeffs
ishappy = false(1,m);
cutOff = zeros(1,m);
for k = 1:m
    [ishappy(k), cutOff(k)] = standardChop(coeffs(:,k), tol);
    if ( ~ishappy(k) )
        % No need to continue if it fails on any column.
        break
    end
end

ishappy = all(ishappy); 
cutOff = max(cutOff);

end


