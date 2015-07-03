function [ishappy, epslevel, cutoff] = standardCheck(f, values, vscl, pref)
%STANDARDCHECK  This function is a wrapper for Nick's standardChop
%  routine to chop a series of Chebyshev coefficients, see below.

% Grab the coefficients:
coeffs = abs(f.coeffs(end:-1:1,:));
[n,m] = size(coeffs);

% Need to handle odd/even cases separately.
isEven = ~mod(n, 2);
if isEven
    % In this case the negative cofficients have an additional term
    % corresponding to the cos(N/2*x) coefficient.
    coeffs = [coeffs(n,:);coeffs(n-1:-1:n/2+1,:)+coeffs(1:n/2-1,:);coeffs(n/2,:)];
else
    coeffs = [coeffs(n:-1:(n+1)/2+1,:)+coeffs(1:(n+1)/2-1,:);coeffs((n+1)/2,:)];
end
coeffs = flipud(coeffs);

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
    shift = vscl;
end
shift = shift.*f.hscale;
shift = 1./max(shift,1);

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
cutoff = 2*max(cutoff)+1;

end


