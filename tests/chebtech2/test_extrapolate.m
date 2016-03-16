% Test file for chebtech2/extrapolate.m

function pass = test_extrapolate(pref)

% Obtain preferences
if ( nargin < 1 )
    pref = chebtech.techPref();
end
    
% Set a tolerance:
tol = 100*pref.chebfuneps;

% Initialise:
s = 0;

n = 17;
x = chebtech2.chebpts(n);
f = chebtech2;

for k = 1:2
    % Test once for a vector of coeffs and against for a matrix.

    % Interior NaN;
    values = sin(x)./x;
    newValues = f.extrapolate(values);
    pass(1 + s) = all( abs(newValues((n+1)/2,:) - 1) < tol );
    
    % Interior Inf;
    values = sin(x-eps)./x;
    newValues = f.extrapolate(values);
    pass(2 + s) = all( abs(newValues((n+1)/2,:) - 1) < tol );

    % Extrapolate left end:
    values = sign(x+1);
    values(1, :) = NaN;
    newValues = f.extrapolate(values);
    pass(3 + s) = all( abs(newValues(1,:) - 1) < tol );
    
    % Extrapolate right end:
    values = sign(x-1);
    values(end, :) = NaN;
    newValues = f.extrapolate(values);
    pass(4 + s) = all( abs(newValues(end,:) + 1) < tol );

    % Check reverting of endpoint values:
    values = sin(x);
    newValues = f.extrapolate(values);
    pass(5 + s) = ~any(values(:) - newValues(:));
    
    values = sin(x)./x;
    newValues = f.extrapolate(values);
    pass(6 + s) = ~any(values(:) - newValues(:));
   
    % Make x a matrix and repeat:
    x = repmat(x, 1, 2);
    s = length(pass)/2;
    
end

end
