% Test file for chebtech2/extrapolate.m

function pass = test_extrapolate(pref)

% Obtain preferences
if ( nargin < 1 )
    pref = chebtech.pref();
end
    
% Set a tolerance:
tol = 100*pref.eps;

% Initialise:
pass = zeros(1, 12);
s = 0;

n = 17;
x = chebtech2.chebpts(n);
f = chebtech2;

for k = 1:2
    % Test once for a vector of coeffs and against for a matrix.

    % Interior NaN;
    f.values = sin(x)./x;
    newValues = f.extrapolate();
    pass(1 + s) = all( abs(newValues((n+1)/2,:) - 1) < tol );
    
    % Interior Inf;
    f.values = sin(x-eps)./x;
    newValues = f.extrapolate();
    pass(2 + s) = all( abs(newValues((n+1)/2,:) - 1) < tol );

    % Extrapolate left end:
    f.values = sign(x+1);
    f.values(1, :) = NaN;
    newValues = f.extrapolate();
    pass(3 + s) = all( abs(newValues(1,:) - 1) < tol );
    
    % Extrapolate right end:
    f.values = sign(x-1);
    f.values(end, :) = NaN;
    newValues = f.extrapolate();
    pass(4 + s) = all( abs(newValues(end,:) + 1) < tol );

    % Check reverting of endpoint values:
    f.values = sin(x);
    newValues = f.extrapolate();
    pass(5 + s) = ~any(f.values(:) - newValues(:));
    
    f.values = sin(x)./x;
    newValues = f.extrapolate();
    pass(6 + s) = ~any(f.values(:) - newValues(:));
   
    % Make x a matrix and repeat:
    x = repmat(x, 1, 2);
    s = length(pass)/2;
    
end

end
