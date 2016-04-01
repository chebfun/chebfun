% Test file for chebtech1/extrapolate.m

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
x = chebtech1.chebpts(n);
f = chebtech1;

for k = 1:2
    % Test once for a vector of coeffs and against for a matrix.

    % Interior NaN;
    values = sin(x-x(4))./(x-x(4));
    newValues = f.extrapolate(values);
    pass(1 + s) = all( abs(newValues(4,:) - 1) < tol );
    
    % Interior Inf;
    values = sin(x-x(4)-eps)./(x-x(4));
    newValues = f.extrapolate(values);
    pass(2 + s) = all( abs(newValues(4,:) - 1) < tol );
    
    % Make x a matrix and repeat:
    x = repmat(x, 1, 2);
    s = length(pass)/2;
    
end

end
