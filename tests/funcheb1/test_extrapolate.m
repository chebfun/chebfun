function pass = test_extrapolate(varargin)

% Obtain preferences
if ( nargin == 1 )
    pref = varargin{1};
else
    pref = funcheb.pref;
end
    
% Set a tolerance:
tol = 100*pref.funcheb.eps;

% Initialise:
pass = zeros(1,4);
s = 0;

n = 17;
x = funcheb1.chebpts(n);

for k = 1:2
    % Test once for a vector of coeffs and against for a matrix.

    % Interior NaN;
    values = sin(x-x(4))./(x-x(4));
    newValues = funcheb1.extrapolate(values);
    pass(1 + s) = all( abs(newValues(4,:) - 1) < tol );
    
    % Interior Inf;
    values = sin(x-x(4)-eps)./(x-x(4));
    newValues = funcheb1.extrapolate(values);
    pass(2 + s) = all( abs(newValues(4,:) - 1) < tol );
    
    % Make x a matrix and repeat:
    x = repmat(x, 1, 2);
    s = length(pass)/2;
    
end

end
