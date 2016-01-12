% Test file for singfun/restrict.m

function pass = test_restrict(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% The order of the exponents:
a = 0.64;
b = -0.64;
c = 1.28;
d = -1.28;

% Fractional root at the left endpoint and the smooth part has no roots in 
% [-1 1]. We restrict f to a subinterval which is far enough from the endpoint
% singularities or roots.

op  = @(x) (1+x).^a.*exp(x);
data.exponents = [a 0];
data.singType = {'sing', 'none'};
f = singfun(op, data, pref);
pass(1) = test_spotcheck_restrict(f, op, [-0.2 0.1], [0, 0], pref);

% fractional pole at the left endpoint and note that the left endpoint is not a
% root, even though the smooth part has a root there. We strict f to a
% subinterval which includes the pole at the left endpoint.

op = @(x) (1+x).^d.*sin(50*pi*x);
data.exponents = [d+1 0];
data.singType = {'sing', 'none'};
f = singfun(op, data, pref);
pass(2) = test_spotcheck_restrict(f, op, [-1 0.3], [1, 0], pref);

% fractional root at the right endpoint and the smooth part has no roots in 
% [-1 1]. We restrict f to multiple subintervals.

op  =@(x) (1-x).^c.*cos(x);
data.exponents = [0 c];
data.singType = {'none', 'root'};
f = singfun(op, data, pref);
pass(3) = test_spotcheck_restrict(f, op, [-1 -0.7 1], [0, 0], pref);

% fractional pole at the right endpoint. We restrict f to multiple
% subintervals.

op = @(x) (1-x).^b;
data.exponents = [0 b];
data.singType = {'none', 'sing'};
f = singfun(op, data, pref);
pass(4) = test_spotcheck_restrict(f, op, [-0.9 -0.3 0.7 1], [0, 1], pref);

% Two fractional poles at both endpoints. We restrict f to multiple
% subintervals.

op = @(x) (1+x).^b.*sin(x).*(1-x).^d;
data.exponents = [b d];
data.singType = {'sing', 'sing'};
f = singfun(op, data, pref);
pass(5) = test_spotcheck_restrict(f, op, [-1 -0.9 0.5 0.7 1], [1, 1], pref);

% Check the case with roots close to endpoints.
p = 1e-4;
op = @(x) (1+x).^b.*sin(x).*(1-x).^(3*c);
data.exponents = [b b];
data.singType = {'sing', 'root'};
f = singfun(op, data, pref);
pass(6) = test_spotcheck_restrict(f, op, [-1+p 1-p], [0, 0], pref);

end

% Spot-check restriction of a given function to a given subinterval.
function result = test_spotcheck_restrict(f, fun_op, subint, shunlr, pref)

% Perform restriction:
g = restrict(f, subint);

% Number of subintervals:
NumInts = length(subint) - 1;

% Preallocate the storage for the result in each subinterval:
result = zeros(1, NumInts);

% If there's only one subinterval, then put it in a cell.
if ( NumInts == 1 )
    g = {g};
end

% Set the length of the step we step away from the singularities:
step = 1e-4;

% Loop over each subinterval:
for j = 1: NumInts
    
    % Construct mapping from restricted subinterval to [-1, 1]:
    a = subint(j);
    b = subint(j+1);
    map = @(t) (2/(b - a))*(t - a) - 1;
    
    % Sample on a grid of 100 points and check for accuracy:
    % Note that 'shunlr' tells us if we should step away from the endpoints of 
    % the original domain due to infinite values.
    if ( j == 1 && shunlr(1) == 1 )
        x = linspace(a + step, b, 100).';
    elseif ( j == NumInts && shunlr(2) == 1 )
        x = linspace(a, b - step, 100).';
    else
        x = linspace(a, b, 100).';
    end
    
    % Compute the exact function values and the approximated ones:
    y_exact = fun_op(x);
    y_approx = feval(g{j}, map(x));
    
    % The absolute value of the values
    absy = abs(y_exact);
    absy(absy < 1) = 1;
    
    % check the error
    err = norm(abs(y_exact - y_approx), inf);
    tol = norm(5e4*absy*eps, inf);
    result(j) = err < tol;
end

result = all(result);

end
