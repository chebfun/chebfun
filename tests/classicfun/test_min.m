% Test file for @classicfun/min.m

function pass = test_min(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

singPref = pref;
singPref.blowup = true;

% Set a domain for BNDFUN.
data.domain = [-2 7];
    
%%
% Spot-check the extrema for a few BNDFUN.
pass(1) = test_spotcheck_min(@(x) -sin(10*x), data, -1, pref);
pass(2) = test_spotcheck_min(@(x) -real(airy(x)), data, ...
    -0.535656656015700, pref);
pass(3) = test_spotcheck_min(@(x) 1./(1 + x.^2), data, ...
    0.02, pref);
pass(4) = test_spotcheck_min(@(x) -(x/10).^3.*cosh(x/10), ...
    data, -0.7^3*cosh(0.7), pref);
    
%%
% Check operation for array-valued BNDFUN inputs.
fun_op = @(x) -[sin(10*x) real(airy(x)) (x/10).^3.*cosh(x/10)];
f = bndfun(fun_op, data, pref);
[y, x] = min(f);
exact_max = -[1 0.535656656015700 0.7^3*cosh(0.7)];
fx = -[sin(10*x(1)) airy(x(2)) (x(3)/10).^3.*cosh(x(3)/10)];
tol = 10*get(f, 'vscale').*get(f, 'epslevel');
pass(5) = (all(abs(y - exact_max) < tol) && ...
    all(abs(fx - exact_max) < tol));
    
%%
% Test for complex-valued BNDFUN.
pass(6) = test_spotcheck_min( ...
    @(x) (x/2).*(exp(1i*(x/2))+1i*sin(x/2)), data, ...
    0, pref);
fun_op = @(x) [((x-2).^2/4+1).*exp(1i*(x/2)) ... 
    -((x+1).^2/4+1).*exp(1i*(x/2))];
f = bndfun(fun_op, data, pref);
[y, x] = min(f);
exact_min = [exp(1i) -exp(-1i/2)];
fx = fun_op(x); 
fx = fx([1 4]);
tol = 10*max(get(f, 'vscale').*get(f, 'epslevel'));
pass(7) = (all(abs(y - exact_min) < 10*tol) && ...
    all(abs(fx - exact_min) < 10*tol));

%% Test for UNBNDFUN:

% Functions on [-inf b]:

% Set the domain:
data.domain = [-Inf -3*pi];

% A blow-up function:
op = @(x) x.*(5+exp(x.^3))./(data.domain(2)-x);
singData = data;
singData.exponents = [0 -1];
f = unbndfun(op, singData, singPref);
[y, x] = min(f);
yExact = -Inf;
xExact = data.domain(2);
errX = x - xExact;
pass(8) = ( norm(errX, inf) < get(f,'epslevel')*get(f,'vscale') ) && ...
    ( y == yExact );

end

%%
% Spot-check the results for a given BNDFUN.
function result = test_spotcheck_min(fun_op, data, exact_min, pref)

f = bndfun(fun_op, data, pref);
[y, x] = min(f);
fx = fun_op(x);
err1 = abs(y - exact_min);
err2 = abs(fx - exact_min);
tol = 2e3*get(f, 'vscale')*get(f, 'epslevel');
result = (err1 < tol) && (err2 < tol);

end
