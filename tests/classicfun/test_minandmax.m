% Test file for @classicfun/minandmax.m

function pass = test_minandmax(pref)

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
pass(1) = test_spotcheck_minmax(@(x) sin(10*x), data, -1, 1, pref);
pass(2) = test_spotcheck_minmax(@(x) real(airy(x)), data, 7.492128863997157e-07, 0.535656656015700, ...
    pref);
pass(3) = test_spotcheck_minmax(@(x) -1./(1 + x.^2), data, -1, ...
    -0.02, pref);
pass(4) = test_spotcheck_minmax(@(x) (x/10).^3.*cosh(x/10), ...
    data, (-.2)^3*cosh(-.2), 0.7^3*cosh(0.7), pref);

%%
% Check operation for array-valued inputs.
fun_op = @(x) [sin(10*x) real(airy(x)) (x/10).^3.*cosh(x/10)];
f = bndfun(fun_op, data, pref);
[y, x] = minandmax(f);
y_exact = [-1 7.492128863997157e-07  (-.2)^3*cosh(-.2);
    1 0.535656656015700 0.7^3*cosh(0.7)];
pass(5) = all(abs(y(:) - y_exact(:)) < 10*max(get(f, 'epslevel')));

% Check that the points x are indeed extreme points of the function 
% operator.
pass(6) = 1;
for k = 1:1:size(f, 2)
    fx = fun_op(x(:, k));
    max(abs(fx(:, k)) - y_exact(:, k));
    if ( max(abs(fx(:, k) - y_exact(:, k))) > 10*get(f, 'epslevel') )
        pass(6) = 0;
        break;
    end
end

%%  
% Test complex-array-valued BNDFUN objects.
f = bndfun(@(x) [exp(sin(2*x)), 1i*cos(20*x)], data);
[vals, pos] = minandmax(f);
f1 = bndfun(@(x) exp(sin(2*x)), data);
[vals1, pos1] = minandmax(f1);
f2 = bndfun(@(x) 1i*cos(20*x), data);
[vals2, pos2] = minandmax(f2);
pass(7) = norm(abs(vals) - abs([vals1 vals2]), inf) < ...
    10*max(get(f, 'vscale').*get(f, 'epslevel'));
    
%% 
% Test on singular BNDFUN.
pow = -0.5;
op = @(x) (x - data.domain(1)).^pow.*(sin(x).^2);
singData = data;
singData.exponents = [pow 0];
pass(8) = test_spotcheck_minmax(op, singData, 0, Inf, singPref);

%%
% Test complex-array-valued BNDFUN.
f = bndfun(@(x) [exp(sin(2*x)), 1i*cos(20*x)], data);
[vals, pos] = minandmax(f);
f1 = bndfun(@(x) exp(sin(2*x)), data);
[vals1, pos1] = minandmax(f1);
f2 = bndfun(@(x) 1i*cos(20*x), data);
[vals2, pos2] = minandmax(f2);
pass(9) = norm(abs(vals) - abs([vals1 vals2]), inf) < ...
    10*max(get(f, 'vscale').*get(f, 'epslevel'));
    
%% Tests for UNBNDFUN:

% Doubly-infinite domain:

% Set the domain:
data.domain = [-Inf Inf];

op = @(x) (1-exp(-x.^2))./x;
f = unbndfun(op, data);
[vals, pos] = minandmax(f);
% These exact solutions are obtained using Mathematica:
vExact = [-0.6381726863389515 ; 0.6381726863389515];
pExact = [-1.120906422778534 ; 1.120906422778534];
errV = vals - vExact;
errP = pos - pExact;
pass(10) = ( norm(errV, inf) < get(f,'epslevel')*get(f,'vscale') ) && ...
    ( norm(errP, inf) < 1e1*get(f,'epslevel')*get(f,'vscale') );

end

%% 
% Spot-check the results for a given BNDFUN.
function result = test_spotcheck_minmax(fun_op, data, exact_min, ...
    exact_max, pref)

    f = bndfun(fun_op, data, pref);
    [y, x] = minandmax(f);
    y_exact = [exact_min ; exact_max];
    fx = fun_op(x);
    result = ((max(abs(y - y_exact)) < 10*get(f, 'epslevel')) && ... 
        (max(abs(fx - y_exact)) < 10*get(f, 'epslevel')));
end
