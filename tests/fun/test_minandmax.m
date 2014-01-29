% Test file for fun/minandmax.m

function pass = test_minandmax(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebpref();
end

% Set a domain for BNDFUN.
dom = [-2 7];

%%
% Spot-check the extrema for a few BNDFUN.
pass(1) = test_spotcheck_minmax(@(x) sin(10*x), dom, -1, 1, pref);
pass(2) = test_spotcheck_minmax(@(x) real(airy(x)), dom, 7.492128863997157e-07, 0.535656656015700, ...
    pref);
pass(3) = test_spotcheck_minmax(@(x) -1./(1 + x.^2), dom, -1, ...
    -0.02, pref);
pass(4) = test_spotcheck_minmax(@(x) (x/10).^3.*cosh(x/10), ...
    dom, (-.2)^3*cosh(-.2), 0.7^3*cosh(0.7), pref);

%%
% Check operation for array-valued inputs.
fun_op = @(x) [sin(10*x) real(airy(x)) (x/10).^3.*cosh(x/10)];
f = bndfun(fun_op, dom, [], [], pref);
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
f = bndfun(@(x) [exp(sin(2*x)), 1i*cos(20*x)], dom);
[vals, pos] = minandmax(f);
f1 = bndfun(@(x) exp(sin(2*x)), dom);
[vals1, pos1] = minandmax(f1);
f2 = bndfun(@(x) 1i*cos(20*x), dom);
[vals2, pos2] = minandmax(f2);
pass(7) = norm(abs(vals) - abs([vals1 vals2]), inf) < ...
    10*max(get(f, 'vscale').*get(f, 'epslevel'));
    
%% 
% Test on singular BNDFUN.
pow = -0.5;
op = @(x) (x - dom(1)).^pow.*(sin(x).^2);
pref.singPrefs.exponents = [pow 0];
pass(8) = test_spotcheck_minmax(op, dom, 0, Inf, pref);

%%
% Test complex-array-valued BNDFUN.
f = bndfun(@(x) [exp(sin(2*x)), 1i*cos(20*x)], dom);
[vals, pos] = minandmax(f);
f1 = bndfun(@(x) exp(sin(2*x)), dom);
[vals1, pos1] = minandmax(f1);
f2 = bndfun(@(x) 1i*cos(20*x), dom);
[vals2, pos2] = minandmax(f2);
pass(9) = norm(abs(vals) - abs([vals1 vals2]), inf) < ...
    10*max(get(f, 'vscale').*get(f, 'epslevel'));
    
%% Tests for UNBNDFUN:

% Doubly-infinite domain:

% Set the domain:
dom = [-Inf Inf];

op = @(x) (1-exp(-x.^2))./x;
f = unbndfun(op, dom);
[vals, pos] = minandmax(f);
% These exact solutions are obtained using Mathematica:
vExact = [-0.6381726863389515 ; 0.6381726863389515];
pExact = [-1.120906422778534 ; 1.120906422778534];
errV = vals - vExact;
errP = pos - pExact;
pass(10) = ( norm(errV, inf) < get(f,'epslevel')*get(f,'vscale') ) && ...
    ( norm(errP, inf) < get(f,'epslevel')*get(f,'vscale') );

end

%% 
% Spot-check the results for a given BNDFUN.
function result = test_spotcheck_minmax(fun_op, dom, exact_min, ...
    exact_max, pref)

    f = bndfun(fun_op, dom, [], [], pref);
    [y, x] = minandmax(f);
    y_exact = [exact_min ; exact_max];
    fx = fun_op(x);
    result = ((max(abs(y - y_exact)) < 10*get(f, 'epslevel')) && ... 
        (max(abs(fx - y_exact)) < 10*get(f, 'epslevel')));
end