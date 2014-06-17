% Test file for bndfun/restrict.m

function pass = test_restrict(pref)

% Get preferences.
if (nargin < 1)
    pref = chebfunpref();
end

% Set the domain
dom = [-2 7];

%%
% Check behavior for empty inputs.
f = bndfun();
f = restrict(f, [-0.5 0.5]);
pass(1) = isempty(f);

%%
% Check behaviour for non-subinterval inputs.
f = bndfun(@(x) sin(x), struct('domain', dom), pref);
g = restrict(f, dom);
pass(2) = isequal(f, g);

try
    g = restrict(f, [dom(1) - 1, dom(2) + 1]); %#ok<NASGU>
    pass(3) = 0;
catch ME
    pass(3) = strcmp(ME.identifier, 'CHEBFUN:BNDFUN:restrict:badInterval');
end

try
    g = restrict(f, [-Inf, 1]); %#ok<NASGU>
    pass(4) = 0;
catch ME
    pass(4) = strcmp(ME.identifier, 'CHEBFUN:BNDFUN:restrict:badInterval');
end

try
    g = restrict(f, [-1 -0.25 0.3 0.1 1]); %#ok<NASGU>
    pass(5) = 0;
catch ME
    pass(5) = strcmp(ME.identifier, 'CHEBFUN:BNDFUN:restrict:badInterval');
end

%%
% Check whether restriction actually results in a BNDFUN on the correct domain.
g = restrict(f,[2, 3]);
pass(6) = all(g.domain == [2, 3]);

%%
% Spot-check a few functions
pass(7) = test_spotcheck_restrict(@(x) exp(x) - 1, dom, [-2 4], [], pref);
pass(8) = test_spotcheck_restrict(@(x) 1./(1 + x.^2), dom, [-0.7 0.9], [], ...
    pref);
pass(9) = test_spotcheck_restrict(@(x) cos(1e3*x), dom, [0.1 0.5], [], pref);
pass(10) = test_spotcheck_restrict(@(t) sinh(t*exp(2*pi*1i/6)), dom, ...
    [-0.4 1], [], pref);

%%
% Check multiple subinterval restriction.
f = bndfun(@(x) sin(x) + sin(x.^2), struct('domain', dom), pref);
g = restrict(f, [-1.7 2.3 6.8]);
h1 = restrict(f, [-1.7 2.3]);
h2 = restrict(f, [2.3 6.8]);
x = linspace(-1, 1, 100).';
err1 = norm(feval(g{1} - h1, x), inf);
err2 = norm(feval(g{2} - h2, x+4), inf);
tol = 10*get(f, 'epslevel');
pass(11) = err1 < tol && err2 < tol;

%%
% Check whether restriction actually results in a BNDFUN on the correct domain.
g = restrict(f,[2, 3, 5]);
pass(12) = all(g{1}.domain == [2, 3] & g{2}.domain == [3,5]);

%%
% Check operation for array-valued functions.
pass(13) = test_spotcheck_restrict(@(x) [sin(x) cos(x) exp(x)], dom, ...
    [-1 -0.7], [], pref);

f = bndfun(@(x) [sin(x) cos(x)], struct('domain', dom), pref);
g = restrict(f, [-0.6 0.1 1]);
h1 = restrict(f, [-0.6 0.1]);
h2 = restrict(f, [0.1 1]);
x = linspace(-.5, 0, 100).';
err1 = norm(feval(g{1} - h1, x), inf);
err2 = norm(feval(g{2} - h2, x+.7), inf);
tol = 10*max(get(f, 'epslevel'));
pass(14) = err1 < tol && err2 < tol;

%% Test on singular function:

pow = -0.5;
op = @(x) (x - dom(1)).^pow.*sin(x);
pref.blowup = true;
data.exponents = [pow 0];
pass(15) = test_spotcheck_restrict(op, dom, [-1 -0.7], data, pref);

end

% Spot-check restriction of a given function to a given subinterval.
function result = test_spotcheck_restrict(fun_op, dom, subint, data, pref)
% Perform restriction.
data.domain = dom;
f = bndfun(fun_op, data, pref);
g = restrict(f, subint);

a = subint(1);
b = subint(2);

% Sample on a grid of 100 points and check for accuracy.
x = linspace(a, b, 100).';
y_exact = fun_op(x);
y_approx = feval(g, x);

result = norm(y_exact - y_approx, Inf) < ...
    10*max(get(f, 'vscale').*get(f, 'epslevel'));
end
