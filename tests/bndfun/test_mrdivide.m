% Test file for bndfun/mrdivide.m

function pass = test_mrdivide(pref)

% Get preferences.
if (nargin < 1)
    pref = chebfunpref();
end

% Set a domain
dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

% Random number to use as a scalar constant.
alpha = -0.194758928283640 + 0.075474485412665i;

%%
% Check division of a BNDFUN object by a numeric array.

f_op = @(x) [sin(x) cos(x)];
f = bndfun(f_op, struct('domain', dom), pref);
pass(1) = isnan(f / 0);

g = f / alpha;
g_exact = @(x) f_op(x) ./ alpha;
pass(2) = norm(feval(g, x) - g_exact(x), inf) < ...
    10*max(get(g, 'vscale')*eps);

% A "least-squares" case where the solution is obvious.
I = eye(2);
g = f / I;
err = g*I - f;
pass(3) = max(max(abs(feval(err, x)))) < ...
    1e2*max(get(g, 'vscale')*eps);

% A less trivial least-squares case for which we still know the answer.
A = [1 1];
g = f / A;
g_exact = @(x) (sin(x) + cos(x))/2;
pass(4) = norm(feval(g, x) - g_exact(x), inf) < ...
    1e2*max(get(g, 'vscale')*eps);

%%
% Check division of a numeric array by a BNDFUN object.

f = bndfun(@(x) sin(x), struct('domain', dom), pref);
g = alpha / f;
pass(5) = abs(innerProduct(f, g) - alpha) < ...
    10*get(g, 'vscale')*eps;

f = bndfun(@(x) [sin(2*pi*x) cos(2*pi*x)], struct('domain', dom), pref);
g = [1 1]/f;
g_exact = @(x) (2/9)*(sin(2*pi*x) + cos(2*pi*x));
pass(6) = norm(feval(g, x) - g_exact(x), inf) < ...
    1e2*max(get(g, 'vscale')*eps);
    

%%
% Check error conditions.

% Catch dimension mismatch errors.
try
    g = f / [1 2 3]; %#ok<NASGU>
    pass(7) = false;
catch ME
    pass(7) = strcmp(ME.identifier, 'CHEBFUN:BNDFUN:mrdivide:size');
end

% Can't do f/g if both f and g are BNDFUN objects.
try
    f = bndfun(@(x) sin(x), struct('domain', dom), pref);
    g = bndfun(@(x) cos(x), struct('domain', dom), pref);
    h = f / g; %#ok<NASGU>
    pass(8) = false;
catch ME
    pass(8) = strcmp(ME.identifier, ...
        'CHEBFUN:BNDFUN:mrdivide:bndfunDivBndfun');
end

% Can't call mrdivide on a bndfund and a non-BNDFUN or non-double
% object.
try
    f = bndfun(@(x) sin(x), struct('domain', dom), pref);
    g = f / true; %#ok<NASGU>
    pass(9) = false;
catch ME
    pass(9) = strcmp(ME.identifier, 'CHEBFUN:BNDFUN:mrdivide:badArg');
end

end
