% Test file for @chebfun/rdivide.m.

function pass = test_rdivide(pref)

if ( nargin < 1 )
    pref = chebpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Check the empty cases.
f = chebfun(@sin, pref);
g = chebfun();
pass(1) = isempty(f./g) && isempty(g./f);

% Check zero-numerator case.
pass(2) = iszero(0./f);

% Check a few simple examples.
f = chebfun(@(x) exp(x), [-1 -0.5 0 0.5 1], pref);
g = 1./f;
g_exact = @(x) exp(-x);
pass(3) = norm(feval(g, x) - g_exact(x), inf) < 10*vscale(g)*epslevel(g);

f = chebfun(@(x) exp(x), [-1 -0.5 0 0.5 1], pref);
g = chebfun(@(x) exp(-x), [-1 1], pref);
h = f./g;
h_exact = @(x) exp(2*x);
pass(4) = norm(feval(h, x) - h_exact(x), inf) < 10*vscale(h)*epslevel(h);

f = chebfun(@(x) [exp(x) exp(-x)], [-1 -0.5 0 0.5 1], pref);
g = 1./f;
g_exact = @(x) [exp(-x) exp(x)];
err = feval(g, x) - g_exact(x);
pass(5) = norm(err(:), inf) < 10*vscale(g)*epslevel(g);

f = chebfun(@(x) [exp(x) exp(-x)], [-1 -0.5 0 0.5 1], pref);
g = chebfun(@(x) [exp(-x) exp(x)], [-1 1], pref);
h = f./g;
h_exact = @(x) [exp(2*x) exp(-2*x)];
err = feval(h, x) - h_exact(x);
pass(6) = norm(err(:), inf) < 10*vscale(h)*epslevel(h);

ft = f.';
gt = g.';
h = ft./gt;
h_exact = @(x) [exp(2*x) exp(-2*x)].';
err = feval(h, x) - h_exact(x);
pass(7) = norm(err(:), inf) < 10*vscale(h)*epslevel(h);

% Check error conditions.
try
    h = f./0;
    pass(8) = false;
catch ME
    pass(8) = strcmp(ME.identifier, 'CHEBFUN:rdivide:DivisionByZero');
end

try
    h = f./chebfun(@(x) 0*x);
    pass(9) = false;
catch ME
    pass(9) = strcmp(ME.identifier, 'CHEBFUN:rdivide:DivisionByZeroChebfun');
end

try
    f./gt;
    pass(10) = false;
catch ME
    pass(10) = strcmp(ME.identifier, 'CHEBFUN:rdivide:dim');
end

try
    f = chebfun(@(x) exp(x), [-1 1]);
    g = chebfun(@(x) exp(x), [0 2]);
    h = f./g;
    pass(11) = false;
catch ME
    pass(11) = strcmp(ME.identifier, 'CHEBFUN:rdivide:domain');
end

%% Integration with singfun:
dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

pow1 = -0.5;
pow2 = -0.3;
op1 = @(x) (x - dom(2)).^pow1.*sin(x);
op2 = @(x) (x - dom(2)).^pow2.*(cos(x).^2+1);
pref.singPrefs.exponents = [0 pow1];
f = chebfun(op1, dom, pref);
pref.singPrefs.exponents = [0 pow2];
g = chebfun(op2, dom, pref);
h = f./g;
vals_h = feval(h, x);
pow = pow1-pow2;
op = @(x)  (x - dom(2)).^pow.*(sin(x)./(cos(x).^2+1));
h_exact = op(x);
pass(12) = ( norm(vals_h-h_exact, inf) < 1e1*max(get(f, 'epslevel'), get(g, 'epslevel'))*...
    norm(h_exact, inf) );

end
