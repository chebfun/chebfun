% Test file for @chebfun/rdivide.m.

function pass = test_rdivide(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

%%

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

%% SCALAR-VALUED

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
pass(3) = norm(feval(g, x) - g_exact(x), inf) < 10*vscale(g)*eps;

f = chebfun(@(x) exp(x), [-1 -0.5 0 0.5 1], pref);
g = chebfun(@(x) exp(-x), [-1 1], pref);
h = f./g;
h_exact = @(x) exp(2*x);
pass(4) = norm(feval(h, x) - h_exact(x), inf) < 10*vscale(h)*eps;

%% ARRAY-VALUED

f = chebfun(@(x) [exp(x) exp(-x)], [-1 -0.5 0 0.5 1], pref);
g = 1./f;
g_exact = @(x) [exp(-x) exp(x)];
err = feval(g, x) - g_exact(x);
pass(5) = norm(err(:), inf) < 10*vscale(g)*eps;

f = chebfun(@(x) [exp(x) exp(-x)], [-1 -0.5 0 0.5 1], pref);
g = chebfun(@(x) [exp(-x) exp(x)], [-1 1], pref);
h = f./g;
h_exact = @(x) [exp(2*x) exp(-2*x)];
err = feval(h, x) - h_exact(x);
pass(6) = norm(err(:), inf) < 10*max(vscale(h)*eps);

ft = f.';
gt = g.';
h = ft./gt;
h_exact = @(x) [exp(2*x) exp(-2*x)].';
err = feval(h, x) - h_exact(x);
pass(7) = norm(err(:), inf) < 10*vscale(h)*eps;

%% QUASIMATRIX

f = chebfun(@(x) [exp(x) exp(-x)], [-1 -0.5 0 0.5 1], pref);
fq = quasimatrix(@(x) [exp(x) exp(-x)], [-1 -0.5 0 0.5 1], pref);
g = 1./fq;
g_exact = @(x) [exp(-x) exp(x)];
err = feval(g, x) - g_exact(x);
pass(8) = norm(err(:), inf) < 10*vscale(g)*eps;

g = chebfun(@(x) [exp(-x) exp(x)], [-1 1], pref);
gq = quasimatrix(@(x) [exp(-x) exp(x)], [-1 1], pref);
h_exact = @(x) [exp(2*x) exp(-2*x)];

h = fq./g;
err = feval(h, x) - h_exact(x);
pass(9) = norm(err(:), inf) < 10*vscale(h)*eps;

h = f./gq;
err = feval(h, x) - h_exact(x);
pass(10) = norm(err(:), inf) < 10*vscale(h)*eps;

h = fq./gq;
err = feval(h, x) - h_exact(x);
pass(11) = norm(err(:), inf) < 10*vscale(h)*eps;

fqt = fq.';
gqt = gq.';
h = fqt./gqt;
h_exact = @(x) [exp(2*x) exp(-2*x)].';
err = feval(h, x) - h_exact(x);
pass(12) = norm(err(:), inf) < 10*vscale(h)*eps;

%% Check error conditions.
try
    h = f./0;
    pass(13) = false;
catch ME
    pass(13) = strcmp(ME.identifier, ...
        'CHEBFUN:CHEBFUN:rdivide:columnRdivide:divisionByZero');
end

try
    h = chebfun(@(x) 1+0*x)./chebfun(@(x) 0*x);
    pass(13) = false;
catch ME
    pass(13) = strcmp(ME.identifier, ...
        'CHEBFUN:CHEBFUN:rdivide:columnRdivide:divisionByZeroChebfun');
end

try
    f./gt;
    pass(14) = false;
catch ME
    pass(14) = strcmp(ME.identifier, ...
        'CHEBFUN:CHEBFUN:rdivide:columnRdivide:dim');
end

try
    f = chebfun(@(x) exp(x), [-1 1]);
    g = chebfun(@(x) exp(x), [0 2]);
    h = f./g;
    pass(15) = false;
catch ME
    pass(15) = strcmp(ME.identifier, ...
        'CHEBFUN:CHEBFUN:rdivide:columnRdivide:domain');
end

%% Test on singular function: piecewise smooth chebfun - splitting on.
dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

pow1 = -0.5;
pow2 = -0.3;
op1 = @(x) (x - dom(2)).^pow1.*sin(100*x);
op2 = @(x) (x - dom(2)).^pow2.*(cos(300*x).^2+1);
f = chebfun(op1, dom, 'exps', [0 pow1], 'splitting', 'on');
g = chebfun(op2, dom, 'exps', [0 pow2], 'splitting', 'on');
h = f./g;
vals_h = feval(h, x);
pow = pow1-pow2;
op = @(x) (x - dom(2)).^pow.*(sin(100*x)./(cos(300*x).^2+1));
h_exact = op(x);
pass(16) = ( norm(vals_h-h_exact, inf) < 1e5*max(eps, ...
    eps)*norm(h_exact, inf) );


%% Test for function defined on unbounded domain:

% Functions on [2 inf]:

% Set the domain:
dom = [2 Inf];
domCheck = [2 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

opf = @(x) exp(-x.^2);
opg = @(x) x.^2;
oph = @(x) exp(-x.^2).*x.^-2;
f = chebfun(opf, dom);
g = chebfun(opg, dom, 'exps', [0 2]);
h = f./g;
hVals = feval(h, x);
hExact = oph(x);
err = hVals - hExact;
pass(17) = norm(err, inf) < 1e1*eps*get(f,'vscale');

%%
% Test that rdivide works correctly when the two triguns involved have 
% widely different vscales.
k = 10;
z = chebfun('4*exp(1i*s)',[0 2*pi],'trig');
f = z./((exp(z)-1));
B10 = factorial(k)*sum((f./z.^(k+1)).*diff(z))/(2i*pi);  % 10 Bernoulli num
exact = 5/66;
pass(18) = abs(exact-B10) < 1e2*eps*get(f,'vscale');


%% Test division between a CHEBFUN and a TRIGFUN.

dom = [0 pi 2*pi];

% 1. One column case.
f = chebfun(@(x) x + x.^2, dom, pref);
g = chebfun(@(x) 2+cos(x), [dom(1) dom(end)], 'periodic');
h1 = f./g;
% We want the result to use the same tech as the one used by f.
pass(19) = strcmpi(func2str(get(h1.funs{1}.onefun, 'tech')), ...
                   func2str(get(f.funs{1}.onefun, 'tech')));
h2 = chebfun(@(x) (x + x.^2)./(2+cos(x)), dom, pref);
pass(20) = norm(h1-h2, inf) < 1e2*eps*get(h2,'vscale');


% 2. Quasimatrix case.
f = chebfun(@(x) [cos(x), sin(x)], [dom(1) dom(end)], 'periodic');
g = chebfun(@(x) [2+x, 2+x.^3], dom, pref);
h1 = f./g;
% We want the result to use the same tech as the one used by g.
pass(21) = strcmpi(func2str(get(h1(:,1).funs{1}.onefun, 'tech')), ...
                   func2str(get(g(:,1).funs{1}.onefun, 'tech')));
pass(22) = strcmpi(func2str(get(h1(:,2).funs{1}.onefun, 'tech')), ...
                   func2str(get(g(:,2).funs{1}.onefun, 'tech')));
h2 = chebfun(@(x) [cos(x)./(2+x), sin(x)./(2+x.^3)], dom, pref);
pass(23) = norm(h1-h2, inf) < 1e2*eps*get(h2,'vscale');

%% Avoid introducing NaNs/Infs when dividing by zero:

y1 = chebfun('abs(sin(5*x))','splitting','on');
y2 = chebfun('abs(sin(4*x))','splitting','on');
y3 = y2 - y1;
g = (exp(2*y3) - exp(y3))./y3;
mg = merge(g);
pass(24) = (vscale(g) < inf) && (numel(g.funs) == 6);

%%
% #1111
try
    f = chebfun(@(x) exp(x));
    g = 0;
    f./g;
    pass(24) = false;
catch ME
    pass(24) = strcmp(ME.identifier, ...
        'CHEBFUN:CHEBFUN:rdivide:columnRdivide:divisionByZero');
end

try
    f = chebfun(@(x) [exp(x) exp(-x)]);
    g = [0 0];
    f./g;
    pass(25) = false;
catch ME
    pass(25) = strcmp(ME.identifier, ...
        'CHEBFUN:CHEBFUN:rdivide:columnRdivide:divisionByZero');
end

try
    f = chebfun(@(x) [exp(x) exp(-x) sin(x)]);
    g = [1 0 1];
    f./g;
    pass(26) = false;
catch ME
    pass(26) = strcmp(ME.identifier, ...
        'CHEBFUN:CHEBFUN:rdivide:columnRdivide:divisionByZero');
end

end
