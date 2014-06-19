function pass = test_sign(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% Initialise random vector:
seedRNG(6178);
x = 2 * rand(100, 1) - 1;
hvsde = @(x) .5*(sign(x) + 1);

%% Simple tests
pref.splitting = 0;
f = chebfun('x.^2', pref);   
tol = get(f, 'epslevel')*get(f, 'hscale');
f1 = sign(f);
pass(1,1) = all(feval(f1, x) == 1);

% Test if pointValues are dealt with correctly: 
f = chebfun(@(x) cos(pi*x) + 2, sort(x)', pref);
f.pointValues = -pi*ones(length(f.ends), 1);
f1 = sign(f);
pass(1,2) = all(f1.pointValues == -1);

% Test also on longer intervals:
f = chebfun(@(x) x, [-1, 100], pref);
f1 = sign(f);
h = chebfun({-1, 1}, [-1, 0, 100]);
pass(1,3) = normest(h - f1) < tol;

% Test functions with breakpoints:
f = chebfun(@(x) x, [-1, 1e-10, 1], pref);
f1 = sign(f);
h = restrict(h,[-1, 1]);
pass(1,4) = normest(f1 - h) < tol;

% Real, imaginary and complex CHEBFUN objects:
f = chebfun({-1, 1, -1, 1, -1, 1}, -3:3);

%% Real CHEBFUN
gHandle1 = @(x) sin(pi*x);
g = chebfun(@(x) gHandle1(x), -3:3, pref);
g1 = sign(g);
pass(2,1) = length(g1.funs) == 6;
pass(2,2) = normest(f - g1) < tol;
pref.splitting = 1;
h1 = chebfun(@(x) sign(gHandle1(x)), -3:3, pref);
pass(2,3) = length(h1.funs) == 6;
pass(2,4) = normest(f - h1) < 100*tol;

%% Array-valued CHEBFUN
f1 = chebfun(@(x) feval(f, [x, x]) , -3:3, pref);
gHandle2 = @(x) cos(pi*(x-.5));
g = chebfun(@(x) [gHandle1(x), gHandle2(x)] , -3:3, pref);
g4 = sign(g);
pass(3,1) = length(g4.funs) == 6;
pass(3,2) = normest(f1 - g4) < tol;
h4 = chebfun(@(x) sign([gHandle1(x), gHandle2(x)]), -3:3, pref);
pass(3,3) = length(h4.funs) == 6;
pass(3,4) = normest(f1 - h4) < tol;

%% A more complicated function:
f = chebfun(@(x) sin(1i*x).*(1i*x + exp(5i*x)));
g = chebfun(@(x) sign(sin(1i*x).*(1i*x + exp(5i*x))),[-1 0 1], ...
    'extrapolate', 'on');
h = sign(f);
pass(4,:) = normest(g - h) < 200*get(h, 'epslevel')*length(h);

%% Test sign() for a complex-valued CHEBFUN.
f = chebfun(@(x) exp(2*pi*1i*x)./(1 + (x - 0.1).^2), [-1 1]);
h = sign(f);
h_exact = @(x) exp(2*pi*1i*x);
pass(5,:) = norm(feval(h, x) - h_exact(x), inf) < 10*vscale(h)*epslevel(h);

%% Test on singular function: a real case

% Set a domain:
dom = [-2 7];

% Generate a few random points to use as test values:
x = diff(dom) * rand(100, 1) + dom(1);

pow = -0.5;
op = @(x) (dom(2)-x).^pow.*( sin(10*x).^2 );
f = chebfun(op, dom, 'exps', [0 pow], 'splitting', 'on');
s = sign(f);
pass(6,:) = ( norm(feval(s-1, x), inf) < eps );

%% Test on singular function: a complex case

% Set a domain:
dom = [-1 1];

% Generate a few random points to use as test values:
x = diff(dom) * rand(100, 1) + dom(1);

pow = -0.5;
op = @(x) (dom(2)-x).^pow.*exp(2*pi*1i*x);
f = chebfun(op, dom, 'exps', [0 pow], 'splitting', 'on');
s = sign(f);
s_exact = @(x) exp(2*pi*1i*x);
vals_s = feval(s, x);
vals_exact = feval(s_exact, x);
err = vals_s - vals_exact;
pass(7,:) = ( norm(err, inf) < 1e1*epslevel(s).*get(s, 'vscale') );

%% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];
domCheck = [-1e2 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) (1-exp(-x.^2))./x;
f = chebfun(op, dom);
s = sign(f);
sVals = feval(s, x);
op = @(x) 2*hvsde(x) - 1;
sExact = op(x);
err = sVals - sExact;
pass(8,:) = all( ~norm(err, inf) );

end
