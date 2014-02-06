function pass = test_mrdivide(pref)

if ( nargin == 0 )
    pref = chebpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
xr = 2 * rand(100, 1) - 1;

T = restrict(chebpoly(0:3), [-1 -0.5 0 0.5 1]);
L = restrict(legpoly(0:3), [-1 0 1]);

%% Scalar A:
pass(1) = normest(T/2 - .5*T) < 10*epslevel(T);

%% Numeric A:
B = T;
A = (1:4);
x = A/B;
x0 = feval(x, 0);
x0_true = -2.625;
pass(2) = abs(x0 - x0_true) < 10*epslevel(x);

A = eye(4);
X = A/L;
X0 = feval(X, 0);
pass(3) = norm(X0 - [.5 0 -1.25 0].') < 10*epslevel(X);

%% ~Transposed A:
A = T;
B = (1:4);
x = A/B;
x0 = feval(x, 0);
x0_true = -1/15;
pass(4) = abs(x0 - x0_true) < 10*epslevel(x);

%% Else:
X = T.'/L.';
C = diag([1 1 4/3 8/5]); 
C(3, 1) = -1/3;
C(4,2) = -3/5;
pass(5) = norm(X - C, inf) < 1e2*epslevel(L);

%% Test on SINGFUN:

% [INF x 1] * scalar = [INF x 1] => column SINGFUN/scalar:
f = chebfun(@(x)sin(20*x)./(x+1), 'exps', [-1 0], 'splitting', 'on');
g = f/3;
op = @(x) sin(20*x)./(3*(x+1));
g_vals = feval(g, xr);
g_exact = op(xr);
err = g_vals - g_exact;
pass(6) = norm(err, inf) < 1e4*vscale(g)*epslevel(g);

% [1 x INF] * [INF x 1] = scalar => scalar/column SINGFUN:
f = chebfun(@(x)(sin(100*x).^2+1)./(1-x).^0.4, 'exps', [0 -0.4], 'splitting', 'on');
g = 5/f;
err = g*f - 5;
pass(7) = abs(err) < vscale(g)*epslevel(g);

% SCALAR * [1 x INF] = [1 x INF] => row SINGFUN/row SINGFUN:
f = chebfun(@(x)(sin(30*x).^2+1)./(1-x).^0.3, 'exps', [0 -0.3], 'splitting', 'on');
f = f.';
g = chebfun(@(x)9*(sin(30*x).^2+1)./(1-x).^0.3, 'exps', [0 -0.3], 'splitting', 'on');
g = g.';
h = g/f;
err = h - 9;
pass(8) = abs(err) < 3*vscale(f)*epslevel(f);

end