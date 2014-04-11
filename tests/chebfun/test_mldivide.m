function pass = test_mldivide(pref)

if ( nargin == 0 )
    pref = chebpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
xr = 2 * rand(100, 1) - 1;

T = restrict(chebpoly(0:3), [-1 -0.5 0 0.5 1]);
L = restrict(legpoly(0:3), [-1 0 1]);

%% Scalar A:
pass(1) = normest(2\T - .5*T) < 10*epslevel(T);

%% Numeric A:
B = T.';
A = (1:4).';
x = A\B;
x0 = feval(x, 0);
x0_true = -1/15;
pass(2) = abs(x0 - x0_true) < 10*epslevel(x);

A = eye(4);
X = A\B;
X0 = feval(X, 0);
pass(3) = norm(X0 - [1 0 -1 0].') < 10*epslevel(X);

%% Transposed A:
A = T.';
B = (1:4).';
X = A\B;
X0 = feval(X, 0);
X0_true = -2.625;
pass(4) = abs(X0 - X0_true) < 10*epslevel(X);

%% Else

X = T\L;
C = diag([1 1 .75 .625]); 
C(1, 3) = .25;
C(2, 4) = .375;
pass(5) = norm(X - C, inf) < 10*epslevel(L);

%% Test on SINGFUN:

% scalar * [1 x INF] = [1 x INF] => scalar\row SINGFUN:
f = chebfun(@(x)sin(20*x)./(x+1), 'exps', [-1 0], 'splitting', 'on');
f = f.';
g = 3\f;
op = @(x) sin(20*x)./(3*(x+1));
g_vals = feval(g, xr);
g_exact = op(xr).';
err = g_vals - g_exact;
pass(6) = norm(err, inf) < 1e4*vscale(g)*epslevel(g);

% [1 x INF] * [INF x 1] = scalar => row SINGFUN\scalar:
f = chebfun(@(x)(sin(100*x).^2+1)./((x+1).^0.25), 'exps', [-0.25 0], 'splitting', 'on');
f = f.';
g = f\3;
err = f*g - 3;
pass(7) = abs(err) < vscale(g)*epslevel(g);

% [INF x 1] * SCALAR = [INF x 1] => column SINGFUN\column SINGFUN:

f = chebfun(@(x)3*(x.^2+3)./(x+1).^0.4, 'exps', [-0.4 0], 'splitting', 'on');
g = chebfun(@(x)(x.^2+3)./(x+1).^0.4, 'exps', [-0.4 0], 'splitting', 'on');
h = f\g;
err = h - 1/3;
pass(8) = norm(err, inf) < vscale(f)*epslevel(f);

% [TODO]: Revive the following test.

%% Test for function defined on unbounded domain:

% % Set the domain:
% dom = [-Inf 3*pi];
% domCheck = [-1e6 3*pi];
% 
% % Generate a few random points to use as test values:
% x = diff(domCheck) * rand(100, 1) + domCheck(1);
% 
% % A*X = B, where both A and B are UNBNDFUNs ==> X = A\B
% 
% opA = @(x) [exp(x) x.*exp(x) (1-exp(x))./x];
% A = chebfun(opA, dom);
% opB = @(x) [(2*x+1).*exp(x) exp(x) 2*(1-exp(x))./x];
% B = chebfun(opB, dom);
% X = A\B;
% res = A*X - B;
% err = feval(res, x);
% pass(9) = norm(err(:), inf) < max([get(A,'epslevel').*get(A,'vscale') ...
%     get(B,'epslevel').*get(B,'vscale')]);

end