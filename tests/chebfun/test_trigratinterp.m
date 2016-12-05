% Test file for trigratinterp
function pass = test_trigratinterp(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebtech.techPref();
end

%%
% Generate a few random points to use as test values.
seedRNG(6178);
tol = 1e-14;

m = 1;
n = 1;
N = (m+n);
th = sort(-1 + 2*rand(2*N+1,1));
fh = @(t) tan(pi*t);
fk = fh(th);

[p, q, r] = trigratinterp(fk, m, n, [], th, 0);
pass(1) = length(p) == 2*m+1;
pass(2) = length(q) == 2*n+1;
pass(3) = norm(p(th)./q(th) - fh(th), inf) < tol;
tt = 2 * rand(100, 1) - 1;
pass(4) = norm((p(tt)./q(tt)-fh(tt)), inf) < 1e3*tol;

%% Check least squres
N = 101;
[p, q, r] = trigratinterp(fh, m, n, N, [], 0);
pass(5) = length(p) == 2*m+1;
pass(6) = length(q) == 2*n+1;
pass(7) = norm(p(th)./q(th) - fh(th), inf) < tol;
tt = 2 * rand(100, 1) - 1;
pass(8) = norm((p(tt)./q(tt)-fh(tt)), inf) < 1e3*tol;


%%
m = 10;
n = 15;
N = m+n;
t = chebfun('t');
f = sin(pi*t)./cos(pi*t);
th = trigpts(2*N+1);
fk = f(th);
[p, q, r] = trigratinterp(fk, m, n);
pass(9) = length(p) == 3;
pass(10) = length(q) == 3;
th = trigpts(2*N+1);
pass(11) = norm(p(th)./q(th) - fh(th), inf) < 1e3*tol;
tt = 2 * rand(100, 1) - 1;
pass(12) = norm((p(tt)./q(tt)-fh(tt)), inf) < 1e3*tol;

%%
m = 5;
n = 5;
N = m+n;
fh = @(t) exp(-4*sin(pi*t/2).^2);
f = chebfun(fh, 'trig');
th = -1+2*rand(2*N+1, 1);
[p, q, r] = trigratinterp(f, m, n, [], th, 0);
pass(13) = length(p) == 2*m+1;
pass(14) = length(q) == 2*n+1;
th = trigpts(2*N+1);
pass(15) = norm(p(th)./q(th) - fh(th), inf) < 5e-10;

%%
m = 10;
n = 15;
N = m+n;
fh = @(t) cos(pi*t)./cos(2*pi*t);
th = trigpts(2*N+1);
fk = fh(th);
[p, q, r] = trigratinterp(fk, m, n);
pass(16) = length(p) == 3;
pass(6) = length(q) == 5;
th = trigpts(2*N+1);
pass(17) = norm(p(th)./q(th) - fh(th), inf) < 1e2*tol;
tt = 2 * rand(100, 1) - 1;
pass(18) = norm((p(tt)./q(tt)-fh(tt)), inf) < 1e3*tol;

fh = @(x) exp(sin(x));
a = 2; b = 2 + 2*pi;
f = chebfun(fh, [a, b], 'trig');
[p, q, r] = trigratinterp(f, 6, 6, [], 2 + 2*pi/(51)*(0:50));  
xx = linspace(a, b, 10001);
pass(19) = norm(fh(xx) - r(xx), inf) < 10 * tol;

% test a simple discrete case
fi = [2, 3, 1];
xi = [-1, 0.5, 0.8];
[p, q, r ] = trigratinterp(fi, 1, 0, 3, xi);
pass(20) = norm(p(xi)./q(xi) - fi, inf) < 10*tol;
pass(21) = norm(r(xi) - fi, inf) < 10*tol;
pass(22) = length(p) == 3;
pass(23) = length(q) == 1;


end
