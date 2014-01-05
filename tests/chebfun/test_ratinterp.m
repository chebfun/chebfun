% Test for @chebfun/ratinterp.m.
% Based on the Chebfun v4 test written by Pedro Gonnet, Jul. 2011.

function pass = test_ratinterp(pref) 

if (nargin == 0)
    pref = chebpref();
end


% create the input function
f = @(x) (x.^4 - 3) ./ ((x+0.2) .* (x-2.2));
d = [0, 2];
x = chebfun(@(x) x, d);
cf = (x.^4 - 3) ./ ((x+0.2) .* (x-2.2));
N = 100;
pp = [-.2;2.2];

tol = 1e-10;

% Approximate on different point sets using the anonymous function
[p, q, r, mu, nu, poles] = ratinterp(d, f, 10, 10, [], 'type0');
pass(1) = (mu == 4) & (nu == 2) & (max(abs(sort(poles)-pp)) < tol);

[p, q, r, mu, nu] = ratinterp(d, f, 10, 10, [], 'type1');
pass(2) = (mu == 4) & (nu == 2) & (max(abs(sort(roots(q,'all'))-pp)) < tol);

[p, q, r, mu, nu] = ratinterp(d, f, 10, 10, [], 'type2');
pass(3) = (mu == 4) & (nu == 2) & (max(abs(sort(roots(q,'all'))-pp)) < tol);

[p, q, r, mu, nu] = ratinterp(d, f, 10, 10, [], 'equi');
pass(4) = (mu == 4) & (nu == 2) & (max(abs(sort(roots(q,'all'))-pp)) < tol);


% Approximate on different point sets using the chebfun
% Not using type0 since bary does not extend well to the unit circle.
[p, q, r, mu, nu] = ratinterp(cf, 10, 10, [], 'type1');
pass(5) = (mu == 4) & (nu == 2) & (max(abs(sort(roots(q,'all'))-pp)) < tol);

[p, q, r, mu, nu] = ratinterp(cf, 10, 10, [], 'type2');
pass(6) = (mu == 4) & (nu == 2) & (max(abs(sort(roots(q,'all'))-pp)) < tol);

[p, q, r, mu, nu] = ratinterp(cf, 10, 10, [], 'equi');
pass(7) = (mu == 4) & (nu == 2) & (max(abs(sort(roots(q,'all'))-pp)) < tol);


% Approximate on different point sets using a vector of function values
ff = f(1+exp(2i*pi*(0:N-1)'/N));
[p, q, r, mu, nu, poles] = ratinterp(d, ff, 10, 10, N, 'type0');
pass(8) = (mu == 4) & (nu == 2) & (max(abs(sort(poles)-pp)) < tol);

[p, q, r, mu, nu] = ratinterp(d, f(chebpts(N,d,1)), 10, 10, N, 'type1');
pass(9) = (mu == 4) & (nu == 2) & (max(abs(sort(roots(q,'all'))-pp)) < tol);

[p, q, r, mu, nu] = ratinterp(d, f(chebpts(N,d,2)), 10, 10, N, 'type2');
pass(10) = (mu == 4) & (nu == 2) & (max(abs(sort(roots(q,'all'))-pp)) < tol);

[p, q, r, mu, nu] = ratinterp(d, f(linspace(0,2,N)), 10, 10, N, 'equi');
pass(11) = (mu == 4) & (nu == 2) & (max(abs(sort(roots(q,'all'))-pp)) < tol);

% Some further tests:
x = chebfun(@(x) x, [1, 3]);
f = abs(exp(x) - 5);
[p,q,r] = ratinterp(f, 2, 3, [], chebpts(6, [1, 3], 2));
pass(12) = ( (length(p) == 3) & (length(q) == 4) );
pass(13) = norm(f - p./q,inf) < 0.6;
xx = linspace(1, 3, 300);
pass(14) = max(abs((f(xx) - r(xx)))) < 0.6;
[p,q,r] = ratinterp(f, 2, 3, [], 'type2');
pass(15) = norm(f - p./q,inf) < 0.6;
pass(16) = max(abs((f(xx) - r(xx)))) < 0.6;

end