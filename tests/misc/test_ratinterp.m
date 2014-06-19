% Test for ratinterp.m.
% Based on the Chebfun v4 test written by Pedro Gonnet, Jul. 2011.

function pass = test_ratinterp(pref)

% Create the input function.
f = @(x) (x.^4 - 3)./((x + 0.2).*(x - 2.2));
dom = [0 2];
x = chebfun(@(x) x, dom);
cf = f(x);
N = 100;

% Approximate on different point sets using the anonymous function.
[p, q, r, mu, nu, poles, res] = ratinterp(dom, f, 10, 10, [], 'type0');
pass(1) = (mu == 4) && (nu == 2) && ...
    (max(abs(sort(poles) - [-0.2 ; 2.2])) < 1e-10);

[p, q, r, mu, nu] = ratinterp(dom, f, 10, 10, [], 'type1');
pass(2) = (mu == 4) && (nu == 2) && ...
    (max(abs(sort(roots(q, 'all')) - [-0.2 ; 2.2])) < 1e-10);

[p, q, r, mu, nu] = ratinterp(dom, f, 10, 10, [], 'type2');
pass(3) = (mu == 4) && (nu == 2) && ...
    (max(abs(sort(roots(q, 'all')) - [-0.2 ; 2.2])) < 1e-10);

[p, q, r, mu, nu ] = ratinterp(dom, f, 10, 10, [], 'equi');
pass(4) = (mu == 4) && (nu == 2) && ...
    (max(abs(sort(roots(q, 'all')) - [-0.2 ; 2.2])) < 1e-10);

% Approximate on different point sets using the chebfun.
% (Don't type0, since bary does not extend well to the unit circle.)
[p, q, r, mu, nu] = ratinterp(cf, 10, 10, [], 'type1');
pass(5) = (mu == 4) && (nu == 2) && ...
    (max(abs(sort(roots(q, 'all')) - [-0.2 ; 2.2])) < 1e-10);

[p, q, r, mu, nu] = ratinterp(cf, 10, 10, [], 'type2');
pass(6) = (mu == 4) && (nu == 2) && ...
    (max(abs(sort(roots(q, 'all')) - [-0.2 ; 2.2])) < 1e-10);

[p, q, r, mu, nu] = ratinterp(cf, 10, 10, [], 'equi');
pass(7) = (mu == 4) && (nu == 2) && ...
    (max(abs(sort(roots(q, 'all')) - [-0.2 ; 2.2])) < 1e-10);

% Approximate on different point sets using a vector of function values.
[p, q, r, mu, nu, poles] = ratinterp(dom, f(1 + exp(2i*pi*(0:N-1)'/N)), ...
    10, 10, N, 'type0');
pass(8) = (mu == 4) && (nu == 2) && ...
    (max(abs(sort(poles) - [-0.2 ; 2.2])) < 1e-10);

[p, q, r, mu, nu] = ratinterp(dom, f(chebpts(N, dom, 1)), 10, 10, N, 'type1');
pass(9) = (mu == 4) && (nu == 2) && ...
    (max(abs(sort(roots(q, 'all')) - [-0.2 ; 2.2])) < 1e-10);

[p, q, r, mu, nu] = ratinterp(dom, f(chebpts(N, dom, 2)), 10, 10, N, 'type2');
pass(10) = (mu == 4) && (nu == 2) && ...
    (max(abs(sort(roots(q, 'all')) - [-0.2 ; 2.2])) < 1e-10);

[p, q, r, mu, nu] = ratinterp(dom, f(linspace(0, 2, N)), 10, 10, N, 'equi');
pass(11) = (mu == 4) && (nu == 2) && ...
    (max(abs(sort(roots(q, 'all')) - [-0.2 ; 2.2])) < 1e-10);

% Make sure that the example from the documentation works.
[p, q, r, mu, nu] = ratinterp([-1 1], @(x) 1./(x - 0.2), 10, 10, [], 'type2');
pass(12) = (mu == 0) && (nu == 1) && (abs(roots(q) - 0.2) < 1e-10);

% Some further tests:
x = chebfun(@(x) x, [1, 3]);
f = abs(exp(x) - 5);

[p, q, r] = ratinterp(f, 2, 3, [], chebpts(6, [1, 3], 2));
pass(13) = (length(p) == 3) && (length(q) == 4);
pass(14) = norm(f - p./q,inf) < 0.6;
xx = linspace(1, 3, 300);
pass(15) = max(abs((f(xx) - r(xx)))) < 0.6;

[p, q, r] = ratinterp(f, 2, 3, [], 'type2');
pass(16) = norm(f - p./q,inf) < 0.6;
pass(17) = max(abs((f(xx) - r(xx)))) < 0.6;

% Robustness will get rid of the denominator in this scenario.
f = @(x) 1./(x - 1.7);
[p, q, r, mu, nu, pol] = ratinterp(f, 128, 1, [], 'type0', 1e-14);
pass(18) = (nu == 0);

% Check interpolation in arbitrary complex nodes.  (Use roots of unity but
% don't inform ratinterp().)
f = @(x) sin(x)./(x - 0.1);
zk = exp(2*pi*1i*(0:1:31).'/32);
[p, q, r, mu, nu, pol] = ratinterp(f, 30, 1, [], zk, 1e-14);
pass(19) = abs(pol - 0.1) < 1e-10;

% Check that ratinterp accepts a domain input: (#841)
[p1,q1] = ratinterp([-1,1], @exp, 1, 1);
warnState = warning('off', 'CHEBFUN:ratinterp:domainDeprecated');
[p2,q2] = ratinterp(domain(-1,1), @exp, 1, 1);
warning(warnState);
pass(20) = norm(p1-p2) == 0 && norm(q1-q2) == 0;

%% Test the examples from the m-file:
f = @(x) 1./(x - 0.2);
xx = linspace(-1,1,20);
[p, q, r] = ratinterp([-1 1], f, 10, 10, [], 'type2', 0);
err = norm(r(xx) - f(xx));
pass(21) = err < 1e-10;
[p, q, r] = ratinterp([-1 1], f, 10, 10, [], 'type2', 0);
err = norm(r(xx) - f(xx));
pass(22) = err < 1e-10;

end
