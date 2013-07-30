function pass = test_besselj(pref)

% Grab some preferences
if ( nargin == 0 )
    pref = chebfun.pref();
end

tol = 1000*pref.chebfun.eps;

% Choose a domain:
a = -1;
b = 5;
% Make a linear chebfun:
x = chebfun(@(x) x, [a b], pref);
xx = linspace(a, b, 100);
c = 1 + 1i;

% Make a test space:
F = @(x) besselj(0, c*x);
f = F(x);
pass(1) = norm(f(xx) - F(xx), inf) < tol;
F = @(x) besselj(0, c*x, 1);
f = F(x);
pass(2) = norm(f(xx) - F(xx), inf) < tol;

F = @(x) besselj(1, c*x);
f = F(x);
pass(3) = norm(f(xx) - F(xx), inf) < tol;
F = @(x) besselj(1, c*x, 1);
f = F(x);
pass(4) = norm(f(xx) - F(xx), inf) < tol;

F = @(x) besselj(.5, c*x+.2);
f = F(x);
pass(5) = norm(f(xx) - F(xx), inf) < tol;
F = @(x) besselj(.5, c*x+.2, 1);
f = F(x);
pass(6) = norm(f(xx) - F(xx), inf) < tol;

% [TODO]: Requires SINGFUN.
% F = @(x) besselj(.5, c*x);
% f = F(x);
% pass(7) = norm(f(xx) - F(xx), inf) < tol;
% F = @(x) besselj(.5, c*x, 1);
% f = F(x);
% pass(8) = norm(f(xx) - F(xx), inf) < tol;

end


