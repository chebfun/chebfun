function pass = test_gmres(pref)
% Test the Chebfun implementation of GMRES for solving Lu = f, where L is an
% operator, and both f and u are chebfuns
%  Sheehan Olver

if ( nargin == 0 )
    pref = chebfunpref();
end

tol = pref.chebfuneps;

d = [-1 1];
x = chebfun('x', d, pref);
f = chebfun('exp(x)', d, pref);
w = 100;
L = @(u) diff(u) + 1i*w*u;

[u, flag] = gmres(L, f);
pass(1) = ~flag;
err = abs(sum(f.*exp(1i.*w.*x)) - (u(1).*exp(1i.*w) - u(-1).*exp(-1i.*w)));
pass(2) = err < 100*tol;

end
