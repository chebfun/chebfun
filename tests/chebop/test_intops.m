function pass = test_intops(pref)
% Test integral operators
% Toby Driscoll 28 May 2009
% Nick Hale 6 Jan 2012

% NOTE: Taken from V4 test chebop_intops.

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 1e-10;

% From http://en.wikipedia.org/wiki/Integro-differential_equation
d = [0 5];
x = chebfun('x', d);
N = chebop(@(u) diff(u) + 2*u + 5*cumsum(u), d);
N.bc = @(x,u) u(0);
u = mldivide(N, 1, pref);
u_exact = 0.5*exp(-x).*sin(2*x);
err(1) = norm(u - u_exact);

% Fredholm
d = [0, 1];
x = chebfun(@(x) x, d);
K = @(x,y) sin(2*pi*(x-y));
A = chebop(@(x,u) u + fred(K,u), d);
u = x.*exp(x);
f = A*u;
err(2) = norm(u - mldivide(A, f, pref));

% Volterra
d = [0, pi];
x = chebfun(@(x) x, d);
K = @(x,y) x.*y;
A = chebop(@(x,u) u - volt(K,u), d);
f = x.^2.*cos(x) + (1-x).*sin(x);
u = mldivide(A, f, pref);
err(3) = norm( u - sin(x) );
err(4) = norm( A*u - f );

pass = err < tol;

%%
% From #2122
K = @(x,t) exp(t.*(x-t));
N = chebop(@(x,y) diff(y) + y - x.*(1+2*x).*volt(K, y) - 1 - 2*x);
N.lbc = 1;
y = N\0;
pass(5) = norm(N.op(y)) + abs(y(0)-1);

%%
% From #2179
K = @(x,y) exp(-(x-y));
A = chebop(@(x,u) diff(u) + fred(K,u));
A.lbc = 0; 
u = A\1;
pass(6) = norm(A*u - 1, inf) < tol;

%%
% From https://groups.google.com/forum/#!topic/chebfun-users/0dpsggp1RyA . 
% Also, see #2210.

L = chebop(@(u) diff(u)+sum(u), [0 1]);
L.lbc = 1; 
u = L\0;
pass(7) = norm(u - (1- 2*chebfun('x',[0 1])/3), inf) < tol;

%%
% From #2185 
A = chebop(@(x,u) diff(u,2) + sin(2*pi*x)*u + sum(u).*u);
A.bc = 'periodic'; 
u = A\1;
pass(8) = norm(A*u-1) < tol;

f = chebfun(@(x) sin(5*pi*x));
K = @(x,y) cos(-pi*(x-y));
A = chebop(@(x,u) 1/25*diff(u,2) + fred(K, u));
A.bc = 'periodic'; 
u = A\f;
pass(9) = norm(A*u-f) < tol;


%%
% From #2360.
K = @(x,y) 1+0*x;
N = chebop(@(x,u) - u + exp(x) + 1/3*x*(1-exp(3*x)) + x*volt(K, u^3), [0, 1]);
u = N\0;
pass(10) = norm(u-chebfun(@exp, [0 1])) < tol;

end

