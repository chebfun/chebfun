function pass = test_intops(pref)
% Test integral operators
% Toby Driscoll 28 May 2009
% Nick Hale 6 Jan 2012

% NOTE: Taken from V4 test chebop_intops.

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 1e-10;

% Fredholm
d = domain(0, 1); 
x = chebfun(@(x) x, d);
F = fred(@(x, y) sin(2*pi*(x-y)), d);
A = eye(d) + F;
u = x.*exp(x);
f = A*u;
res = u-mldivide(A, f, pref);
pass(1) = norm(res{1}) < 1e6*tol;

% Volterra
d = domain(0,pi);
x = chebfun(@(x) x, d);
V = volt(@(x, y) x.*y, d);
f = x.^2.*cos(x) + (1-x).*sin(x);
u = mldivide(1-V, f, pref);
res1 = u - sin(x);
pass(2) = norm( res1{1}) < 1e6*tol;
res2 = (1-V)*u - f;
pass(3) = norm( res2{1} ) < 1e4*tol;


%% Now available as chebops!

% Fredholm
d = [0, 1];
x = chebfun(@(x) x, d);
K = @(x,y) sin(2*pi*(x-y));
A = chebop(@(x,u) u + fred(K,u), d);
u = x.*exp(x);
f = A*u;
pass(4) = norm(u-A\f) < 1e6*tol;

% Volterra
d = [0, pi];
x = chebfun(@(x) x, d);
K = @(x,y) x.*y;
A = chebop(@(x,u) u - volt(K,u), d);
f = x.^2.*cos(x) + (1-x).*sin(x);
u = mldivide(A, f, pref);

pass(5) = norm( u - sin(x) ) < 1e6*tol;
pass(6) = norm( A*u - f ) < 1e4*tol;

end

