function pass = test_paramODE(pref)
% Test solving a parameter dependent ODE. 
% Nick Hale, August 2011

% NOTE: Taken from V4 test chebop_paramODE.

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 1e-10;

%% Simple problem

x = chebfun('x');

% Natural setup
N = chebop(@(x, u, a) x.*u + .001*diff(u,2) + a);
N.lbc = @(u, a) [u + a + 1 ; diff(u)];
N.rbc = @(u, a) u - 1;
u = mldivide(N, chebfun(0), pref);
res1 = N(x, u);

% Forced setup using a system
N = chebop(@(x, u, a) [x.*u + .001*diff(u,2) + a ; diff(a)]);
N.lbc = @(u, a) [u + a + 1 ; diff(u)];
N.rbc = @(u, a) u - 1;
v = mldivide(N, [chebfun(0) ; chebfun(0)], pref);
res2 = N(x, v);

err(1) = norm(res1) + norm(res2{1}) + norm(u{1}-v{1}) + norm(u{2}-v{2});

%% More complicated (piecewise system + 2 params)
x = chebfun('x');
N = chebop(@(x,u,a,b,v) [x.*v + .001*diff(u,2) + a + 2*b ; diff(v)-u], [-1 0 1]);
N.lbc = @(u, a, b, v) [u + a ; diff(u)];
N.rbc = @(u, a, b, v) [u ; diff(u)-a ; v-b];
rhs = [sin(x) ; 1];
u = mldivide(N, rhs, pref);

res = N(x, u) - rhs;
res1 = res{1};
res2 = res{2};
err(2) = norm(res1, inf) + norm(res2, inf);

%%

pass = err < tol;

end

