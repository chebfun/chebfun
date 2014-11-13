function pass = test_paramODE(pref)
% Test solving a parameter dependent ODE. 
% Nick Hale, August 2011

% NOTE: Taken from V4 test chebop_paramODE.

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 10*pref.errTol;

%% Simple problem

x = chebfun('x');

% Natural setup
N = chebop(@(x, u, a) x.*u + .001*diff(u,2) + a);
% N = chebop(@(x, u, a) x.*u + .001*diff(u,2) + a.*diff(u));
N.lbc = @(u, a) [u + a + 1 ; diff(u)];
N.rbc = @(u, a) u - 1;
% N.init = [chebfun(0) ; 0];
u = mldivide(N, chebfun(0), pref);
res1 = N(x, u);

%% Linear, simple problem, forced setup using a system

N = chebop(@(x, u, a) [x.*u + .001*diff(u,2) + a ; diff(a)]);
% N = chebop(@(x, u, a) [x.*u + .001*diff(u,2) + a.*diff(u) ; diff(a)]);
N.lbc = @(u, a) [u + a + 1 ; diff(u)];
N.rbc = @(u, a) u - 1;
v = mldivide(N, [chebfun(0) ; chebfun(0)], pref);
res2 = N(x, v);
err(1) = norm(res1) + norm(res2{1}) + norm(u{1}-v{1}) + norm(u{2}-v{2});

%% Linear: More complicated (piecewise system + 2 params)

x = chebfun('x');
N = chebop(@(x,u,v,a,b) [x.*v + .001*diff(u,2) + a + 2*b ; a + diff(v) - u], [-1 1]);
N.lbc = @(u, v, a, b) [u + a ; diff(u)];
N.rbc = @(u, v, a, b) [u ; diff(u)-a ; v-b];
% N.init = [chebfun(0) ; 0 ; 0 ; chebfun(0)];
rhs = [sin(x) ; x];
u = mldivide(N, rhs, pref);

res = N(x, u) - rhs;
res1 = res{1};
res2 = res{2};
err(2) = norm(res1, inf) + norm(res2, inf);

%%

pass = err < tol;

end
