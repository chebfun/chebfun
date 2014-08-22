function pass = test_linearInit(pref)
% A linear CHEBOP test. This test tests that if we pass an initial guess to a
% linear CHEBOP N that the solution is actually correct.

%% Setup
dom = [0 pi];
if ( nargin == 0 )
    pref = cheboppref;
end
pref.errTol = 1e-11;

%% Simple scalar problem
N = chebop(@(x,u) diff(u,2) + x.*u, dom);
N.lbc = 2; 
N.rbc = 3;

x = chebfun(@(x) x, dom);
rhs = sin(x);

%% Simple
% Start with chebcolloc2
pref.discretization = @chebcolloc2;
N.init = sin(20*x);
u1 = solvebvp(N, rhs, pref);
err(1) = norm(N(u1) - rhs);
err(2) = abs(feval(N.lbc(u1), dom(1))) + abs(feval(N.rbc(u1), dom(end)));
%% Add a breakpoint, change discretization
pref.discretization = @ultraS;
N.domain = [0 1 pi];
u2 = solvebvp(N, rhs, pref);
err(3) = norm(N(u2) - rhs);
err(4) = abs(feval(N.lbc(u1), dom(1))) + abs(feval(N.rbc(u1), dom(end)));
%% System:
d = [-pi pi];
A = chebop(@(x,u,v) [u - diff(v) ; diff(u) + v],d);
A.lbc = @(u,v) u + 1;
A.rbc = @(u,v) v;
x = chebfun('x',d);
N.init = [cos(20*x); sin(20*x)];
rhs = [exp(-x.^2); 1+exp(-x.^2)];
uv = mldivide(A, rhs, pref);
err(5) = norm(A(uv)-rhs);
err(6) = abs(feval(A.lbc(uv{1}, uv{2}), d(1))) + ...
    abs(feval(A.rbc(uv{1}, uv{2}), d(end)));
%% Happy?
tol = 100*pref.errTol;
pass = err < tol;
