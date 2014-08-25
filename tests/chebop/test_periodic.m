function pass = test_periodic(pref)
% Test 'periodic' syntax for CHEBOP

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 1e-8;

%% A simple ODE:
N = chebop(@(u) diff(u, 2) - u);
% Periodic RHS:
rhs = chebfun(@(x) sin(pi*x));

% Apply periodic BCs manually:
N.bc = @(u) [u(-1) - u(1) ; feval(diff(u),-1) - feval(diff(u), 1)];
u = mldivide(N, rhs, pref);

% Apply via 'periodic':
N.bc = 'periodic';
v = mldivide(N, rhs, pref);

% Compare results:
err(1) = norm(u - v);

%% A periodic piecewise system:
d = [-pi, 0, pi];
A = chebop(d);
A.op = @(x, u, v) [u-diff(v) ; diff(u,2)+v];
x = chebfun('x',d);
f = [chebfun(0, d) ; cos(x)];
A.bc = 'periodic';
uv = mldivide(A, f, pref);

trueSoln = [cos(x+3*pi/4)/sqrt(2) ; cos(x+pi/4)/sqrt(2)];
err(2) = norm(uv - trueSoln);

%% Eigenvalue problem:

d = [-pi, 0, pi];
A = chebop(d);
A.op = @(x, u, v) [u-diff(v) ; diff(u,2)+v];
A.bc = 'periodic';
B = chebop(d);
B.op = @(x, u, v) [v + u ; diff(v)];

[V, D] = eigs(A, B, 5, 0, pref);
e = diag(D);

% Sort the eigenvalues to ensure things will work on all machines. Eigs() does
% not appear to return the eigenvalues in a consistent way for different
% machines. We expect five pair of eigenvalues to appear, check whether they are
% all there
% `sort(x)` sorts complex values by abs() and then by angle(). In order have
% consistent sorting, we move all the eigenvalues up into the first quadrant
% before sorting them.
[ignored, idx] = sort(real(e));
e = e(idx);

e12 = e(1:2);
e35 = e(3:5);
[ignored, idx] = sort(imag(e12));
e12 = e12(idx);
[ignored, idx] = sort(imag(e35));
e35 = e35(idx);

e = [e12; e35];

err(3) = norm(real(e) - [0 0 1 1 1].', inf) + ...
    norm(imag(e) - [-1 1 -1 0 1].', inf);
err(4) = norm(V{1}(pi) - V{1}(pi), inf) + norm(V{2}(pi) - V{2}(pi), inf);

%% Nonlinear ODE.
% u'' - sin(u) = cos(2x).

% Domain.
dom = [-pi pi];

% Define the rhs, and the intial guess.
f = chebfun(@(x) cos(2*x), dom);
u0 = f;

% Define the non-linear operator.
N = chebop(@(u) diff(u, 2) - sin(u), dom);

% Solve using the periodic tag.
N.bc = 'periodic';
N.init = u0;
u = N \ f;

% Solve imposing directly the periodic boundary condition.
N.bc = @(x, u) [u(dom(2)) - u(dom(1)); ...
   feval(diff(u), dom(2)) - feval(diff(u), dom(1))];
N.init = u0;
v = N \ f;

% Compare.
err(5) = norm(u - v, inf);

%% System of nonlinear ODEs.
%  u - v' + v = 0, u'' - cos(v) = cos(x).

% Define the domain.
dom = [-pi pi];

% Define the rhs, and the intial guesses.
f = [ chebfun(0, dom) ; chebfun(@(x) cos(x), dom) ];
u0 = f;

% Define the non-linear operator.
N = chebop(@(x, u, v) [ u - diff(v) + v ; diff(u, 2) - cos(v) ], dom);

% Solve using the periodic tag.
N.bc = 'periodic';
N.init = u0;
u = N \ f;

% Solve imposing directly the periodic boundary condition.
N.bc = @(x, u, v) [ u(dom(1)) - u(dom(2)); ...
                    feval(diff(u), dom(1)) - feval(diff(u), dom(2)) ; ...
                    v(dom(1)) - v(dom(2)) ];
N.init = u0;
v = N \ f;

% Compare.
err(6) = norm(u - v, inf);

%%

pass = err < tol;

end
