function pass = test_periodic_nonlin(pref)
% Test 'periodic' syntax for nonlinear ODEs.

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 1e-8;

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
pref.discretization = @chebcolloc2;
u = solvebvp(N, f, pref);

% Solve imposing directly the periodic boundary condition.
N.bc = @(x, u) [u(dom(2)) - u(dom(1)); ...
   feval(diff(u), dom(2)) - feval(diff(u), dom(1))];
N.init = u0;
v = N \ f;

% Compare.
pass(1) = norm(u - v, inf) < tol;

%% Test the TRIGCOLLOC class. FIRST ORDER: 
%  u' - u*cos(u) = cos(x), on [-pi pi].

% Set up domain, rhs f, and intial guesses u0 and v0.
dom = [-pi pi];
f = chebfun(@(x) cos(x), dom);
u0 = f;

% Define the non-linear operator.
N = chebop(@(u) diff(u) - u.*cos(u), dom);
N.bc = 'periodic';
N.init = u0;

% Solve with TRIGCOLLOC.
u = N \ f;

pass(2) = norm(N*u - f) < tol;
pass(3) = abs(u(dom(1)) - u(dom(2))) < tol;
pass(4) = isequal(get(u.funs{1}, 'tech'), @trigtech);

%% Test the TRIGCOLLOC class. SECOND ORDER: 
%  u'' - u^2*cos(u) = cos(x), on [-pi pi].

% Set up domain, rhs f, and intial guesses u0 and v0.
dom = [-pi pi];
f = chebfun(@(x) cos(x), dom);
u0 = f;

% Define the non-linear operator.
N = chebop(@(u) diff(u, 2) - u.^2.*cos(u), dom);
N.bc = 'periodic';
N.init = u0;

% Solve with TRIGCOLLOC.
u = N \ f;

pass(5) = norm(N*u - f) < tol;
pass(6) = abs(u(dom(1)) - u(dom(2))) < tol;
pass(7) = abs(feval(diff(u), dom(1)) - feval(diff(u), dom(2))) < tol;
pass(8) = isequal(get(u.funs{1}, 'tech'), @trigtech);

%% Test the TRIGSPEC class. SECOND ORDER: 
%  u'' - u^2*cos(u) = cos(x), on [-pi pi].

% Set up domain, rhs f, and intial guesses u0 and v0.
dom = [-pi pi];
f = chebfun(@(x) cos(x), dom);
u0 = f;

% Define the non-linear operator.
N = chebop(@(u) diff(u, 2) - u.^2.*cos(u), dom);
N.bc = 'periodic';
N.init = u0;

% Solve with TRIGSPEC.
pref.discretization = @trigspec;
u = solvebvp(N, f, pref);

pass(9) = norm(N*u - f) < 1e3*tol;
pass(10) = abs(u(dom(1)) - u(dom(2))) < tol;
pass(11) = abs(feval(diff(u), dom(1)) - feval(diff(u), dom(2))) < tol;
pass(12) = isequal(get(u.funs{1}, 'tech'), @trigtech);
pass(13) = isreal(u);

end
