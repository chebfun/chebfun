function pass = test_periodic_nonlin(pref)
% Test 'periodic' syntax for nonlinear ODEs.

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 1e-10;

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
pref.discretization = 'chebcolloc2';
u = solvebvp(N, f, pref);

% Solve imposing directly the periodic boundary condition.
N.bc = @(x, u) [u(dom(2)) - u(dom(1)); ...
   feval(diff(u), dom(2)) - feval(diff(u), dom(1))];
N.init = u0;
v = N \ f;

% Compare.
err(1) = norm(u - v, inf);

%% Test the FOURCOLLOC class. FIRST ORDER: 
%  u' - u*cos(u) = cos(x), on [-2*pi 2*pi].

% Set up domain, rhs f, and intial guesses u0 and v0.
dom = [-2*pi 2*pi];
f = chebfun(@(x) cos(x), dom);
u0 = f;

% Define the non-linear operator.
N = chebop(@(u) diff(u) - u.*cos(u), dom);
N.bc = 'periodic';
N.init = u0;

% Solve with FOURIER technology.
u = N \ f;

% Solve with CHEBYSHEV technology.
pref.discretization = 'chebcolloc2';
v = solvebvp(N, f, pref);

% Comparison.
tol = pref.errTol;
xx = linspace(dom(1), dom(2), 100);
err(2) = norm(u(xx) - v(xx), inf);

%% Test the FOURCOLLOC class. SECOND ORDER: 
%  u'' - u^2*cos(u) = cos(x), on [-pi pi].

% Set up domain, rhs f, and intial guesses u0 and v0.
dom = [-pi pi];
f = chebfun(@(x) cos(x), dom);
u0 = f;

% Define the non-linear operator.
N = chebop(@(u) diff(u, 2) - u.^2.*cos(u), dom);
N.bc = 'periodic';
N.init = u0;

% Solve with FOURIER technology.
u = N \ f;

% Solve with CHEBYSHEV technology.
pref.discretization = 'chebcolloc2';
v = solvebvp(N, f, pref);

% Comparison.
tol = pref.errTol;
xx = linspace(dom(1), dom(2), 100);
err(3) = norm(u(xx) - v(xx), inf);

%%
pass = err < tol;

end
