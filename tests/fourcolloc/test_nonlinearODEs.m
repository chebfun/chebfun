% Test file for nonlinear ODEs using the FOURCOLLOC class.

function pass = test_systemODEs(pref)

% Get the preferences.
if ( nargin < 1 )
    pref = cheboppref;
end

%% FIRST ORDER: 
%  u' - sin(u) = cos(x), on [-pi pi].

% Set up domain, rhs f, and intial guesses u0 and v0.
dom = [-pi pi];
f = chebfun(@(x) cos(x), dom);
u0 = @(x) cos(x);
v0 = u0;

% Define the non-linear operator.
N = chebop(@(u) diff(u) - sin(u), dom);

% Solve with FOURIER technology.
u0 = chebfun(u0, dom, 'periodic');
N.init = u0;
pref.discretization = 'periodic';
u = solvebvp(N, f, pref);

% Solve with CHEBYSHEV technology.
v0 = chebfun(v0, dom);
N.init = v0;
N.bc = @(x, v) v(dom(1)) - v(dom(2));
pref.discretization = 'colloc2';
v = solvebvp(N, f, pref);

% Comparison.
tol = pref.errTol;
xx = linspace(dom(1), dom(2), 100);
pass(1) = norm(u(xx) - v(xx), inf) < tol;

%% FIRST ORDER: 
%  u' - u*cos(u) = cos(x), on [-2*pi 2*pi].

% Set up domain, rhs f, and intial guesses u0 and v0.
dom = [-2*pi 2*pi];
f = chebfun(@(x) cos(x), dom);
u0 = @(x) cos(x);
v0 = u0;

% Define the non-linear operator.
N = chebop(@(u) diff(u) - u.*cos(u), dom);

% Solve with FOURIER technology.
u0 = chebfun(u0, dom, 'periodic');
N.init = u0;
pref.discretization = 'periodic';
u = solvebvp(N, f, pref);

% Solve with CHEBYSHEV technology.
v0 = chebfun(v0, dom);
N.init = v0;
N.bc = @(x, v) v(dom(1)) - v(dom(2));
pref.discretization = 'colloc2';
v = solvebvp(N, f, pref);

% Comparison.
tol = pref.errTol;
xx = linspace(dom(1), dom(2), 100);
pass(2) = norm(u(xx) - v(xx), inf) < tol;

%% SECOND ORDER: 
%  u'' - sin(u) = cos(2x), on [-pi pi].

% Set up domain, rhs f, and intial guesses u0 and v0.
dom = [-pi pi];
f = chebfun(@(x) cos(2*x), dom);
u0 = @(x) cos(2*x);
v0 = u0;

% Define the non-linear operator.
N = chebop(@(u) diff(u, 2) - sin(u), dom);

% Solve with FOURIER technology.
u0 = chebfun(u0, dom, 'periodic');
N.init = u0;
pref.discretization = 'periodic';
u = solvebvp(N, f, pref);

% Solve with CHEBYSHEV technology.
v0 = chebfun(v0, dom);
N.init = v0;
N.bc = @(x, v) [ v(dom(1)) - v(dom(2)) ; ...
             feval(diff(v), dom(1)) - feval(diff(v), dom(2)) ];
pref.discretization = 'colloc2';
v = solvebvp(N, f, pref);

% Comparison.
tol = pref.errTol;
xx = linspace(dom(1), dom(2), 100);
pass(3) = norm(u(xx) - v(xx), inf) < tol;

%% SECOND ORDER: 
%  u'' - u^2*cos(u) = cos(x), on [-pi pi].

% Set up domain, rhs f, and intial guesses u0 and v0.
dom = [-pi pi];
f = chebfun(@(x) cos(x), dom);
u0 = @(x) cos(x);
v0 = u0;

% Define the non-linear operator.
N = chebop(@(u) diff(u, 2) - u.^2.*cos(u), dom);

% Solve with FOURIER technology.
u0 = chebfun(u0, dom, 'periodic');
N.init = u0;
pref.discretization = 'periodic';
u = solvebvp(N, f, pref);

% Solve with CHEBYSHEV technology.
v0 = chebfun(v0, dom);
N.init = v0;
N.bc = @(x, v) [ v(dom(1)) - v(dom(2)) ; ...
             feval(diff(v), dom(1)) - feval(diff(v), dom(2)) ];
pref.discretization = 'colloc2';
v = solvebvp(N, f, pref);

% Comparison.
tol = pref.errTol;
xx = linspace(dom(1), dom(2), 100);
pass(4) = norm(u(xx) - v(xx), inf) < tol;

end