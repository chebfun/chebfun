function pass = test_periodic_system(pref)
% Test 'periodic' syntax for system of ODEs.

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 1e-10;

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
pref.discretization = 'chebcolloc2';
u = solvebvp(N, f, pref);

% Solve imposing directly the periodic boundary condition.
N.bc = @(x, u, v) [ u(dom(1)) - u(dom(2)); ...
                    feval(diff(u), dom(1)) - feval(diff(u), dom(2)) ; ...
                    v(dom(1)) - v(dom(2)) ];
N.init = u0;
v = N \ f;

% Compare.
err(1) = norm(u - v, inf);

%% Test the FOURCOLLOC class. FIRST AND SECOND ORDER LINEAR ODEs.
% u - v' = 0, u'' + v = cos(x), on [-pi pi].

% Set domain, operator L, and rhs f.
dom = [-pi, pi];
L = chebop(@(x, u, v) [ u - diff(v) ; diff(u, 2) + v ], dom);
L.bc = 'periodic';
F = [ chebfun(0, dom) ; chebfun(@(x) cos(x), dom) ];

% Solve with FOURIER technology.
U = L \ F;

% Solve with CHEBYSHEV technology.
pref.discretization = 'chebcolloc2';
V = solvebvp(L, F, pref);

% Comparison.
xx = linspace(dom(1), dom(2), 100);
tol = pref.errTol;
err(2) = norm(U{1}(xx) - V{1}(xx), inf);
err(3) = norm(U{2}(xx) - V{2}(xx), inf);

%% Test the FOURCOLLOC class. FIRST AND SECOND ORDER NONLINEAR ODEs.
%  u - v' + v = 0, u'' - cos(v) = cos(x).

% Set domain, operator L, and rhs f.
dom = [-pi, pi];
N = chebop(@(x, u, v) [ u - diff(v) + v ; diff(u, 2) - cos(v) ], dom);
N.bc = 'periodic';
F = [ chebfun(0, dom) ; chebfun(@(x) cos(x), dom) ];
N.init = F;

% Solve with FOURIER technology.
U = N \ F;

% Solve with CHEBYSHEV technology.
pref.discretization = 'chebcolloc2';
V = solvebvp(N, F, pref);

% Comparison.
xx = linspace(dom(1), dom(2), 100);
tol = pref.errTol;
err(4) = norm(U{1}(xx) - V{1}(xx), inf);
err(5) = norm(U{2}(xx) - V{2}(xx), inf);

%%
pass = err < tol;

end
