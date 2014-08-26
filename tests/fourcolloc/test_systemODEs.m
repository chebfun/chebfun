% Test file for system of ODEs using the FOURCOLLOC class.

function pass = test_systemODEs(pref)

% Get the preferences.
if ( nargin < 1 )
    pref = cheboppref;
end

%% FIRST AND SECOND ORDER LINEAR ODEs.
% u - v' = 0, u'' + v = cos(x), on [-pi pi].

% Set domain, operator L, and rhs f.
dom = [-pi, pi];
L = chebop(@(x, u, v) [ u - diff(v) ; diff(u, 2) + v ], dom);
F = [ chebfun(0, dom) ; chebfun(@(x) cos(x), dom) ];

% Solve with FOURIER technology.
pref.discretization = 'periodic';
U = solvebvp(L, F, pref);

% Solve with CHEBYSHEV technology.
pref.discretization = 'chebcolloc2';
L.bc = 'periodic';
V = solvebvp(L, F, pref);

% Comparison.
xx = linspace(dom(1), dom(2), 100);
tol = pref.errTol;
pass(1) = norm(U{1}(xx) - V{1}(xx), inf) < tol && norm(U{2}(xx) - V{2}(xx), inf);

%% FIRST AND SECOND ORDER NONLINEAR ODEs.
%  u - v' + v = 0, u'' - cos(v) = cos(x).

% Set domain, operator L, and rhs f.
dom = [-pi, pi];
N = chebop(@(x, u, v) [ u - diff(v) + v ; diff(u, 2) - cos(v) ], dom);
F = [ chebfun(0, dom) ; chebfun(@(x) cos(x), dom) ];

% Solve with FOURIER technology.
pref.discretization = 'periodic';
N.init = [ chebfun(0, dom, 'periodic') ; chebfun(@(x) cos(x), dom, 'periodic') ];
U = solvebvp(N, F, pref);

% Solve with CHEBYSHEV technology.
pref.discretization = 'chebcolloc2';
N.bc = 'periodic';
N.init = F;
V = solvebvp(N, F, pref);

% Comparison.
xx = linspace(dom(1), dom(2), 100);
tol = pref.errTol;
pass(2) = norm(U{1}(xx) - V{1}(xx), inf) < tol && norm(U{2}(xx) - V{2}(xx), inf);

end
