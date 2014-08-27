% Test file for linear ODEs using the FOURCOLLOC class.

function pass = test_systemODEs(pref)

% Get the preferences.
if ( nargin < 1 )
    pref = cheboppref;
end

%% FIRST ORDER AND CONSTANT COEFFICIENTS: 
%  u' + u = cos(x), on [0 2*pi].

% Set domain, operator L, and rhs f.
dom = [0 2*pi];
L = chebop(@(u) diff(u) + u, dom);
f = chebfun(@(x) cos(x), dom);

% Solve with FOURIER technology.
L.bc = 'periodic';
u = L \ f;

% Compare with exact solution.
tol = pref.errTol;
exact = chebfun(@(x) 1/2*cos(x) + 1/2*sin(x), dom, 'periodic');
pass(1) = norm(u - exact, inf) < tol;

%% FIRST ORDER AND VARIABLES COEFFICIENTS: 
%  u' + (1+cos(x))u = cos(2x), on [-2*pi 2*pi].

% Set domain, variable coefficient bb, and rhs f.
dom = [-2*pi 2*pi];
bb = @(x) 1 + cos(x);
f = chebfun(@(x) cos(2*x), dom);

% Solve with FOURIER technology.
b = chebfun(bb, dom, 'periodic');
L = chebop(@(u) diff(u) + b.*u, dom);
L.bc = 'periodic';
u = L \ f;

% Solve with CHEBYSHEV technology.
b = chebfun(bb, dom);
L = chebop(@(v) diff(v) + b.*v, dom);
L.bc = 'periodic';
pref.discretization = 'chebcolloc2';
v = solvebvp(L, f, pref);

% Comparison.
tol = pref.errTol;
xx = linspace(dom(1), dom(2), 100);
pass(2) = norm(u(xx) - v(xx), inf) < tol;

%% SECOND ORDER AND CONSTANT COEFFICIENTS: 
%  u'' + 10u' + 5u = cos(x), on [-2*pi 2*pi].

% Set domain, constants coefficients a and b, operator L,
% and rhs f.
dom = [-2*pi 2*pi];
a = 10;
b = 5;
L = chebop(@(u) diff(u, 2) + a*diff(u) + b*u, dom); 
L.bc = 'periodic';
f = chebfun(@(x) cos(x), dom);

% Solve with FOURIER technology.
u = L \ f;

% Compare with exact solution.
tol = pref.errTol;
exact = chebfun(@(x) 1/29*cos(x) + 5/58*sin(x), dom, 'periodic');
pass(3) = norm(u - exact, inf) < tol;

%% SECOND ORDER AND VARIABLE COEFFICIENTS: 
%  (2+cos(4x))u'' + sin(cos(2x))u' + exp(cos(x))u = cos(x), on [-pi pi].

% Set domain, variable coefficients aa, bb and cc, and rhs f.
dom = [-pi pi];
aa = @(x) 2 + cos(4*x);
bb = @(x) sin(cos(2*x));
cc = @(x) exp(cos(x));
f = chebfun(@(x) cos(x), dom);

% Solve with FOURIER technology.
a = chebfun(aa, dom, 'periodic');
b = chebfun(bb, dom, 'periodic');
c = chebfun(cc, dom, 'periodic');
L = chebop(@(u) a.*diff(u, 2) + b.*diff(u) + c.*u, dom);
L.bc = 'periodic';
u = L \ f;

% Solve with CHEBYSHEV technology.
a = chebfun(aa, dom);
b = chebfun(bb, dom);
c = chebfun(cc, dom);
L = chebop(@(v) a.*diff(v, 2) + b.*diff(v) + c.*v, dom);
L.bc = 'periodic';
pref.discretization = 'chebcolloc2';
v = solvebvp(L, f, pref);

% Comparison.
tol = 10*pref.errTol;
xx = linspace(dom(1), dom(2), 100);
pass(4) = norm(u(xx) - v(xx), inf) < tol;

%% THIRD ORDER AND VARIABLE COEFFICIENTS: 
%  (2+cos(x))u''' + sin(cos(2x))u'' + exp(cos(x))u' + sin(x)u = cos(x),
%  on [-pi pi].

% Set domain, variable coefficients aa, bb, cc and dd, and rhs f.
dom = [-pi pi];
aa = @(x) 2 + cos(x);
bb = @(x) sin(cos(2*x));
cc = @(x) exp(cos(x));
dd = @(x) sin(x);
f = chebfun(@(x) cos(x), dom);

% Solve with FOURIER technology.
a = chebfun(aa, dom, 'periodic');
b = chebfun(bb, dom, 'periodic');
c = chebfun(cc, dom, 'periodic');
d = chebfun(dd, dom, 'periodic');
L = chebop(@(u) a.*diff(u, 3) + b.*diff(u, 2) + c.*diff(u) + d.*u, dom);
L.bc = 'periodic';
u = L \ f;

% Solve with CHEBYSHEV technology.
a = chebfun(aa, dom);
b = chebfun(bb, dom);
c = chebfun(cc, dom);
d = chebfun(dd, dom);
L = chebop(@(v) a.*diff(v, 3) + b.*diff(v, 2) + c.*diff(v) + d.*v, dom);
L.bc = 'periodic';
pref.discretization = 'chebcolloc2';
v = solvebvp(L, f, pref);

% Comparison.
tol = 100*pref.errTol;
xx = linspace(dom(1), dom(2), 100);
pass(5) = norm(u(xx) - v(xx), inf) < tol;

%% FOURTH ORDER AND VARIABLE COEFFICIENTS: 
%  (2+cos(x))u'''' + sin(cos(2x))u''' + exp(cos(x))u'' + ... 
%  sin(x)u' + sin(2*x)u = cos(10*x), on [-pi pi].

% Set domain, variable coefficients aa, bb, cc, dd and ee, and rhs f.
dom = [-pi pi];
aa = @(x) 2 + cos(x);
bb = @(x) sin(cos(2*x));
cc = @(x) exp(cos(x));
dd = @(x) sin(x);
ee = @(x) sin(2*x);
f = chebfun(@(x) cos(10*x), dom);

% Solve with FOURIER technology.
a = chebfun(aa, dom, 'periodic');
b = chebfun(bb, dom, 'periodic');
c = chebfun(cc, dom, 'periodic');
d = chebfun(dd, dom, 'periodic');
e = chebfun(ee, dom, 'periodic');
L = chebop(@(u) a.*diff(u, 4) + b.*diff(u, 3) + c.*diff(u, 2) + ...
    d.*diff(u) + e.*u, dom);
L.bc = 'periodic';
u = L \ f;

% Solve with CHEBYSHEV technology.
a = chebfun(aa, dom);
b = chebfun(bb, dom);
c = chebfun(cc, dom);
d = chebfun(dd, dom);
e = chebfun(ee, dom);
L = chebop(@(v) a.*diff(v, 4) + b.*diff(v, 3) + c.*diff(v, 2) + ...
    d.*diff(v) + e.*v, dom);
L.bc = 'periodic';
pref.discretization = 'chebcolloc2';
v = solvebvp(L, f, pref);

% Comparison.
tol = 100*pref.errTol;
xx = linspace(dom(1), dom(2), 100);
pass(6) = norm(u(xx) - v(xx), inf) < tol;

end
