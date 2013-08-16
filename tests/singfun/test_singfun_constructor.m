% Test file for SINGFUN constructor.

function pass = test_singfun_constructor(pref)

% Get preferences:
if ( nargin < 1 )
    pref = singfun.pref.singfun;
end
% Set the tolerance:
tol = 1e3*pref.eps;

pass = zeros(1, 4); % Pre-allocate pass matrix

%%
% Test calling syntax
a = 1;
b = .5;
fh = @(x) sin(x)./((1+x).^a.*(1-x).^b);
f = singfun(fh, [-a, -b]);
g = singfun(fh, [-a, -b], {'pole', 'sing'});
pass(1) = isequal(f,g);

%%
% Select some random points as sample points
% These random points are in [-.999, .999]
seedRNG(777)
%x = -.999 + 1.998*rand(1, 100);
x = -1 + 2*rand(1, 100);

fh = @(x) sin(x)./((1-x).*(1+x));
f = singfun(fh, [], {'pole', 'pole'});
pass(2) = norm(feval(f, x) - feval(fh, x), inf) < tol
norm(feval(f, x) - feval(fh, x), inf)

%%
%
fh = @(x) sin(exp(cos(10*pi*x)))./((1+x).^a.*(1-x).^b);
f = singfun(fh, [-a, -b]);
pass(3) = norm(feval(f, x) - feval(fh, x), inf) < tol
norm(feval(f, x) - feval(fh, x), inf)


%%