% Test file for singfun/plus.m

function pass = test_feval(pref)

% Get preferences.
if ( nargin < 1 )
    pref = singfun.pref.singfun;
end

% Generate a few random points to use as test values.
seedRNG(786);
x = -1 + 2*rand(100, 1);

pass = zeros(1, 4); % Pre-allocate pass vector
tol = 100*pref.eps;

%% 
% Check feval on empty set of points
f = singfun(@(x) x, [0, 0], {'none', 'none'});
pass(1) = isempty(feval(f, []));
%%
% Check feval on a SINGFUN without exponents
fh = @(x) sin(cos(10*x.^2));
f = singfun(fh, [0, 0], {'none', 'none'});
pass(2) = norm(feval(f,x) - feval(fh,x), inf) < tol;

%%
% Check feval on a SINGFUN with negative exponents
a = 1 + rand();
b = 1 + rand();
fh = @(x) sin(cos(10*x.^2))./((1+x).^a.*(1-x).^b);
f = singfun(fh, [-a, -b], {'sing', 'sing'});
pass(3) = norm(feval(f,x) - feval(fh,x), inf) < tol;

%%
% Check feval on a SINGFUN with positive exponents
a = rand();
b = rand();
fh = @(x) sin(cos(10*x.^2)).*(1+x).^a.*(1-x).^b;
f = singfun(fh, [a, b], {'branch', 'branch'});
pass(4) = norm(feval(f,x) - feval(fh,x), inf) < tol;

end