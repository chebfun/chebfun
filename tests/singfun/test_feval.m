% Test file for singfun/feval.m

function pass = test_feval(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(786);
x = -1 + 2*rand(100, 1);

%% 
% Check feval on empty set of points
data.exponents = [0, 0];
data.singType = {'none', 'none'};
f = singfun(@(x) x, data, pref);
pass(1) = isempty(feval(f, []));
%%
% Check feval on a SINGFUN without exponents
fh = @(x) sin(cos(10*x.^2));
data.exponents = [0, 0];
data.singType = {'none', 'none'};
f = singfun(fh, data, pref);
pass(2) = norm(feval(f,x) - feval(fh,x), inf) < get(f, 'epslevel');

%%
% Check feval on a SINGFUN with negative exponents
a = 1 + rand();
b = 1 + rand();
fh = @(x) sin(cos(10*x.^2))./((1+x).^a.*(1-x).^b);
data.exponents = [-a, -b];
data.singType = {'sing', 'sing'};
f = singfun(fh, data, pref);
err = norm(feval(f,x) - feval(fh,x), inf);
tol = 1e3*get(f, 'epslevel');
pass(3) = err < tol;

%%
% Check feval on a SINGFUN with positive exponents
a = rand();
b = rand();
fh = @(x) sin(cos(10*x.^2)).*(1+x).^a.*(1-x).^b;
data.exponents = [a, b];
data.singType = {'root', 'root'};
f = singfun(fh, data, pref);
pass(4) = norm(feval(f,x) - feval(fh,x), inf) < get(f, 'epslevel');

end
