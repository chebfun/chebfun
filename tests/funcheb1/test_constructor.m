% Test file for funcheb1 constructor.

function pass = test_constructor(pref)

% Get preferences:
if ( nargin < 1 )
    pref = funcheb.pref;
end
% Set the tolerance:
tol = 100*pref.funcheb.eps;

%%
% Test on a scalar-valued function:
pref.funcheb.refinementFunction = 'default';
f = @(x) sin(x);
g = populate(funcheb1, f, [], [], pref);
x = funcheb1.chebpts(length(g.values));
pass(1) = norm(f(x) - g.values, inf) < tol;

%%
% Test on a vector-valued function:
pref.funcheb.refinementFunction = 'default';
f = @(x) [sin(x) cos(x) exp(x)];
g = populate(funcheb1, f, [], [], pref);
x = funcheb1.chebpts(length(g.values));
pass(2) = norm(f(x) - g.values, inf) < tol;

%%
% Some other tests:

% This should fail with an error:
try
    f = @(x) x + NaN;
    populate(funcheb1, f, [], [], pref);
    pass(3) = false;
catch ME
    pass(3) = strcmp(ME.message, 'Too many NaNs to handle.');
end

% As should this:
try
    f = @(x) x + Inf;
    populate(funcheb1, f, [], [], pref);
    pass(4) = false;
catch ME
    pass(4) = strcmp(ME.message, 'Too many NaNs to handle.');
end
