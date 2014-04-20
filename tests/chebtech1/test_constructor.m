% Test file for chebtech1 constructor.

function pass = test_constructor(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebtech.techPref();
end

%%
% Test on a scalar-valued function:
pref.refinementFunction = 'default';
f = @(x) sin(x);
g = populate(chebtech1, f, [], [], pref);
x = chebtech1.chebpts(length(g.coeffs));
values = g.coeffs2vals(g.coeffs);
pass(1) = norm(f(x) - values, inf) < 10*g.vscale.*g.epslevel;

%%
% Test on an array-valued function:
pref.refinementFunction = 'default';
f = @(x) [sin(x) cos(x) exp(x)];
g = populate(chebtech1, f, [], [], pref);
x = chebtech1.chebpts(length(g.coeffs));
values = g.coeffs2vals(g.coeffs);
pass(2) = norm(f(x) - values, inf) < 10*max(g.vscale.*g.epslevel);

%%
% Some other tests:

% This should fail with an error:
try
    f = @(x) x + NaN;
    populate(chebtech1, f, [], [], pref);
    pass(3) = false;
catch ME
    pass(3) = strcmp(ME.message, 'Too many NaNs/Infs to handle.');
end

% As should this:
try
    f = @(x) x + Inf;
    populate(chebtech1, f, [], [], pref);
    pass(4) = false;
catch ME
    pass(4) = strcmp(ME.message, 'Too many NaNs/Infs to handle.');
end

% Check that things don't crash if pref.minPoints and pref.maxPoints are equal.
try
    pref.minPoints = 8;
    pref.maxPoints = 8;
    populate(chebtech1, @sin, [], [], pref);
    pass(5) = true;
catch
    pass(5) = false;
end

