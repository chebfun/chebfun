% Test file for chebtech1 constructor.

function pass = test_constructor(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebtech.techPref();
end

%%
% Test on a scalar-valued function:
pref.refinementFunction = 'nested';
f = @(x) sin(x);
g = populate(chebtech1, f, [], [], pref);
x = chebtech1.chebpts(length(g.values));
pass(1) = norm(f(x) - g.values, inf) < 10*g.vscale.*g.epslevel;

% Test on an array-valued function:
pref.refinementFunction = 'nested';
f = @(x) [sin(x) cos(x) exp(x)];
g = populate(chebtech1, f, [], [], pref);
x = chebtech1.chebpts(length(g.values));
pass(2) = norm(f(x) - g.values, inf) < 10*max(g.vscale.*g.epslevel);

%%
% Test on a scalar-valued function:
pref.refinementFunction = 'resampling';
f = @(x) sin(x);
g = populate(chebtech1, f, [], [], pref);
x = chebtech1.chebpts(length(g.values));
pass(3) = norm(f(x) - g.values, inf) < 10*g.vscale.*g.epslevel;

% Test on an array-valued function:
pref.refinementFunction = 'resampling';
f = @(x) [sin(x) cos(x) exp(x)];
g = populate(chebtech1, f, [], [], pref);
x = chebtech1.chebpts(length(g.values));
pass(4) = norm(f(x) - g.values, inf) < 10*max(g.vscale.*g.epslevel);

%%
% Some other tests:

% This should fail with an error:
try
    f = @(x) x + NaN;
    populate(chebtech1, f, [], [], pref);
    pass(5) = false;
catch ME
    pass(5) = strcmp(ME.message, 'Too many NaNs/Infs to handle.');
end

% As should this:
try
    f = @(x) x + Inf;
    populate(chebtech1, f, [], [], pref);
    pass(6) = false;
catch ME
    pass(6) = strcmp(ME.message, 'Too many NaNs/Infs to handle.');
end

% Check that things don't crash if pref.minPoints and pref.maxPoints are equal.
try
    pref.minPoints = 8;
    pref.maxPoints = 8;
    populate(chebtech1, @sin, [], [], pref);
    pass(7) = true;
catch
    pass(7) = false;
end