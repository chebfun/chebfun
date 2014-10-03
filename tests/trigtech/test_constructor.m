% Test file for trigtech constructor.

function pass = test_constructor(pref)

% Get preferences:
if ( nargin < 1 )
    pref = trigtech.techPref();
end

%%
% Test on a scalar-valued function:
pref.refinementFunction = 'nested';
f = @(x) tanh(sin(pi*x));
g = populate(trigtech, f, [], [], pref);
x = trigtech.trigpts(length(g.coeffs));
pass(1) = norm(f(x) - g.values, inf) < 10*g.vscale.*g.epslevel;

% Test on an array-valued function:
pref.refinementFunction = 'nested';
f = @(x) [exp(sin(pi*x)) sin(cos(4*pi*x)) cos(pi*x)];
g = populate(trigtech, f, [], [], pref);
x = trigtech.trigpts(length(g.coeffs));
pass(2) = norm(f(x) - g.values, inf) < 10*max(g.vscale.*g.epslevel);

%%
% Test on a scalar-valued function:
pref.refinementFunction = 'resampling';
f = @(x) tanh(sin(pi*x));
g = populate(trigtech, f, [], [], pref);
x = trigtech.trigpts(length(g.coeffs));
pass(3) = norm(f(x) - g.values, inf) < 10*g.vscale.*g.epslevel;

% Test on an array-valued function:
pref.refinementFunction = 'resampling';
f = @(x) [exp(sin(pi*x)) sin(cos(4*pi*x)) cos(pi*x)];
g = populate(trigtech, f, [], [], pref);
x = trigtech.trigpts(length(g.coeffs));
pass(4) = norm(f(x) - g.values, inf) < 10*max(g.vscale.*g.epslevel);

%%
% Some other tests:

% This should fail with an error:
try
    f = @(x) sin(pi*x) + NaN;
    populate(trigtech, f, [], [], pref);
    pass(5) = false;
catch ME
    pass(5) = strcmp(ME.message, 'Cannot handle functions that evaluate to Inf or NaN.');
end

% As should this:
try
    f = @(x) sin(pi*x) + Inf;
    populate(trigtech, f, [], [], pref);
    pass(6) = false;
catch ME
    pass(6) = strcmp(ME.message, 'Cannot handle functions that evaluate to Inf or NaN.');
end

% Check that things don't crash if pref.minSamples and pref.maxLength are equal.
try
    pref.minSamples = 8;
    pref.maxLength = 8;
    populate(trigtech, @(x) sin(pi*x), [], [], pref);
    pass(7) = true;
catch
    pass(7) = false;
end

% Test logical-valued functions:
f = trigtech(@(x) x > -2);
g = trigtech(1);
pass(8) = normest(f - g) < f.epslevel;

f = trigtech(@(x) x < -2);
g = trigtech(0);
pass(9) = normest(f - g) < f.epslevel;

end
