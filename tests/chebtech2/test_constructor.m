% Test file for chebtech2 constructor.
% Here, we check populate().  (This function is not user-facing.)

function pass = test_constructor(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebtech.techPref();
end
% Set the tolerance:
tol = 100*pref.eps;

% Initialize with default data:
data = chebtech.parseDataInputs(struct());

%%
% Test on a scalar-valued function:
pref.extrapolate = 0;
pref.refinementFunction = 'nested';
f = @(x) sin(x);
g = populate(chebtech2, f, data, pref);
x = chebtech2.chebpts(length(g.coeffs));
values = g.coeffs2vals(g.coeffs);
pass(1) = norm(f(x) - values, inf) < tol;
pass(2) = abs(vscale(g) - sin(1)) < eps && g.ishappy && eps < tol;

pref.extrapolate = 1;
pref.refinementFunction = 'nested';
g = populate(chebtech2, f, data, pref);
x = chebtech2.chebpts(length(g.coeffs));
values = g.coeffs2vals(g.coeffs);
pass(3) = norm(f(x) - values, inf) < tol;
pass(4) = norm(vscale(g) - sin(1), inf) < tol && logical(eps);

pref.extrapolate = 0;
pref.refinementFunction = 'resampling';
g = populate(chebtech2, f, data, pref);
x = chebtech2.chebpts(length(g.coeffs));
values = g.coeffs2vals(g.coeffs);
pass(5) = norm(f(x) - values, inf) < tol;
pass(6) = abs(vscale(g) - sin(1)) < eps && logical(eps);

pref.extrapolate = 1;
pref.refinementFunction = 'resampling';
g = populate(chebtech2, f, data, pref);
x = chebtech2.chebpts(length(g.coeffs));
values = g.coeffs2vals(g.coeffs);
pass(7) = norm(f(x) - values, inf) < tol;
pass(8) = norm(vscale(g) - sin(1), inf) < tol && logical(eps);

%%
% Test on an array-valued function:
pref.extrapolate = 0;
pref.refinementFunction = 'nested';
f = @(x) [sin(x) cos(x) exp(x)];
g = populate(chebtech2, f, data, pref);
x = chebtech2.chebpts(length(g.coeffs));
values = g.coeffs2vals(g.coeffs);
pass(9) = norm(f(x) - values, inf) < tol;

pref.extrapolate = 1;
pref.refinementFunction = 'nested';
g = populate(chebtech2, f, data, pref);
x = chebtech2.chebpts(length(g.coeffs));
values = g.coeffs2vals(g.coeffs);
pass(10) = norm(f(x) - values, inf) < tol;

pref.extrapolate = 0;
pref.refinementFunction = 'resampling';
g = populate(chebtech2, f, data, pref);
x = chebtech2.chebpts(length(g.coeffs));
values = g.coeffs2vals(g.coeffs);
pass(11) = norm(f(x) - values, inf) < tol;

pref.extrapolate = 1;
pref.refinementFunction = 'resampling';
g = populate(chebtech2, f, data, pref);
x = chebtech2.chebpts(length(g.coeffs));
values = g.coeffs2vals(g.coeffs);
pass(12) = norm(f(x) - values, inf) < tol;

%%
% Some other tests:

% This should fail with an error:
try
    f = @(x) x + NaN;
    populate(chebtech2, f, data, pref);
    pass(13) = false;
catch ME
    pass(13) = strcmp(ME.message, 'Too many NaNs/Infs to handle.');
end

% As should this:
try
    f = @(x) x + Inf;
    populate(chebtech2, f, data, pref);
    pass(14) = false;
catch ME
    pass(14) = strcmp(ME.message, 'Too many NaNs/Infs to handle.');
end

% Test that the extrapolation option avoids endpoint evaluations.
pref.extrapolate = 1;
try
    populate(chebtech2, @(x) [F(x) F(x)], data, pref);
    pass(15) = true;
catch ME %#ok<NASGU>
    pass(15) = false;
end
pref.extrapolate = 0;

    function y = F(x)
        if ( any(abs(x) == 1) )
            error('Extrapolate should prevent endpoint evaluation.');
        end
        y = sin(x);
    end

% Check that things don't crash if pref.minSamples and pref.maxLength are equal.
try
    pref.minSamples = 8;
    pref.maxLength = 8;
    populate(chebtech2, @sin, data, pref);
    pass(16) = true;
catch
    pass(16) = false;
end

%%
% Test logical-valued functions:
f = chebtech2(@(x) x > -2);
g = chebtech2(1);
pass(17) = normest(f - g) < eps;

f = chebtech2(@(x) x < -2);
g = chebtech2(0);
pass(18) = normest(f - g) < eps;

end
