% Test file for chebtech2 constructor.

function pass = test_constructor(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebtech.techPref();
end
% Set the tolerance:
tol = 100*pref.eps;

%%
% Test on a scalar-valued function:
pref.extrapolate = 0;
pref.refinementFunction = 'nested';
f = @(x) sin(x);
g = populate(chebtech2, f, [], [], pref);
x = chebtech2.chebpts(length(g.values));
pass(1) = norm(f(x) - g.values, inf) < tol;
pass(2) = g.vscale == sin(1) && g.ishappy && g.epslevel < tol;

pref.extrapolate = 1;
pref.refinementFunction = 'nested';
g = populate(chebtech2, f, [], [], pref);
x = chebtech2.chebpts(length(g.values));
pass(3) = norm(f(x) - g.values, inf) < tol;
pass(4) = norm(g.vscale - sin(1), inf) < tol && logical(g.epslevel);

pref.extrapolate = 0;
pref.refinementFunction = 'resampling';
g = populate(chebtech2, f, [], [], pref);
x = chebtech2.chebpts(length(g.values));
pass(5) = norm(f(x) - g.values, inf) < tol;
pass(6) = g.vscale == sin(1) && logical(g.epslevel);

pref.extrapolate = 1;
pref.refinementFunction = 'resampling';
g = populate(chebtech2, f, [], [], pref);
x = chebtech2.chebpts(length(g.values));
pass(7) = norm(f(x) - g.values, inf) < tol;
pass(8) = norm(g.vscale - sin(1), inf) < tol && logical(g.epslevel);

%%
% Test on an array-valued function:
pref.extrapolate = 0;
pref.refinementFunction = 'nested';
f = @(x) [sin(x) cos(x) exp(x)];
g = populate(chebtech2, f, [], [], pref);
x = chebtech2.chebpts(length(g.values));
pass(9) = norm(f(x) - g.values, inf) < tol;
pass(10) = norm(g.vscale - [sin(1) cos(0) exp(1)], inf) < 10*g.epslevel ...
    && logical(g.epslevel);

pref.extrapolate = 1;
pref.refinementFunction = 'nested';
g = populate(chebtech2, f, [], [], pref);
x = chebtech2.chebpts(length(g.values));
pass(11) = norm(f(x) - g.values, inf) < tol;
pass(12) = norm(g.vscale - [sin(1) cos(0) exp(1)], inf) < tol && logical(g.epslevel);

pref.extrapolate = 0;
pref.refinementFunction = 'resampling';
g = populate(chebtech2, f, [], [], pref);
x = chebtech2.chebpts(length(g.values));
pass(13) = norm(f(x) - g.values, inf) < tol;
pass(14) = norm(g.vscale - [sin(1) cos(0) exp(1)], inf) < 10*g.epslevel ...
    && logical(g.epslevel);

pref.extrapolate = 1;
pref.refinementFunction = 'resampling';
g = populate(chebtech2, f, [], [], pref);
x = chebtech2.chebpts(length(g.values));
pass(15) = norm(f(x) - g.values, inf) < tol;
pass(16) = norm(g.vscale - [sin(1) cos(0) exp(1)], inf) < tol && logical(g.epslevel);

%%
% Some other tests:

% This should fail with an error:
try
    f = @(x) x + NaN;
    populate(chebtech2, f);
    pass(17) = false;
catch ME
    pass(17) = strcmp(ME.message, 'Too many NaNs/Infs to handle.');
end

% As should this:
try
    f = @(x) x + Inf;
    populate(chebtech2, f);
    pass(18) = false;
catch ME
    pass(18) = strcmp(ME.message, 'Too many NaNs/Infs to handle.');
end

% Test that the extrapolation option avoids endpoint evaluations.
pref.extrapolate = 1;
try 
    populate(chebtech2, @(x) [F(x) F(x)], [], [], pref);
    pass(19) = true;
catch ME %#ok<NASGU>
    pass(19) = false;
end

    function y = F(x)
        if ( any(abs(x) == 1) )
            error('Extrapolate should prevent endpoint evaluation.');
        end
        y = sin(x);
    end

end
