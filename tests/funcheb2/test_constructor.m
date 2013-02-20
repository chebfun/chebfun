function pass = test_constructor(pref)

% Get preferences:
if ( nargin < 1 )
    pref = funcheb2.pref;
end
% Set the tolerance:
tol = 100*pref.funcheb2.eps;

%%
% Test on a scalar-valued function:
pref.funcheb2.extrapolate = 0;
pref.funcheb2.refinementFunction = 'nested';
f = @(x) sin(x);
[values, coeffs, vscale, ishappy, epslevel] = funcheb2.constructor(f, [], pref);
x = funcheb2.chebpts(length(values));
pass(1) = norm(f(x) - values, inf) < tol;
pass(2) = vscale == sin(1) && ishappy && epslevel < tol;

pref.funcheb2.extrapolate = 1;
pref.funcheb2.refinementFunction = 'nested';
[values, coeffs, vscale, epslevel] = funcheb2.constructor(f, [], pref);
x = funcheb2.chebpts(length(values));
pass(3) = norm(f(x) - values, inf) < tol;
pass(4) = norm(vscale - sin(1), inf) < tol && logical(epslevel);

pref.funcheb2.extrapolate = 0;
pref.funcheb2.refinementFunction = 'resampling';
[values, coeffs, vscale, epslevel] = funcheb2.constructor(f, [], pref);
x = funcheb2.chebpts(length(values));
pass(5) = norm(f(x) - values, inf) < tol;
pass(6) = vscale == sin(1) && logical(epslevel);

pref.funcheb2.extrapolate = 1;
pref.funcheb2.refinementFunction = 'resampling';
[values, coeffs, vscale, epslevel] = funcheb2.constructor(f, [], pref);
x = funcheb2.chebpts(length(values));
pass(7) = norm(f(x) - values, inf) < tol;
pass(8) = norm(vscale - sin(1), inf) < tol && logical(epslevel);

%%
% Test on a vector-valued function:
pref.funcheb2.extrapolate = 0;
pref.funcheb2.refinementFunction = 'nested';
f = @(x) [sin(x) cos(x) exp(x)];
[values, coeffs, vscale, epslevel] = funcheb2.constructor(f, [], pref);
x = funcheb2.chebpts(length(values));
pass(9) = norm(f(x) - values, inf) < tol;
pass(10) = all(vscale == [sin(1) cos(0) exp(1)]) && logical(epslevel);

pref.funcheb2.extrapolate = 1;
pref.funcheb2.refinementFunction = 'nested';
[values, coeffs, vscale, epslevel] = funcheb2.constructor(f, [], pref);
x = funcheb2.chebpts(length(values));
pass(11) = norm(f(x) - values, inf) < tol;
pass(12) = norm(vscale - [sin(1) cos(0) exp(1)], inf) < tol && logical(epslevel);

pref.funcheb2.extrapolate = 0;
pref.funcheb2.refinementFunction = 'resampling';
[values, coeffs, vscale, epslevel] = funcheb2.constructor(f, [], pref);
x = funcheb2.chebpts(length(values));
pass(13) = norm(f(x) - values, inf) < tol;
pass(14) = all(vscale == [sin(1) cos(0) exp(1)]) && logical(epslevel);

pref.funcheb2.extrapolate = 1;
pref.funcheb2.refinementFunction = 'resampling';
[values, coeffs, vscale, epslevel] = funcheb2.constructor(f, [], pref);
x = funcheb2.chebpts(length(values));
pass(15) = norm(f(x) - values, inf) < tol;
pass(16) = norm(vscale - [sin(1) cos(0) exp(1)], inf) < tol && logical(epslevel);

%%
% Some other tests:

% This should fail with an error:
try
    f = @(x) x + NaN;
    funcheb2.constructor(f);
    pass(17) = false;
catch ME
    pass(17) = strcmp(ME.message, 'Too many NaNs to handle.');
end

% As should this:
try
    f = @(x) x + Inf;
    funcheb2.constructor(f);
    pass(18) = false;
catch ME
    pass(18) = strcmp(ME.message, 'Too many NaNs to handle.');
end

% Test that the extrapolation option avoids endpoint evaluations.
pref.funcheb2.extrapolate = 1;
funcheb2.constructor(@(x) [F(x) F(x)], [], pref);
pass(19) = true;

    function y = F(x)
        if ( any(abs(x) == 1) )
            error('Extrapolate should prevent endpoint evaluation.');
        end
        y = sin(x);
    end

end