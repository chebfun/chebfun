% Test file for funcheb2/mldivide.

function pass = test_mldivide(pref)

% Get preferences.
if (nargin < 1)
    pref = funcheb2.pref();
end

% Set a tolerance.  (pref.eps does not matter here.)
tol = 10*eps;

%%
% Basic correctness checks.

% We get a known exact solution in this case.
f = funcheb2(@(x) sin(x), pref);
x = f \ f;
err = f - x*f;
pass(1) = abs(x - 1) < tol;
pass(2) = max(abs(err.values(:))) < tol;

% Same here.
f = funcheb2(@(x) [sin(x) cos(x)], pref);
g = funcheb2(@(x) sin(x + pi/4), pref);
x = f \ g;
err = g - f*x;
pass(3) = max(abs(x - [1/sqrt(2) ; 1/sqrt(2)])) < tol;
pass(4) = max(abs(err.values(:))) < tol;

% A known least-squares solution.
f = funcheb2(@(x) [ones(size(x)) x x.^2 x.^3], pref);
g = funcheb2(@(x) x.^4 + x.^3 + x + 1, pref);
x = f \ g;
err = g - f*x;
pass(5) = max(abs(x - [32/35 ; 1 ; 6/7 ; 1])) < tol;

%%
% Check error conditions.

% mldivide doesn't work between a FUNCHEB2 and a non-FUNCHEB2.
try
    f = funcheb2(@(x) [sin(x) cos(x) exp(x)], pref);
    f \ 2;
    pass(6) = 0;
catch ME
    pass(6) = strcmp(ME.identifier, ...
        'CHEBFUN:FUNCHEB2:mldivide:funcheb2MldivideUnknown');
end

end

