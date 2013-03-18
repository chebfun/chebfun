% Test file for funcheb/mldivide.

function pass = test_mldivide(pref)

% Get preferences.
if (nargin < 1)
    pref = funcheb.pref;
end

% Set a tolerance.  (pref.eps does not matter here.)
tol = 10*eps;

for ( n = 1:2 )
    if ( n == 1 )
        testclass = funcheb1();
    else 
        testclass = funcheb2();
    end

    %%
    % Basic correctness checks.

    % We get a known exact solution in this case.
    f = testclass.make(@(x) sin(x), [], [], pref);
    x = f \ f;
    err = f - x*f;
    pass(n, 1) = abs(x - 1) < tol;
    pass(n, 2) = max(abs(err.values(:))) < tol;

    % Same here.
    f = testclass.make(@(x) [sin(x) cos(x)], [], [], pref);
    g = testclass.make(@(x) sin(x + pi/4), [], [], pref);
    x = f \ g;
    err = g - f*x;
    pass(n, 3) = max(abs(x - [1/sqrt(2) ; 1/sqrt(2)])) < tol;
    pass(n, 4) = max(abs(err.values(:))) < tol;

    % A known least-squares solution.
    f = testclass.make(@(x) [ones(size(x)) x x.^2 x.^3], [], [], pref);
    g = testclass.make(@(x) x.^4 + x.^3 + x + 1, [], [], pref);
    x = f \ g;
    pass(n, 5) = max(abs(x - [32/35 ; 1 ; 6/7 ; 1])) < tol;

    %%
    % Check error conditions.

    % mldivide doesn't work between a FUNCHEB and a non-FUNCHEB.
    try
        f = testclass.make(@(x) [sin(x) cos(x) exp(x)], [], [], pref);
        f \ 2;
        pass(n, 6) = 0;
    catch ME
        pass(n, 6) = strcmp(ME.identifier, ...
            'CHEBFUN:FUNCHEB:mldivide:funchebMldivideUnknown');
    end
end

end

