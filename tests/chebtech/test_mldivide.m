% Test file for chebtech/mldivide.m

function pass = test_mldivide(pref)

% Get preferences.
if (nargin < 1)
    pref = chebtech.techPref();
end

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    %%
    % Basic correctness checks.

    % We get a known exact solution in this case.
    f = testclass.make(@(x) sin(x), [], pref);
    x = f \ f;
    err = f - x*f;
    pass(n, 1) = abs(x - 1) < 10*f.vscale.*f.epslevel;
    pass(n, 2) = max(abs(err.coeffs(:))) < 10*f.vscale.*f.epslevel;

    % Same here.
    f = testclass.make(@(x) [sin(x) cos(x)], [], pref);
    g = testclass.make(@(x) sin(x + pi/4), [], pref);
    tol_f = 10*max(f.vscale.*f.epslevel);
    tol_g = 10*max(f.vscale.*f.epslevel);
    x = f \ g;
    err = g - f*x;
    pass(n, 3) = max(abs(x - [1/sqrt(2) ; 1/sqrt(2)])) < max(tol_f, tol_g);
    pass(n, 4) = max(abs(err.coeffs(:))) < max(tol_f, tol_g);

    % A known least-squares solution.
    f = testclass.make(@(x) [ones(size(x)) x x.^2 x.^3], [], pref);
    g = testclass.make(@(x) x.^4 + x.^3 + x + 1, [], pref);
    tol_f = 10*max(f.vscale.*f.epslevel);
    tol_g = 10*max(f.vscale.*f.epslevel);
    x = f \ g;
    pass(n, 5) = max(abs(x - [32/35 ; 1 ; 6/7 ; 1])) < max(tol_f, tol_g);

    %%
    % Check error conditions.

    % mldivide doesn't work between a CHEBTECH and a non-CHEBTECH.
    try
        f = testclass.make(@(x) [sin(x) cos(x) exp(x)], [], pref);
        f \ 2; %#ok<VUNUS>
        pass(n, 6) = 0;
    catch ME
        pass(n, 6) = strcmp(ME.identifier, ...
            'CHEBFUN:CHEBTECH:mldivide:chebtechMldivideUnknown');
    end
end

end

