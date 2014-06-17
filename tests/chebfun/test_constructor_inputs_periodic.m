function pass = test_constructor_inputs_periodic(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
xx = 2 * rand(100, 1) - 1;

% Generate a few random points to use as test values.
seedRNG(6178);
xx = 2 * rand(100, 1) - 1;

% Test the vectorise flag:
try
    f = chebfun(@(x) cos(pi*x)^2, pref, 'periodic', 'vectorize'); %#ok<NASGU>
    pass(1) = true;
catch
    pass(1) = false;
end

% Test that construction over a non-smooth domain fails:
try
    f = chebfun(@(x) cos(pi*x), [-1 0 1], 'periodic');
    pass(2) = false;
catch ME
    pass(2) = strcmp(ME.identifier,'CHEBFUN:parseInputs:periodic');
end

% Test fixed-length construction:
f = chebfun(@(x) exp(sin(pi*x)), 'periodic', 10);
pass(3) = numel(f.funs) == 1 && size(f.funs{1},1) == 10;

% Test equispaced construction:
f = chebfun(rand(10,3), 'periodic', 'equi');
pass(4) = numel(f.funs) == 1 && size(f.funs{1}, 1) == 10;

% Test construction from coefficients.
f = chebfun([0.5 ; 0 ; 0.5], 'periodic', 'coeffs');
pass(5) = isequal(f.funs{1}.onefun.coeffs, [0.5 ; 0 ; 0.5]);

% Test construction from numeric string.
f = chebfun('1', 'periodic');
pass(6) = all(feval(f, linspace(-1, 1, 10)) == 1);

% Test 'trunc', flag.
f = chebfun(@(x) 1 + sin(pi*sin(10*pi*x)), 'periodic', 'trunc', 11);
c = get(f, 'coeffs');
pass(7) = abs(1 - c(6)) < get(f, 'epslevel');

end
