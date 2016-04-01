function pass = test_constructor_inputs(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% [TODO]: This test needs to be updated to include more exotic input options.

% Generate a few random points to use as test values.
seedRNG(6178);
xx = 2 * rand(100, 1) - 1;

% Test the vectorise flag:
try
    f = chebfun(@(x) x^2, pref, 'vectorize'); %#ok<NASGU>
    pass(1) = true;
catch
    pass(1) = false;
end

% Test the vectorise flag:
f = chebfun(@(x) sin(x), pref, 5);
pass(2) = ( length(f.funs{1}) == 5 ); 
% [TODO]: Change this once @CHEBFUN/LENGTH is implemented.

% Test the 'splitting on' flag.
f = chebfun(@(x) abs(x), 'splitting', 'on');
pass(3) = all(size(f.funs) == [1, 2]);

% Test the 'extrapolate' flag.
f = chebfun(@(x) sign(x), [0, 1], 'extrapolate', 'on');
pass(4) = get(f, 'ishappy');

% Test construction from an array of FUN objects:
f = chebfun({0, 1, 2, 3, 4, 5}, [0, 1, 2, 3, 4, 5, 6]);
g = chebfun(f.funs);
pass(5) = all(size(g.funs) == [1, 6]) && all(g.domain == 0:6);

% Test array-valued construction from an array of FUN objects:
f = chebfun({[0, 1], [2, 3], [4, 5]}, [0, 1, 2, 3]);
g = chebfun(f.funs);
pass(6) = all(size(g.funs) == [1, 3]) && all(g.domain == 0:3);

% Test fixed-length construction:
f = chebfun(@sin, 10);
pass(7) = numel(f.funs) == 1 && size(f.funs{1},1) == 10;

% Test equispaced construction:
f = chebfun(rand(10,3), 'equi');
pass(8) = numel(f.funs) == 1 && size(f.funs{1}, 1) > 10;

% Test construction from coefficients.
f = chebfun([1 ; 2 ; 3], 'coeffs');
pass(9) = isequal(f.funs{1}.onefun.coeffs, [1 ; 2 ; 3]);

% Test 'chebkind' and 'kind' flags.
f1 = chebfun(@(x) x, 'chebkind', '1st');
f3 = chebfun(@(x) x, 'chebkind', 1);
pass(10) = isa(f1.funs{1}.onefun, 'chebtech1') && ...
    isa(f3.funs{1}.onefun, 'chebtech1');

% Test construction from numeric string.
f = chebfun('1');
pass(11) = all(feval(f, linspace(-1, 1, 10)) == 1);

% Test 'trunc', flag.
f = chebfun(@abs, 'trunc', 10, 'splitting', 'on');
c = get(f, 'coeffs');
pass(12) = abs(-4/63/pi - c(9)) < 10*eps;

% Test construction from cells of strings:
f = chebfun({'x','x-1'}, [0 1 2]);
pass(13) = norm(feval(f, [.5, 1.5]) - .5) < eps;

% Test construction from a chebfun:
x = chebfun('x', [0, 5]);
f = chebfun(x);
pass(13) = all(domain(f) == [0 5] );

f = chebfun(x, [0 2]);
pass(14) = all(domain(f) == [0 2] );

% Test construction from a piecewise chebfun:
f = chebfun(chebfun({'x','x'}, [-1 0 1]));
x = [-.5, .5];
pass(15) = numel(f.funs) == 1  && norm(feval(f, x) - x) < eps;

% Test 'minSamples' flag.
f_op = @(x) -x - x.^2 + exp(-(50*(x - .5)).^4);
f1 = chebfun(f_op, 'minSamples', 17);
err1 = norm(feval(f1, xx) - f_op(xx), inf);
f2 = chebfun(f_op, 'minSamples', 33);
err2 = norm(feval(f2, xx) - f_op(xx), inf);
pass(16) = (err1 > 1e-3) && (err2 < 1e2*vscale(f2)*eps);

% Test support for "legacy" preferences.
f_op = @(x) sin(200*x);

f = chebfun(f_op, 'resampling', 'on');
pass(17) = ishappy(f);

warnstate = warning('off','CHEBFUN:CHEBFUN:constructor:notResolved');
f = chebfun(f_op, 'maxdegree', 129, 'tech', @chebtech2);
pass(18) = ~ishappy(f) && (length(f) == 129);
warning(warnstate);

f = chebfun(f_op, 'splitting', 'on', 'splitdegree', 65);
pass(19) = ishappy(f) && all(cellfun(@(fk) length(fk) <= 65, f.funs));

warnstate = warning('off','CHEBFUN:CHEBFUN:constructor:funNotResolved');
f = chebfun(f_op, 'splitting', 'on', 'splitLength', 65, 'splitMaxLength', 200);
pass(20) = ~ishappy(f) && (length(f) <= 65*4);
warning(warnstate);

% Test construction with a mixture of preference object and keyword inputs.
p = pref;
p.tech = @chebtech1;
f = chebfun(@(x) 1./x, [0 1], 'exps', [-1 0], p);
pass(21) = get(f, 'ishappy') && isa(f.funs{1}.onefun.smoothPart, 'chebtech1');
f = chebfun(@(x) 1./x, [0 1], p, 'exps', [-1 0]);
pass(22) = get(f, 'ishappy') && isa(f.funs{1}.onefun.smoothPart, 'chebtech1');

p = struct();
p.tech = @chebtech1;
f = chebfun(@(x) 1./x, [0 1], 'exps', [-1 0], p);
pass(23) = get(f, 'ishappy') && isa(f.funs{1}.onefun.smoothPart, 'chebtech1');
f = chebfun(@(x) 1./x, [0 1], p, 'exps', [-1 0]);
pass(24) = get(f, 'ishappy') && isa(f.funs{1}.onefun.smoothPart, 'chebtech1');

% Test the vectorise flag on an array-valued function:
try
    f = chebfun(@(x) [x^2 sin(x)*cos(x)], pref, 'vectorize'); %#ok<NASGU>
    pass(25) = true;
catch
    pass(25) = false;
end

% Test the 'eps' preference.
f = chebfun(@sin);
g = chebfun(@sin, 'eps', 1e-6);
pass(26) = length(g) < length(f);

end
