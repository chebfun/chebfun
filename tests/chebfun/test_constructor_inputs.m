function pass = test_constructor_inputs(pref)

if ( nargin == 0 )
    pref = chebpref();
end

% [TODO]: This test needs to be updated to include more exotic input options.

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
f2 = chebfun(@(x) x, 'kind', '2nd');
f3 = chebfun(@(x) x, 'chebkind', 1);
pass(10) = isa(f1.funs{1}.onefun, 'chebtech1') && ...
    isa(f2.funs{1}.onefun, 'chebtech2') && isa(f3.funs{1}.onefun, 'chebtech1');

% Test construction from numeric string.
f = chebfun('1');
pass(11) = all(feval(f, linspace(-1, 1, 10)) == 1);

% Test 'trunc', flag.
f = chebfun(@abs, 'trunc', 10, 'splitting', 'on');
c = get(f, 'coeffs');
pass(12) = abs(-4/63/pi - c{1}(2)) < get(f, 'epslevel');

% Test construction from cells of strings:
f = chebfun({'x','x-1'}, [0 1 2]);
pass(13) = norm(feval(f, [.5, 1.5]) - .5) < get(f, 'epslevel');

% Test construction from a piecewise chebfun:
f = chebfun(chebfun({'x','x'}, [-1 0 1]));
x = [-.5, .5];
pass(14) = numel(f.funs) == 1  && norm(feval(f, x) - x) < get(f, 'epslevel');

end
