% Test file for @chebfun/subsref.m.

function pass = test_subsref(pref)

% [TODO]: test calls of the form f(:,1) and f(1,:), etc.

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

% Test () syntaxes.
f = chebfun(@(x) sin(x - 0.1));
pass(1) = isequal(feval(f, xr), f(xr));
pass(2) = isequal(feval(f, xr, 'left'), f(xr, 'left'));
pass(3) = isequal(feval(f, xr, 'right'), f(xr, 'right'));
g = chebfun(@(x) cos(x + 0.2));
pass(4) = isequal(compose(g, f), f(g));
pass(5) = isequal(f, f(:));

try
    y = f('X');
    pass(6) = false;
catch ME
    pass(6) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:subsref:nonnumeric');
end

x_mtx = reshape(xr, [100 10]);
pass(7) = isequal(feval(f, x_mtx), f(x_mtx));
x_3mtx = reshape(xr, [5 20 10]);
pass(8) = isequal(feval(f, x_3mtx), f(x_3mtx));
x_4mtx = reshape(xr, [5 4 5 10]);
pass(9) = isequal(feval(f, x_4mtx), f(x_4mtx));

f.isTransposed = 1;
pass(10) = isequal(feval(f, xr), f(xr)) && isequal(feval(f, x_4mtx), f(x_4mtx));

% Test () syntaxes with an array-valued chebfun.
f = chebfun(@(x) [sin(x - 0.1) cos(x - 0.2)]);
pass(11) = isequal(feval(f, xr), f(xr));
y1 = feval(f, xr);
y1 = y1(:, 2);
y2 = f(xr, 2);
pass(12) = isequal(y1, y2);

x_mtx = reshape(xr, [100 10]);
pass(13) = isequal(feval(f, x_mtx), f(x_mtx));
x_3mtx = reshape(xr, [5 20 10]);
pass(14) = isequal(feval(f, x_3mtx), f(x_3mtx));
x_4mtx = reshape(xr, [5 4 5 10]);
pass(15) = isequal(feval(f, x_4mtx), f(x_4mtx));

f.isTransposed = 1;
pass(16) = isequal(feval(f, xr), f(xr)) && isequal(feval(f, x_4mtx), f(x_4mtx));

% Test {} syntaxes.
f = chebfun(@(x) sin(x - 0.1));
pass(17) = isequal(f, f{:});
err = norm(f{-1, -0.1, 0.2, 1} - restrict(f, [-1 -0.1 0.2 1]));
tol = epslevel(f);
pass(18) = err < tol;

try
    y = f{'X'};
    pass(19) = false;
catch ME
    pass(19) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:subsref:badDomain');
end

try
    index.subs = {[1 2], [3 4]}.';
    index.type = '{}';
    y = subsref(f, index);
    pass(20) = false;
catch ME
    pass(20) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:subsref:dimensions');
end

try
    index.subs = [];
    index.type = '[]';
    y = subsref(f, index);
    pass(21) = false;
catch ME
    pass(21) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:subsref:unexpectedType');
end

% Test {} syntaxes with an array-valued chebfun.
f = chebfun(@(x) [sin(x - 0.1) cos(x - 0.2)]);
err = norm(f{-1, -0.1, 0.2, 1} - restrict(f, [-1 -0.1 0.2 1]), inf);
tol = epslevel(f);
pass(22) = err < 5*tol;

end
