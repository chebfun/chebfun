% Test file for @chebfun/subsref.m.

function pass = test_feval(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfun.pref;
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
    pass(6) = strcmp(ME.identifier, 'CHEBFUN:subsref:nonnumeric');
end

f.isTransposed = 1;
pass(7) = isequal(feval(f, xr), f(xr));

% Test () syntaxes with an array-valued chebfun.
f = chebfun(@(x) [sin(x - 0.1) cos(x - 0.2)]);
pass(8) = isequal(feval(f, xr), f(xr));
y1 = feval(f, xr);
y1 = y1(:, 2);
y2 = f(xr, 2);
pass(9) = isequal(y1, y2);

f.isTransposed = 1;
pass(10) = isequal(feval(f, xr), f(xr));

% Test {} syntaxes.
f = chebfun(@(x) sin(x - 0.1));
pass(11) = isequal(f, f{:});
pass(12) = isequal(f{-1, -0.1, 0.2, 1}, restrict(f, [-1 -0.1 0.2 1]));

try
    y = f{'X'};
    pass(13) = false;
catch ME
    pass(13) = strcmp(ME.identifier, 'CHEBFUN:subsref:baddomain');
end

try
    index.subs = {[1 2], [3 4]}.';
    index.type = '{}';
    y = subsref(f, index);
    pass(14) = false;
catch ME
    pass(14) = strcmp(ME.identifier, 'CHEBFUN:subsref:dimensions');
end

try
    index.subs = [];
    index.type = '[]';
    y = subsref(f, index);
    pass(15) = false;
catch ME
    pass(15) = strcmp(ME.identifier, 'CHEBFUN:subsref:unexpectedType');
end

% Test {} syntaxes with an array-valued chebfun.
f = chebfun(@(x) [sin(x - 0.1) cos(x - 0.2)]);
pass(16) = isequal(f{-1, -0.1, 0.2, 1}, restrict(f, [-1 -0.1 0.2 1]));

end
