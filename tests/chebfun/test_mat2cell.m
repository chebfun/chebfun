% Test file for @chebfun/mat2cell.m.

function pass = test_mat2cell(pref)

% Obtain preferences.
if ( nargin == 0 )
    pref = chebpref();
end

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

% Check empty case.
C = mat2cell(chebfun());
pass(1) = isequal(size(C), [1 1]) && isempty(C{1});

% Check a simple example.
F = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1], pref);
f = chebfun(@sin, [-1 1], pref);
g = chebfun(@cos, [-1 1], pref);
h = chebfun(@exp, [-1 1], pref);
fg = chebfun(@(x) [sin(x) cos(x)], pref);

C = mat2cell(F);
pass(2) = isequal(size(C), [1 3]) && (size(C{1}, 2) == 1) && ...
    (size(C{2}, 2) == 1) && (size(C{3}, 2) == 1);
pass(3) = norm(feval(C{1}, xr) - feval(f, xr), inf) < ...
    10*vscale(C{1})*epslevel(C{1});
pass(4) = norm(feval(C{2}, xr) - feval(g, xr), inf) < ...
    10*vscale(C{2})*epslevel(C{2});
pass(5) = norm(feval(C{3}, xr) - feval(h, xr), inf) < ...
    10*vscale(C{3})*epslevel(C{3});

C = mat2cell(F, [2 1]);
pass(6) = isequal(size(C), [1 2]) && (size(C{1}, 2) == 2) && ...
    (size(C{2}, 2) == 1);
err = feval(C{1}, xr) - feval(fg, xr);
pass(7) = norm(err(:), inf) < 10*vscale(C{1})*epslevel(C{1}) && ...
    norm(feval(C{2}, xr) - feval(h, xr), inf) < 10*vscale(C{2})*epslevel(C{2});

C = mat2cell(F, 1, [2 1]);
pass(8) = isequal(size(C), [1 2]) && (size(C{1}, 2) == 2) && ...
    (size(C{2}, 2) == 1);
err = feval(C{1}, xr) - feval(fg, xr);
pass(9) = norm(err(:), inf) < 10*vscale(C{1})*epslevel(C{1}) && ...
    norm(feval(C{2}, xr) - feval(h, xr), inf) < 10*vscale(C{2})*epslevel(C{2});

% Check that sizes come out right for row chebfuns.
Ft = F.';
ft = f.';
gt = g.';
ht = h.';
fgt = fg.';

C = mat2cell(Ft);
pass(10) = isequal(size(C), [3 1]) && (size(C{1}, 1) == 1) && ...
    (size(C{2}, 1) == 1) && (size(C{3}, 1) == 1);

C = mat2cell(Ft, [2 1]);
pass(11) = isequal(size(C), [2 1]) && (size(C{1}, 1) == 2) && ...
    (size(C{2}, 1) == 1);

C = mat2cell(Ft, [2 1], 1);
pass(12) = isequal(size(C), [2 1]) && (size(C{1}, 1) == 2) && ...
    (size(C{2}, 1) == 1);

% Check error conditions.
try
    C = mat2cell(F, 'bad', [2 1]);
    pass(13) = false;
catch ME
    pass(13) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:mat2cell:size');
end

try
    C = mat2cell(F, 2, [2 1]);
    pass(14) = false;
catch ME
    pass(14) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:mat2cell:size');
end

try
    C = mat2cell(F, 1, [2 2]);
    pass(15) = false;
catch ME
    pass(15) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:mat2cell:size');
end

try
    C = mat2cell(Ft, [2 1], 2);
    pass(16) = false;
catch ME
    pass(16) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:mat2cell:size');
end

try
    C = mat2cell(Ft, [2 2], 1);
    pass(17) = false;
catch ME
    pass(17) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:mat2cell:size');
end

end
